configfile: "config.yaml"

# prepare safe organism name & config paths variables setup
organism_safe = config['organism'].replace(' ', '_')
METADATA_DIR = config['paths']['metadata_dir']
SEQ_DIR = config['paths']['sequences_dir']
AST_DIR = config['paths']['ast_dir']
RESULTS_DIR = config['paths']['results_dir']
db_DIR = config['paths']['reference_db_dir']

# -------------------------------------------------------
# Function: load sample names dynamically after metadata exists
# -------------------------------------------------------
import os

def get_sample_ids(wildcards):
    path = checkpoints.download_sequences.get(**wildcards).output[0]
    if not os.path.exists(path):
        return []
    with open(path) as f:
        header = f.readline().strip().split(',')
        idx = header.index('Run')
        return [line.split(',')[idx].strip() for line in f if line.strip()]

# -------------------------------------------------------
# Rule: all — desired output including downstream targets
# -------------------------------------------------------
rule all:
    input:
        f"{METADATA_DIR}/SraRunInfo_{organism_safe}.csv",
        f"{SEQ_DIR}/download_log.csv",
        expand(f"{RESULTS_DIR}/checkm/{{sample}}", sample=get_sample_ids),
        expand(f"{RESULTS_DIR}/bakta/{{sample}}/{{sample}}.fna", sample=get_sample_ids),
        expand(f"{RESULTS_DIR}/amrfinder/{{sample}}/{{sample}}_amrfinder.tsv", sample=get_sample_ids)

# -------------------------------------------------------
# Rule: fetch_metadata — downloads SraRunInfo file
# -------------------------------------------------------
rule fetch_metadata:
    output:
        out=f"{METADATA_DIR}/SraRunInfo_{organism_safe}.csv"
    params:
        email=config["email"],
        organism=organism_safe,
        retmax=config["retmax"],
        outdir=METADATA_DIR
    shell:
        """
        python3 src/sra_extractor_metadata.py {params.email} {params.organism} {params.retmax}
        """

# -------------------------------------------------------
# Rule: download_sequences — downloads and logs FASTQ data
# -------------------------------------------------------
checkpoint download_sequences:
    input:
        metadata = f"{METADATA_DIR}/SraRunInfo_{organism_safe}.csv"
    output:
        touch(f"{SEQ_DIR}/download_log.csv")
    shell:
        """
        python3 src/sra_extractor.py \
            --email {config[email]} \
            --organism {organism_safe} \
            --retmax {config[retmax]} \
            --download
        """

# -------------------------------------------------------
# Rule: fastp — cleaning sequences & sequence analysis
# -------------------------------------------------------
rule fastp:
    input:
        r1 = f"{SEQ_DIR}/{{sample}}/{{sample}}_1.fastq.gz",
        r2 = f"{SEQ_DIR}/{{sample}}/{{sample}}_2.fastq.gz"
    output:
        r1 = f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_1.fastq.gz",
        r2 = f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_2.fastq.gz",
        html = f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_fastp.html",
        json = f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_fastp.json"
    conda:
        "envs/environment_qc.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -h {output.html} -j {output.json}
        """
# -------------------------------------------------------
# Rule: trimgalore — adapter trimming
# -------------------------------------------------------
rule trim_galore:
    input:
        r1=f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_2.fastq.gz"
    output:
        r1=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz",
        r2=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz"
    conda:
        "envs/environment_qc.yaml"
    shell:
        """
        trim_galore --paired --gzip -o {RESULTS_DIR}/trim_galore/{wildcards.sample} {input.r1} {input.r2}
        """
# -------------------------------------------------------
# Rule: fastqc — quality of sequences
# -------------------------------------------------------
rule fastqc:
    input:
        R1 = f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz",
        R2 = f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz"
    output:
        html1 = f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_1_val_1_fastqc.html",
        zip1  = f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_1_val_1_fastqc.zip",
        html2 = f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_2_val_2_fastqc.html",
        zip2  = f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_2_val_2_fastqc.zip"
    conda:
        "envs/environment_qc.yaml"
    shell:
        "fastqc {input.R1} {input.R2} -o {RESULTS_DIR}/fastqc/{wildcards.sample}"
# -------------------------------------------------------
# Rule: unicycler — assembly to draft genome
# -------------------------------------------------------
rule unicycler_assembly:
    input:
        r1=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz",
        r2=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz"
    output:
        f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta"
    conda:
        "envs/environment_qc.yaml"
    threads: 16
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -o {RESULTS_DIR}/assembly/{wildcards.sample} -t {threads}"
# -------------------------------------------------------
# Rule: Checkm2 — Analyze contamination on draft genome
# -------------------------------------------------------
rule checkm:
    input:
        f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta"
    output:
        directory(f"{RESULTS_DIR}/checkm/{{sample}}")
    conda:
        "envs/environment_checkm.yaml"
    threads: 16
    params:
        db=f"{db_DIR}/checkm2_db/uniref100.KO.1.dmnd"
    shell:
        """
        checkm2 predict --input {input} --output-directory {output} \
                --threads {threads} --database_path {params.db}
        """
# -------------------------------------------------------
# Rule: Bakta — gene annotation
# -------------------------------------------------------
rule bakta:
    input:
        f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta"
    output:
        fna = f"{RESULTS_DIR}/bakta/{{sample}}/{{sample}}.fna"
    conda:
        "envs/environment_bakta.yaml"
    threads: 16
    params:
        db = f"{db_DIR}/bakta_db"
    shell:
        """
        rm -rf results/bakta/{wildcards.sample}
        bakta --db {params.db} \
              --output {RESULTS_DIR}/bakta/{wildcards.sample} \
              --prefix {wildcards.sample} \
              --threads {threads} \
              {input}
        """
# -------------------------------------------------------
# Rule: AMRFinderPlus — Antibiotic gene resistance annotation
# -------------------------------------------------------
rule amrfinder_db:
    output:
        touch("amrfinder_db_ready.txt")
    conda:
        "envs/environment_amr.yaml"
    shell:
        """
        amrfinder --force_update
        """

rule amrfinder:
    input:
        fasta = f"{RESULTS_DIR}/bakta/{{sample}}/{{sample}}.fna",
        db_updated = "amrfinder_db_ready.txt"
    output:
        out = f"{RESULTS_DIR}/amrfinder/{{sample}}/{{sample}}_amrfinder.tsv"
    conda:
        "envs/environment_amr.yaml"
    threads: 8
    shell:
        """
        mkdir -p {RESULTS_DIR}/amrfinder/{wildcards.sample}
        amrfinder --nucleotide {input.fasta} --annotation_format bakta \
                  --plus \
                  --organism Klebsiella_pneumoniae \
                  --output {output.out} \
                  --threads {threads}
        """
