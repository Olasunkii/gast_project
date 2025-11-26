configfile: "config.yaml"

# prepare safe organism name & config paths variables setup
organism_safe = config['organism'].replace(' ', '_')
METADATA_DIR = config['paths']['metadata_dir']
SEQ_DIR = config['paths']['sequences_dir']
AST_DIR = config['paths']['ast_dir']
RESULTS_DIR = config['paths']['results_dir']

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
        expand(f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_1.fastq.gz", sample=get_sample_ids),
        expand(f"{RESULTS_DIR}/fastp/{{sample}}/{{sample}}_2.fastq.gz", sample=get_sample_ids),

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
# Rule: fastp — cleaning sequences
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