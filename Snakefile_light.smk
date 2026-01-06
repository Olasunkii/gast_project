import os
configfile: "configs/config_parameter.yaml"
configfile: "configs/config.yaml"

# Variables setup
ORGANISM_SAFE = config['organism'].replace(' ', '_')
METADATA_DIR = config['paths']['metadata_dir']
HOST_METADATA_DIR = config['paths']['host_metadata_dir']
SEQ_DIR = config['paths']['sequences_dir']
AST_DIR = config['paths']['ast_dir']
RESULTS_DIR = config['paths']['results_dir']
DB_DIR = config['paths']['reference_db_dir']

# -------------------------------------------------------
# Function: load sample id's dynamically after metadata exists
# -------------------------------------------------------
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
        # only keep Metadata
        f"{METADATA_DIR}/SraRunInfo_{ORGANISM_SAFE}.csv",
        # Keep Consistency Checks
        f"{RESULTS_DIR}/consistency_checks/carbapenem_consistency_check.tsv",
        f"{RESULTS_DIR}/consistency_checks/genome_consistency_check.tsv",
        # only keep Final ML Outputs
        f"{RESULTS_DIR}/X_train.csv",
        f"{RESULTS_DIR}/y_train.csv",
        f"{RESULTS_DIR}/X_validation.csv",
        f"{RESULTS_DIR}/y_validation.csv",
        f"{RESULTS_DIR}/X_test.csv",
        f"{RESULTS_DIR}/y_test.csv"
# -------------------------------------------------------
# Rule: fetch_metadata — downloads SraRunInfo file
# -------------------------------------------------------       
rule fetch_metadata:
    output:
        out=f"{METADATA_DIR}/SraRunInfo_{ORGANISM_SAFE}.csv",
        host_metadata = f"{HOST_METADATA_DIR}/host_metadata_all.csv",
        ast_dir = directory(AST_DIR)
    params:
        email=config["email"],
        organism=ORGANISM_SAFE,
        retmax=config["retmax"]
    shell:
        "python3 src/sra_extractor_metadata.py {params.email} {params.organism} {params.retmax}"
# -------------------------------------------------------
# Rule: download_sequences — downloads and logs FASTQ data
# -------------------------------------------------------
checkpoint download_sequences:
    input:
        metadata = f"{METADATA_DIR}/SraRunInfo_{ORGANISM_SAFE}.csv"
    output:
        temp(touch(f"{SEQ_DIR}/download_log.csv")) # Temp: Delete raw downloads after processing
    shell:
        "python3 src/sra_extractor.py --email {config[email]} --organism {ORGANISM_SAFE} --retmax {config[retmax]} --download"
# -------------------------------------------------------
# Rule: trimgalore — adapter trimming
# -------------------------------------------------------
rule trim_galore:
    input:
        r1 = f"{SEQ_DIR}/{{sample}}/{{sample}}_1.fastq.gz",
        r2 = f"{SEQ_DIR}/{{sample}}/{{sample}}_2.fastq.gz"
    output:
        r1=temp(f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz"),
        r2=temp(f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz")
    conda: "envs/environment_qc.yaml"
    shell: "trim_galore --paired --gzip -o {RESULTS_DIR}/trim_galore/{wildcards.sample} {input.r1} {input.r2}"
# -------------------------------------------------------
# Rule: fastqc — quality of sequences
# -------------------------------------------------------
rule fastqc:
    input:
        R1 = f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz",
        R2 = f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz"
    output:
        html1 = temp(f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_1_val_1_fastqc.html"),
        zip1  = temp(f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_1_val_1_fastqc.zip"),
        html2 = temp(f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_2_val_2_fastqc.html"),
        zip2  = temp(f"{RESULTS_DIR}/fastqc/{{sample}}/{{sample}}_2_val_2_fastqc.zip")
    conda: "envs/environment_qc.yaml"
    shell: "fastqc {input.R1} {input.R2} -o {RESULTS_DIR}/fastqc/{wildcards.sample}"
# -------------------------------------------------------
# Rule: unicycler — assembly to draft genome
# -------------------------------------------------------
rule unicycler_assembly:
    input:
        r1=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_1_val_1.fq.gz",
        r2=f"{RESULTS_DIR}/trim_galore/{{sample}}/{{sample}}_2_val_2.fq.gz"
    output:
        temp(f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta")
    conda: "envs/environment_qc.yaml"
    threads: 16
    shell: "unicycler -1 {input.r1} -2 {input.r2} -o {RESULTS_DIR}/assembly/{wildcards.sample} -t {threads}"
# -------------------------------------------------------
# Rule: Checkm2 — Analyze contamination on draft genome
# -------------------------------------------------------
rule checkm:
    input: f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta"
    output: temp(directory(f"{RESULTS_DIR}/checkm/{{sample}}"))
    conda: "envs/environment_checkm.yaml"
    threads: 16
    params: db=f"{DB_DIR}/checkm2_db/uniref100.KO.1.dmnd"
    shell: "checkm2 predict --input {input} --output-directory {output} --threads {threads} --database_path {params.db}"
# -------------------------------------------------------
# Rule: Bakta — gene annotation
# -------------------------------------------------------
rule bakta:
    input: f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta"
    output: fna = temp(f"{RESULTS_DIR}/bakta/{{sample}}/{{sample}}.fna")
    conda: "envs/environment_bakta.yaml"
    threads: 16
    params: db = f"{DB_DIR}/bakta_db"
    shell: "bakta --db {params.db} --output {RESULTS_DIR}/bakta/{wildcards.sample} --prefix {wildcards.sample} --threads {threads} {input}"
# -------------------------------------------------------
# Rule: AMRFinderPlus — Antibiotic gene resistance annotation
# -------------------------------------------------------
rule amrfinder_db:
    output: temp(touch("amrfinder_db_ready.txt"))
    conda: "envs/environment_amr.yaml"
    shell: "amrfinder --force_update"

rule amrfinder:
    input:
        fasta = f"{RESULTS_DIR}/bakta/{{sample}}/{{sample}}.fna",
        db_updated = "amrfinder_db_ready.txt"
    output: out = temp(f"{RESULTS_DIR}/amrfinder/{{sample}}.tsv")
    conda: "envs/environment_amr.yaml"
    threads: 8
    shell: "amrfinder --nucleotide {input.fasta} --annotation_format bakta --plus --organism Klebsiella_pneumoniae --output {output.out} --threads {threads}"
# -------------------------------------------------------
# Rule: Feature Extraction — transforming resistance genes to gene presence/absence
# -------------------------------------------------------
rule amrfinder_transformation:
    input: expand(f"{RESULTS_DIR}/amrfinder/{{sample}}.tsv", sample=get_sample_ids)
    output: temp(f"{RESULTS_DIR}/amrfinder/amr_transformed.tsv")
    conda: "envs/environment_python.yaml"
    shell: "python src/amr_transformer.py --amr_files {input} --output {output}"

#-------------------------------------------
# Rule: Consistency check — phenotypic data vs EUCAST breaking points
# -------------------------------------------------------
rule check_carbapenems:
    input: AST_DIR
    output: f"{RESULTS_DIR}/consistency_checks/carbapenem_consistency_check.tsv"
    conda: "envs/environment_python.yaml"
    shell: "python src/phenotype_checker.py --input {input} --config configs/config.yaml --output {output}"
# -------------------------------------------------------
# Rule: Consistency check — Genome completion
# -------------------------------------------------------
rule check_genome:
    input: expand(f"{RESULTS_DIR}/checkm/{{sample}}", sample=get_sample_ids)
    output: f"{RESULTS_DIR}/consistency_checks/genome_consistency_check.tsv"
    conda: "envs/environment_python.yaml"
    shell: "python src/genome_checker.py --input '{RESULTS_DIR}/checkm/' --config configs/config_parameter.yaml --output {output}"
# -------------------------------------------------------
# Rule: genotypic and phenotypic integration
# -------------------------------------------------------
rule integration:
    input:
        metadata_file = f"{METADATA_DIR}/SraRunInfo_{ORGANISM_SAFE}.csv",
        host_metadata= f"{HOST_METADATA_DIR}/host_metadata_all.csv",
        assembly_dir = expand(f"{RESULTS_DIR}/assembly/{{sample}}/assembly.fasta",sample=get_sample_ids),
        amr_file= f"{RESULTS_DIR}/amrfinder/amr_transformed.tsv"
    output: temp(f"{RESULTS_DIR}/integrated_data/integrated_data.csv")
    conda: "envs/environment_python.yaml"
    shell: "python src/DataIntegrator.py {SEQ_DIR} {input.metadata_file} {input.host_metadata} {AST_DIR} '{input.assembly_dir}' {input.amr_file} {output}"
# -------------------------------------------------------
# Rule: ML preprocessing integrated data
# -------------------------------------------------------
rule preprocessing:
    input: data=f"{RESULTS_DIR}/integrated_data/integrated_data.csv", config="configs/config.yaml"
    output: temp(f"{RESULTS_DIR}/integrated_data/integrated_data_preprocessed.csv")
    conda: "envs/environment_python.yaml"
    shell: "python src/ml_preprocessor.py --input {input.data} --config {input.config} --output {output}"

# -------------------------------------------------------
# Rule: ML prepration - target creation & formatting
# -------------------------------------------------------
rule run_ml_builder:
    input:
        data=f"{RESULTS_DIR}/integrated_data/integrated_data_preprocessed.csv",
        config_main="configs/config.yaml",
        config_params="configs/config_parameter.yaml"
    output:
        X_train=f"{RESULTS_DIR}/X_train.csv",
        y_train=f"{RESULTS_DIR}/y_train.csv",
        X_validation=f"{RESULTS_DIR}/X_validation.csv",
        y_validation=f"{RESULTS_DIR}/y_validation.csv",
        X_test=f"{RESULTS_DIR}/X_test.csv",
        y_test=f"{RESULTS_DIR}/y_test.csv"
    params: outdir=RESULTS_DIR
    conda: "envs/environment_python.yaml"
    shell: "python src/MLBuilder.py --input {input.data} --config {input.config_main} --config-params {input.config_params} --output {params.outdir}"