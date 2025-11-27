# ===============================
# Snakefile (fixed + automated)
# ===============================

configfile: "config.yaml"

metadata_dir = config["paths"]["metadata_dir"]
seq_dir      = config["paths"]["sequences_dir"]
organism_tag = config["organism"].replace(" ", "_")


# -------------------------------------------------------
# FINAL OUTPUTS
# -------------------------------------------------------
rule all:
    input:
        f"{metadata_dir}/SraRunInfo_{organism_tag}.csv",
        f"{seq_dir}/download_log.csv"


# -------------------------------------------------------
# RULE 1 — AUTO-GENERATE METADATA
# -------------------------------------------------------
rule fetch_metadata:
    output:
        f"{metadata_dir}/SraRunInfo_{organism_tag}.csv"
    container:
        "wgs_pipeline:latest"
    shell:
        """
        mkdir -p {metadata_dir}
        python3 /app/scripts/sra_extractor.py \
            --email {config[email]} \
            --organism "{config[organism]}" \
            --retmax {config[retmax]}
        """


# -------------------------------------------------------
# RULE 2 — DOWNLOAD ALL FASTQ FILES
# -------------------------------------------------------
rule download_sequences:
    input:
        metadata = f"{metadata_dir}/SraRunInfo_{organism_tag}.csv"
    output:
        touch(f"{seq_dir}/download_log.csv")
    container:
        "wgs_pipeline:latest"
    shell:
        """
        mkdir -p {seq_dir}
        python3 /app/scripts/sra_extractor.py \
            --email {config[email]} \
            --organism "{config[organism]}" \
            --retmax {config[retmax]} \
            --download
        """