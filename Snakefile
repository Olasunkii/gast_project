# Snakefile
configfile: "config.yaml"

# -------------------------------------------------------
# Rule: all — defines the final expected output(s)
# -------------------------------------------------------
rule all:
    input:
        f"{config['paths']['metadata_dir']}/SraRunInfo_{config['organism'].replace(' ', '_')}.csv",
        f"{config['paths']['sequences_dir']}/download_log.csv"

# -------------------------------------------------------
# Rule: download_sequences — downloads FASTQ files
# -------------------------------------------------------
rule download_sequences:
    input:
        metadata = f"{config['paths']['metadata_dir']}/SraRunInfo_{config['organism'].replace(' ', '_')}.csv"
    output:
        touch(f"{config['paths']['sequences_dir']}/download_log.csv")
    container:
        "wgs_pipeline:latest"
    shell:
         """
        mkdir -p {config[paths][sequences_dir]}
        python3 /app/scripts/sra_extractor.py \
            --email {config[email]} \
            --organism "{config[organism]}" \
            --retmax {config[retmax]} \
            --download
        """
