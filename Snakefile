# Snakefile
configfile: "config.yaml"
organism_safe = config['organism'].replace(' ', '_')

# -------------------------------------------------------
# Rule: all — defines the final expected output(s)
# -------------------------------------------------------
rule all:
    input:
        f"{config['paths']['metadata_dir']}/SraRunInfo_{config['organism'].replace(' ', '_')}.csv",
        f"{config['paths']['sequences_dir']}/download_log.csv"

# -------------------------------------------------------
# Rule: downloading metadata — defines the input for fastq files
# -------------------------------------------------------
rule fetch_metadata:
    output:
        out=f"{config['paths']['metadata_dir']}/SraRunInfo_{config['organism'].replace(' ', '_')}.csv"
    params:
        email=config["email"],
        organism=organism_safe,
        retmax=config["retmax"],
        outdir=config["paths"]["metadata_dir"]
    container:
        "wgs_pipeline:latest"
    shell:
        """
        mkdir -p {params.outdir}
        python3 scripts/sra_extractor_metadata.py {params.email} {params.organism} {params.retmax}
        # mv SraRunInfo_{params.organism}.csv {params.outdir}/
        """
# -------------------------------------------------------
# Rule: download_sequences — downloads FASTQ files
# -------------------------------------------------------
rule download_sequences:
    input:
        metadata = f"{config['paths']['metadata_dir']}/SraRunInfo_{organism_safe}.csv"
    output:
        touch(f"{config['paths']['sequences_dir']}/download_log.csv")
    container:
        "wgs_pipeline:latest"
    shell:
         """
        mkdir -p {config[paths][sequences_dir]}
        python3 scripts/sra_extractor.py \
            --email {config[email]} \
            --organism {organism_safe} \
            --retmax {config[retmax]} \
            --download
        """
