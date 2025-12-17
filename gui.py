import streamlit as st
import yaml
import subprocess

st.title("Snakemake Launcher")


#---GUI inputs (matched to the config_parameter.yaml)---
organism = st.selectbox(
    "Organism",
    [
        "Klebsiella pneumoniae"
    ]
)

email = st.text_input("Email (needed to retrieve data from NCBI)", value="example@mail.com")
retmax = st.number_input("NCBI retmax", min_value=1, max_value=100000, value=1)
cores = st.number_input(
    "Number of CPU cores",
    min_value=1,
    max_value=256,
    value=16
)
latency_wait = st.number_input(
    "Filesystem latency wait (seconds)",
    min_value=0,
    max_value=600,
    value=30
)
# ---- ML section ----
st.subheader("ML configuration")

column_target_name = st.text_input(
    "Target column name",
    value="carbapenem_resistant"
)

scaling_method = st.selectbox(
    "Scaling method",
    ["none", "standard", "robust"]
)

test_size = st.number_input("Test size", min_value=0.05, max_value=0.9, value=0.3)
validation_size = st.number_input("Validation size", min_value=0.05, max_value=0.5, value=0.1)
stratify = st.checkbox("Stratify split", value=True)
random_state = st.number_input("Random state", min_value=0, value=42)

#Settings for Genome consistency check
st.subheader("Genome consistency QC")

min_completeness = st.number_input("Min completeness (%)", value=99)
max_contamination = st.number_input("Max contamination (%)", value=3)
max_contigs = st.number_input("Max contigs", value=500)


#---Write config and run snakemake---
if st.button("Run workflow"):
    config = {
        "organism": organism,
        "email": email,
        "retmax": retmax,
        "ml": {
            "column_target_name": column_target_name,
            "scaling": {
                "method": scaling_method
            },
            "split": {
                "test_size": test_size,
                "validation_size": validation_size,
                "stratify": stratify,
                "random_state": random_state,
            },
        },
        "genome_completeness_qc": {
            "completeness": {
                "min": min_completeness
            },
            "contamination": {
                "max": max_contamination
            },
            "total_contigs": {
                "max": max_contigs
            },
        },
    }

    with open("configs/config_parameter.yaml", "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    #build the snakemake command
    cmd = [
        "snakemake",
        "--snakefile", "Snakefile",
        "--cores", str(cores),
        "--latency-wait", str(latency_wait),
        "--use-conda"
    ]

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    #informs the reader of current snakemake progress
    log = st.empty()
    output = ""

    for line in process.stdout:
        output += line
        log.text_area(
            "Snakemake output",
            output,
            height=400
        )

