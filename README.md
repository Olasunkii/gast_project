# gast_project

This project automates the retrieval of public sequences and metadata for *Klebsiella pneumoniae* and is organized into four phases: data extraction, data cleaning and annotation, data integration, and data preparation for machine-learning.

**Input**
A user-adjusted configuration and sufficient compute resources (minimum eight CPU cores).

**Output**
A curated dataset suitable as input for machine-learning models, integrating draft genomes, detected resistance genes, antibiotic susceptibility data and host metadata.

## Workflow walkthrough

The workflow retrieves host and sample metadata and corresponding raw reads, preprocesses those reads with fastp, trims residual adapters with Trim Galore, evaluates read quality with FastQC, assembles draft genomes with Unicycler, annotates assemblies using Bakta, assesses genome quality with CheckM2, identifies resistance and associated genes with AMRFinderPlus, and continues integrating the draft genome, host metadata , and detected resistance genes into a single integrated CSV produced by a python script. Finally the dataset is prepared for machine learning by standaridizing, formatting and data partitioning. There are multiple checks where the user is able to document which isolates fail to meet standards, focusing on carbapenem phenotype breakpoints and genome completeness.

<p align="center">
<img src="overview_workflow.png" width="600">
</p>

**Overview of bioinformatics tool**

-   fastp
    -   Read preprocessing: quality filtering, adapter trimming, deduplication, and QC reporting.

-   Trim Galore
    -   Adapter and low-quality base trimming.

-   FastQC
    -   Quality assessment of raw reads, providing per-base metrics and contamination indicators.

-   Unicycler
    -   Short-read assembly of bacterial genomes into draft assemblies. Unicycler can also assemble long-read.

-   Bakta
    -   Automated annotation of bacterial draft genomes, generating gene and protein feature sets.

-   CheckM2
    -   Genome quality estimation using completeness and contamination metrics from marker-gene models.

-   AMRFinderPlus
    -   Detection of antimicrobial resistance genes, stress-response determinants, and selected virulence factors.


## Setup

Clone the repository, install the dependencies listed in `requirements.txt`, and place the required reference databases as described in the prerequisite section. The `setup_gast.py` script establishes the directory structure needed for the workflow.
Use at least *eight* CPU cores when processing fewer than ten samples. For larger batches, *sixteen* or more CPU cores are recommended.

Adjust the configuration file before execution. Required fields:

* `retmax`: number of samples to retrieve
* `email`: contact address for NCBI data retrieval
* `organism`: default is *Klebsiella pneumoniae*; modify this value to change to a different organism


### Prerequisites

* Install required packages for automated data retrieval.
* Download the associated databases (Checkm2 and bakta).

**Install required packages**
Use the standard pip installation command to install the required packages stated in the requirements.txt file:

```
pip install -r requirements.txt
```


**Reference Database Preparation**

Install and unpack the Bakta database:

```
wget https://zenodo.org/record/14916843/files/db-light.tar.xz
mkdir -p bakta_db
tar -xf db-light.tar.xz -C bakta_db #unpack reference database
mv bakta_db/db-light/* bakta_db/ #move database to right folder structure
rmdir bakta_db/db-light
```

Install and unpack the CheckM2 database:

```
wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz
mkdir -p checkm2_db
tar -xzf checkm2_database.tar.gz -C checkm2_db #unpack reference database
mv checkm2_db/CheckM2_database/* checkm2_db/ #move database to right folder structure
rmdir checkm2_db/CheckM2_database
```
The AMRFinderPlus database is updated automatically within its Snakemake rule; the rule writes a designated output file (*amrfinder_db_ready.txt*) that marks the database as current, and removing that file forces Snakemake to perform the database update again.


## Running pipeline via gui

Start gui by
```
streamlit run gui.py
```
Fill in the inputs as desired. This will be converter into a configuration file (config_parameter.yaml). Press the button to launch the Snakemake pipeline.

**Start pipeline via command line**

If not, the configuration parameters must be written exactly as required by the user. You can use `config_parameter_example.yaml` as a template; it must be followed strictly for the pipeline to run correctly. Then start the pipeline by


```
snakemake --cores #numb_cores --use-conda --latency-wait 30
```