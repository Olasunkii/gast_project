# ==============================
# Full WGS Pipeline Image
# ==============================
FROM ubuntu:22.04

LABEL maintainer="olasu-bro@omics_research"
LABEL description="Full WGS pipeline container (SRA → FASTQ → QC → Assembly → AMR)"

ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------
# System tools
# ------------------------------
RUN apt-get update && apt-get install -y \
    wget curl git python3 python3-pip \
    gzip pigz \
    sra-toolkit fastp \
    build-essential \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ==============================
# Install SPAdes 3.15.5
# ==============================
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz && \
    tar -xzf SPAdes-3.15.5-Linux.tar.gz && \
    mv SPAdes-3.15.5-Linux /opt/spades && \
    ln -s /opt/spades/bin/spades.py /usr/local/bin/spades.py && \
    ln -s /opt/spades/bin/metaspades.py /usr/local/bin/metaspades.py && \
    rm SPAdes-3.15.5-Linux.tar.gz


# ==============================
# Install AMRFinderPlus
# ==============================
# RUN wget https://github.com/ncbi/amr/releases/download/v4.0.23/amrfinder-4.0.23.tar.gz && \
#     tar -xzf amrfinder-4.0.23.tar.gz && \
#     mv amrfinder-4.0.23 /opt/amrfinder && \
#     ln -s /opt/amrfinder/amrfinder /usr/local/bin/amrfinder && \
#     rm amrfinder-4.0.23.tar.gz && \
#     amrfinder -u

# ------------------------------
# Python dependencies
# ------------------------------
COPY requirements.txt /tmp/
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# ------------------------------
# Copy pipeline scripts
# ------------------------------
WORKDIR /app
COPY scripts/ ./scripts/

RUN mkdir -p data/metadata data/sequences results

ENTRYPOINT ["python3", "scripts/sra_extractor.py"]



#Steps when starting this project
# 0. Start an interactive session
    # docker run -it --rm -v $(pwd):/app --entrypoint bash wgs_pipeline:latest

# 1. If you just restarted your machine, check that Docker is running:
    # sudo systemctl status docker
    # If inactive, start with:
        #sudo systemctl start docker

# 2. When updates are made, build with 
    # docker build -t wgs_pipeline:latest .

# 3. Then confirm everything installed:
    # docker run -it --rm wgs_pipeline:latest prefetch --version
    # docker run -it --rm wgs_pipeline:latest fastp --version

##### IMPORTANTEE ############
#1. Start Docker desktop
#2. Run the following:
    # docker build -t wgs_pipeline:latest .
    # docker run -it --rm -v $(pwd):/app --entrypoint bash wgs_pipeline:latest

# Install sankemake in container and pulp
    # pip install snakemake
    # pip install pulp==2.7.0
