# Use an official Ubuntu base image
FROM ubuntu:22.04

WORKDIR /app

# Install necessary packages
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    bzip2 \
    git \
    build-essential \
    libgl1-mesa-glx \
    libxrender1 \
    default-jre \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Install nextflow
RUN curl -s https://get.nextflow.io | bash

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && /opt/conda/bin/conda clean -a

# Install Mamba
RUN /opt/conda/bin/conda install -c conda-forge mamba
# Add Conda to PATH
ENV PATH=/opt/conda/bin:$PATH
# Copy environment files
COPY envs /app/envs

# Create Conda environments
RUN /opt/conda/bin/mamba env create -f /app/envs/genome-tools.yml
RUN /opt/conda/bin/mamba env create -f /app/envs/typing.yml
RUN /opt/conda/bin/mamba env create -f /app/envs/abricate.yml
RUN /opt/conda/bin/mamba env create -f /app/envs/dragonflye.yml

# Copy SeqSero2 database files
COPY seqsero2_db/antigens.pickle /opt/conda/envs/typing/database/antigens.pickle
COPY seqsero2_db/H_and_O_and_specific_genes.fasta /opt/conda/envs/typing/database/H_and_O_and_specific_genes.fasta

# Copy the script to set up the Abricate database
COPY setup_abricate_db.sh /app/

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]

# Activate the environment and set up the Abricate database
RUN source activate abricate && bash /app/setup_abricate_db.sh

# Copy the logo into the image
COPY images/fastpaslogo.svg /app/images/fastpaslogo.svg

# Copy Nextflow script and config file
COPY bacterial_assembly_pipeline.nf /app/
COPY nextflow.config /app/
# Copy subworkflows folder
COPY modules /app/modules/
COPY subworkflows /app/subworkflows/
# Set the default RUN_NAME environment variable
ENV RUN_NAME="default"
WORKDIR /input_output/${RUN_NAME}
# Set the entrypoint to run Nextflow
ENTRYPOINT ["/app/nextflow"]
