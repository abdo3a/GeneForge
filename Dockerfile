# Dockerfile for GeneForge pipeline
# Base image: Ubuntu 22.04 for compatibility with bioinformatics tools
FROM ubuntu:22.04

# Metadata
LABEL maintainer="Your Name <your.email@example.com>"
LABEL description="Docker image for GeneForge, a Nextflow pipeline for gene prediction and functional annotation"
LABEL version="1.0"

# Set non-interactive frontend to avoid prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    openjdk-11-jre \
    python3 \
    python3-pip \
    python3-dev \
    perl \
    libperl-dev \
    libxml2-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    samtools \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Set Python3 as default
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install Perl modules for BRAKER3 and funannotate
RUN cpan -i XML::Parser && \
    cpan -i YAML && \
    cpan -i Hash::Merge && \
    cpan -i Logger::Simple && \
    cpan -i Parallel::ForkManager && \
    cpan -i MCE::Mutex && \
    cpan -i JSON && \
    cpan -i Bio::Perl

# Install Trimmomatic (v0.39)
RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic && \
    rm Trimmomatic-0.39.zip && \
    chmod +x /opt/trimmomatic/trimmomatic-0.39.jar

# Install STAR (v2.7.11b)
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    cd STAR-2.7.11b/source && \
    make STAR && \
    mv STAR /usr/local/bin/ && \
    cd ../.. && \
    rm -rf STAR-2.7.11b 2.7.11b.tar.gz

# Install StringTie (v2.2.3)
RUN wget https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz && \
    tar -xzf stringtie-2.2.3.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.3.Linux_x86_64/stringtie /usr/local/bin/ && \
    rm -rf stringtie-2.2.3.Linux_x86_64 stringtie-2.2.3.Linux_x86_64.tar.gz

# Install tRNAscan-SE (v2.0.12)
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz && \
    tar -xzf trnascan-se-2.0.12.tar.gz && \
    cd tRNAscan-SE-2.0.12 && \
    ./configure && make && make install && \
    cd .. && \
    rm -rf tRNAscan-SE-2.0.12 trnascan-se-2.0.12.tar.gz

# Install EukHighConfidenceFilter (for tRNAscan-SE)
RUN wget https://github.com/Gaius-Augustus/EukHighConfidenceFilter/archive/refs/heads/master.zip && \
    unzip master.zip && \
    mv EukHighConfidenceFilter-master /opt/EukHighConfidenceFilter && \
    chmod +x /opt/EukHighConfidenceFilter/EukHighConfidenceFilter && \
    ln -s /opt/EukHighConfidenceFilter/EukHighConfidenceFilter /usr/local/bin/EukHighConfidenceFilter && \
    rm master.zip

# Install funannotate (v1.8.18)
RUN pip install --no-cache-dir funannotate==1.8.18 && \
    mkdir -p /opt/funannotate_db && \
    chmod -R a+w /opt/funannotate_db

# Install BUSCO (v5.7.1)
RUN pip install --no-cache-dir busco==5.7.1

# Install Phobius (v1.01)
RUN wget http://phobius.sbc.su.se/data/phobius101_linux.tar.gz && \
    tar -xzf phobius101_linux.tar.gz && \
    mv phobius /opt/phobius && \
    chmod +x /opt/phobius/phobius.pl && \
    ln -s /opt/phobius/phobius.pl /usr/local/bin/phobius.pl && \
    rm phobius101_linux.tar.gz

# Install InterProScan (v5.71-103.0)
RUN wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-103.0/interproscan-5.71-103.0-64-bit.tar.gz && \
    tar -xzf interproscan-5.71-103.0-64-bit.tar.gz && \
    mv interproscan-5.71-103.0 /opt/interproscan && \
    rm interproscan-5.71-103.0-64-bit.tar.gz && \
    chmod +x /opt/interproscan/interproscan.sh && \
    ln -s /opt/interproscan/interproscan.sh /usr/local/bin/interproscan.sh

# Install eggNOG-mapper (v2.1.12)
RUN git clone https://github.com/eggnogdb/eggnog-mapper.git && \
    cd eggnog-mapper && \
    git checkout v2.1.12 && \
    pip install --no-cache-dir -r requirements.txt && \
    mv emapper.py /usr/local/bin/ && \
    mkdir -p /opt/eggnog-mapper && \
    mv * /opt/eggnog-mapper/ && \
    cd .. && \
    rm -rf eggnog-mapper

# Install AGAT (v1.4.0)
RUN pip install --no-cache-dir agat==1.4.0

# Install gffread (v0.12.7)
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar -xzf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin/ && \
    rm -rf gffread-0.12.7.Linux_x86_64 gffread-0.12.7.Linux_x86_64.tar.gz

# Install gtf_genome_to_cdna_fasta.pl
RUN wget https://raw.githubusercontent.com/gpertea/gffread/master/scripts/gtf_genome_to_cdna_fasta.pl && \
    mv gtf_genome_to_cdna_fasta.pl /usr/local/bin/ && \
    chmod +x /usr/local/bin/gtf_genome_to_cdna_fasta.pl

# Set environment variables for funannotate
ENV FUNANNOTATE_DB=/opt/funannotate_db
ENV PATH=/usr/local/bin:/opt/trimmomatic:/opt/phobius:/opt/interproscan:/opt/eggnog-mapper:$PATH

# Create working directory
RUN mkdir -p /workspace
WORKDIR /workspace

# Clean up
RUN apt-get clean && \
    rm -rf /tmp/* /var/tmp/*
