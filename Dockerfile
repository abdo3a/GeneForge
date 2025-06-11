FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    python3 \
    python3-pip \
    perl \
    openjdk-11-jre \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    unzip \
    libcairo2-dev \
    libpango1.0-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar && \
    rm Trimmomatic-0.39.zip

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make STAR && \
    mv STAR /usr/local/bin/ && \
    cd ../.. && \
    rm -rf STAR-2.7.10b 2.7.10b.tar.gz

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar -xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure && make && make install && \
    cd .. && rm -rf samtools-1.19 samtools-1.19.tar.bz2

# Install StringTie
RUN wget https://github.com/gpertea/stringtie/archive/refs/tags/v2.2.1.tar.gz && \
    tar -xzf v2.2.1.tar.gz && \
    cd stringtie-2.2.1 && \
    make release && \
    mv stringtie /usr/local/bin/ && \
    cd .. && rm -rf stringtie-2.2.1 v2.2.1.tar.gz

# Install tRNAscan-SE
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.9.tar.gz && \
    tar -xzf trnascan-se-2.0.9.tar.gz && \
    cd tRNAscan-SE-2.0.9 && \
    ./configure && make && make install && \
    cd .. && rm -rf tRNAscan-SE-2.0.9 trnascan-se-2.0.9.tar.gz

# Install funannotate and dependencies
RUN pip3 install --no-cache-dir \
    biopython \
    requests \
    psutil \
    numpy \
    pandas \
    scikit-learn \
    goatools \
    natsort \
    matplotlib \
    seaborn \
    && wget https://github.com/nextgenusfs/funannotate/archive/refs/tags/1.8.15.tar.gz && \
    tar -xzf 1.8.15.tar.gz && \
    cd funannotate-1.8.15 && \
    python3 setup.py install && \
    cd .. && rm -rf funannotate-1.8.15 1.8.15.tar.gz

# Install BUSCO
RUN pip3 install --no-cache-dir busco==5.4.7

# Install gtf_genome_to_cdna_fasta.pl (from TransDecoder)
RUN wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz && \
    tar -xzf TransDecoder-v5.7.1.tar.gz && \
    mv TransDecoder-TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl /usr/local/bin/ && \
    chmod +x /usr/local/bin/gtf_genome_to_cdna_fasta.pl && \
    rm -rf TransDecoder-TransDecoder-v5.7.1 TransDecoder-v5.7.1.tar.gz

# Install intervaltree for stringtie2utr.py
RUN pip3 install --no-cache-dir intervaltree

# Install AGAT for merging GFF files and fixing overlaps
RUN wget https://github.com/NBISweden/AGAT/archive/refs/tags/v1.2.0.tar.gz && \
    tar -xzf v1.2.0.tar.gz && \
    cd AGAT-1.2.0 && \
    perl Makefile.PL && \
    make && make install && \
    cd .. && rm -rf AGAT-1.2.0 v1.2.0.tar.gz

# Install gffread
RUN wget https://github.com/gpertea/gffread/archive/refs/tags/v0.12.7.tar.gz && \
    tar -xzf v0.12.7.tar.gz && \
    cd gffread-0.12.7 && \
    make && \
    mv gffread /usr/local/bin/ && \
    cd .. && rm -rf gffread-0.12.7 v0.12.7.tar.gz

# Install GenomeTools for gt gff3validator
RUN wget http://genometools.org/pub/genometools-1.6.5.tar.gz && \
    tar -xzf genometools-1.6.5.tar.gz && \
    cd genometools-1.6.5 && \
    make && make install && \
    cd .. && rm -rf genometools-1.6.5 genometools-1.6.5.tar.gz

# Install Phobius
RUN wget http://phobius.sbc.su.se/data/phobius101.tar.gz && \
    tar -xzf phobius101.tar.gz && \
    mv phobius101 /opt/phobius && \
    ln -s /opt/phobius/phobius.pl /usr/local/bin/phobius.pl && \
    rm phobius101.tar.gz

# Install eggNOG-mapper
RUN wget https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.1.12.tar.gz && \
    tar -xzf 2.1.12.tar.gz && \
    cd eggnog-mapper-2.1.12 && \
    python3 setup.py install && \
    cd .. && rm -rf eggnog-mapper-2.1.12 2.1.12.tar.gz

# Set funannotate database environment variable
ENV FUNANNOTATE_DB=/opt/funannotate_db
RUN mkdir -p ${FUNANNOTATE_DB} && \
    funannotate setup -d ${FUNANNOTATE_DB} -b protists

# Set LD_LIBRARY_PATH for GenomeTools
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Set working directory
WORKDIR /workspace
