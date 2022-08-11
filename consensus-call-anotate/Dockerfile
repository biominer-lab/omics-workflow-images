# Consensus annotate snv and indel based on method of PCAWG
# Code from github of SGA Annotate variants
#
# Version 0.1
#

FROM ubuntu:14.04
MAINTAINER Jun Shang <shangjunv@163.com> 
LABEL Description="Runs bam-based annotation on candidate SNV, indel calls" Vendor="OICR" Version="0.1"

VOLUME /data
WORKDIR /data

# get ubuntu packages
RUN apt-get update && \
    apt-get install -y \
    automake \
    perl \
    autotools-dev \
    build-essential \
    cmake \
    git \
    libhts-dev \
    libhts0 \
    libjemalloc-dev \
    libncurses5-dev \
    libsparsehash-dev \
    libz-dev \
    python \
    python-dev \
    python-pip \
    tabix \
    wget \
    zlib1g-dev 

# get pyvcf for mergevcf and annotate_from_readcounts.py
RUN pip install --upgrade pip && \
    pip install pyvcf

# build remaining dependencies:
# bcftools: sort
RUN cd /tmp && \
    wget https://github.com/samtools/bcftools/archive/refs/tags/1.15.1.tar.gz && \
    tar -xzvf 1.15.1.tar.gz && \
    rm -rf 1.15.1.tar.gz && \
    cd bcftools-1.15.1 && \
    make && \
    make install

# mergevcf: merge snv and indel vcf
RUN cd /tmp && \
    wget https://github.com/ljdursi/mergevcf/archive/0.2.tar.gz && \
    tar -xzf 0.2.tar.gz && \
    rm 0.2.tar.gz && \
    cd mergevcf-0.2 && \
    python setup.py install && \
    cd .. && \
    rm -rf mergevcf-0.2

# vt: manipulating and generating VCF files
RUN cd /tmp && \
    wget https://github.com/atks/vt/archive/0.5772.tar.gz && \
    tar -xzf 0.5772.tar.gz && \
    rm 0.5772.tar.gz && \
    cd vt-0.5772 && \
    make &&\
    mv vt /usr/local/bin && \
    cd .. &&\
    rm -rf vt-0.5772

# vcftools: Filtering out snv and indel
RUN cd /tmp && \
    wget -nv https://github.com/vcftools/vcftools/archive/refs/tags/v0.1.16.tar.gz && \
    tar -xzvf v0.1.16.tar.gz && \
    rm -rf v0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install

# bamtools: 
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xzvf v2.4.0.tar.gz && \
    rm -rf v2.4.0.tar.gz && \
    cd bamtools-2.4.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# bam-readcount: for SNV annotation
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/genome/bam-readcount/archive/v0.7.4.tar.gz && \
    tar -xzvf v0.7.4.tar.gz && \
    rm v0.7.4.tar.gz && \
    cd bam-readcount-0.7.4 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# samtools - for indexing reference, etc
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -xjvf samtools-1.15.tar.bz2 && \
    rm samtools-1.15.tar.bz2 && \
    cd samtools-1.15 && \
    make prefix=/usr/local/ install && \
    cd .. && \
    rm -rf samtools-1.15

# vcflib - for tools like vcfbreakmulti
# set a fixed version for reproducibility
RUN mkdir -p /deps && \
    cd /deps && \
    git clone --recursive git://github.com/ekg/vcflib.git && \
    cd vcflib && \
    git checkout d453d91592fe8a74d92b49cd6c7cd73f79a8b70b && \
    make 

# build SGA
RUN mkdir -p /src && \
    cd /src && \
    wget https://github.com/jts/sga/archive/v0.10.14.tar.gz && \
    tar -xzvf v0.10.14.tar.gz && \
    cd sga-0.10.14/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/deps/bamtools-2.4.0 --with-jemalloc=/usr --prefix=/usr/local && \
    make && \
    make install

# Put auxilliary data in /usr/local/share
COPY indel.header /usr/local/share
COPY snv.header /usr/local/share

# Put scripts in /usr/local/bin
COPY indel_annotate.sh /usr/local/bin
COPY snv_annotate.sh /usr/local/bin
COPY annotate.sh /usr/local/bin
COPY annotate_from_readcounts.py /usr/local/bin

ENTRYPOINT ["/usr/local/bin/annotate.sh"]