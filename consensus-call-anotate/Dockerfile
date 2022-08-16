# Consensus annotate snv and indel based on method of PCAWG
# Code from github of SGA Annotate variants
#
# Version 0.1.7
#

FROM ubuntu:14.04
MAINTAINER Jun Shang <shangjunv@163.com> 
LABEL Description="Runs bam-based annotation on candidate SNV, indel calls" Vendor="OICR" Version="0.1.7"

VOLUME /data
WORKDIR /data

# get ubuntu packages
RUN apt-get update && \
    apt-get install -y \
    automake \
    autotools-dev \
    build-essential \
    gcc \
    git \
    perl \
    libhts-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
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
    pkg-config \
    zlib1g-dev 

# build remaining dependencies:
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz && \
    tar -xzvf cmake-3.20.0.tar.gz && \
    rm cmake-3.20.0.tar.gz && \
    cd cmake-3.20.0 && \
    ./bootstrap && \
    make && \
    make install

# vcflib - for tools like vcfbreakmulti
# set a fixed version for reproducibility
# RUN mkdir -p /deps && \
#     cd /deps && \
#     wget https://github.com/vcflib/vcflib/archive/refs/tags/v1.0.0.tar.gz && \
#     tar -xzvf v1.0.0.tar.gz && \
#     rm v1.0.0.tar.gz && \
#     cd vcflib-1.0.0 && \
#     make 

RUN mkdir -p /deps && \
    cd /deps && \
    git clone --recursive https://github.com/vcflib/vcflib.git && \
    cd vcflib && \
    mkdir -p build && cd build && \
    cmake  -DCMAKE_BUILD_TYPE=Debug -DZIG=OFF -DOPENMP=OFF .. && \
    cmake --build . && \
    cmake --install .

# HTSlib: 
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 && \
    tar -xjvf htslib-1.15.1.tar.bz2 && \
    rm -rf htslib-1.15.1.tar.bz2 && \
    cd htslib-1.15.1 && \
    autoreconf -i && \
    ./configure --disable-bz2 --disable-lzma --disable-libgsl && \
    make && \
    make install

# samtools - for indexing reference, etc
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -xjvf samtools-1.15.tar.bz2 && \
    rm samtools-1.15.tar.bz2 && \
    cd samtools-1.15 && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure --with-htslib=/deps/htslib-1.15.1 && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.15

# get pyvcf for mergevcf and annotate_from_readcounts.py
RUN pip install pyvcf -i http://pypi.douban.com/simple/

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
    ./autogen.sh && \
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

RUN chmod a+x /usr/local/bin/annotate.sh
RUN chmod a+x /usr/local/bin/annotate_from_readcounts.py