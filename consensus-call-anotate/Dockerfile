# Consensus annotate snv and indel based on method of PCAWG
#
# Version 0.2.4
#

# The build-stage image:
FROM continuumio/miniconda3 AS build

# Install conda-pack:
RUN conda install -c conda-forge conda-pack mamba

# Install the package as normal:
COPY environment.yml .
COPY condarc ~/.condarc
RUN mamba env create -f environment.yml

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n consensus-call-anotate -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack

# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM debian:latest AS runtime

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

# Copy /venv from the previous stage
COPY --from=build /venv /venv
ENV PATH=/venv/bin:$PATH

WORKDIR /data

# ---------------------------------------------------
# Customized based on your situations.
MAINTAINER "Jun Shang <shangjunv@163.com>"
LABEL Description="This image is used for merging snv and indel"
ENV PYTHON_VERSION 2.7*
ENV PYVCF_VERSION 0.6.8
ENV HTSLIB_VERSION 1.3.*
ENV SAMTOOLS_VERSION 1.6
ENV ZLIB_VERSION 1.2.8*
ENV VCFLIB_VERSION 1.0.0_rc1_16.05.18
ENV MERGEVCF_VERSION 1.0.1
ENV VT_VERSION 2015.11.10
ENV BAM_READCOUNT_VERSION 0.8
ENV VCFTOOLS_VERSION 0.1.15
ENV BAMTOOLS_VERSION 2.4.1
ENV SGA_VERSION 0.10.15

# Put auxilliary data in /usr/local/share
COPY indel.header /usr/local/share/indel.header
COPY snv.header /usr/local/share/snv.header

# Put scripts in /usr/local/bin
COPY indel_annotate.sh /usr/bin
COPY snv_annotate.sh /usr/bin
COPY annotate.sh /usr/bin
COPY annotate_from_readcounts.py /usr/bin

RUN chmod a+x /usr/bin/annotate.sh
RUN chmod a+x /usr/bin/annotate_from_readcounts.py