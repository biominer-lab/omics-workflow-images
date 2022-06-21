# The build-stage image:
FROM continuumio/miniconda3 AS build

# Install conda-pack:
RUN conda install -c conda-forge conda-pack mamba

# Install the package as normal:
COPY environment.yml .
COPY condarc ~/.condarc
RUN mamba env create -f environment.yml
RUN pip install pysam==0.6

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n dnaseq-gatk-snv-workflow -o /tmp/env.tar && \
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
MAINTAINER "Jingcheng Yang <yjcyxky@163.com>; Jun Shang <shangjunv@163.com>"
LABEL Description="This image is used for DNA-Seq SNV (GATK)."
LABEL HOME="https://github.com/biominer-lab/omics-workflows"
ENV FASTQC_VERSION 0.11.9
ENV FASTP_VERSION 0.23.2
ENV FASTQ_SCREEN_VERSION 0.15.2
ENV QUALIMAP_VERSION 2.2.2d
ENV BWA_VERSION 0.7.17
ENV GATK_VERSION 3.8
ENV HAP_PY_VERSION 0.3.10
ENV PICARD_VERSION 2.18.7
