# sv merging mainly based on method of PCAWG
#
# Version 0.1
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
RUN conda-pack -n consensus-sv-merge -o /tmp/env.tar && \
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
LABEL Description="This image is used for merging sv"
ENV PYTHON_VERSION 2.7*
ENV PYVCF_VERSION 0.6.8
ENV BEDTOOLS_VERSION 2.29.2
ENV PANDAS_VERSION 0.24.2
ENV NUMPY_VERSION 1.16.5
ENV PYBEDTOOLS_VERSION 0.7.10
ENV MATPLOTLIB_VERSION 2.2.3
ENV NETWORKX_VERSION 2.1

# Put scripts in /usr/local/bin
COPY dellyVcf2Bedpe.py /usr/bin
COPY dRangerVcf2Bedpe.py /usr/bin
COPY gridssVcf2Bedpe.py /usr/bin
COPY mantaVcf2Bedpe.py /usr/bin
COPY SangerVcf2Bedpe.py /usr/bin
COPY smufinVcf2Bedpe.py /usr/bin
COPY snowmanVcf2Bedpe.py /usr/bin
COPY svabaVcf2Bedpe.py /usr/bin
COPY sv_merge_reorder_pairs.py /usr/bin
COPY sv_merge_graph.py /usr/bin

RUN chmod a+x /usr/bin/dellyVcf2Bedpe.py
RUN chmod a+x /usr/bin/dRangerVcf2Bedpe.py
RUN chmod a+x /usr/bin/gridssVcf2Bedpe.py
RUN chmod a+x /usr/bin/mantaVcf2Bedpe.py
RUN chmod a+x /usr/bin/SangerVcf2Bedpe.py
RUN chmod a+x /usr/bin/smufinVcf2Bedpe.py
RUN chmod a+x /usr/bin/snowmanVcf2Bedpe.py
RUN chmod a+x /usr/bin/svabaVcf2Bedpe.py
RUN chmod a+x /usr/bin/sv_merge_reorder_pairs.py
RUN chmod a+x /usr/bin/sv_merge_graph.py
