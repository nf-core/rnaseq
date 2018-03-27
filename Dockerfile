FROM continuumio/miniconda
MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
LABEL authors="phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for the NGI-RNAseq pipeline"

COPY environment.yml /
RUN conda update -n base conda && \
    conda env create -f /environment.yml && \
    conda clean -a
ENV PATH /opt/conda/envs/nfcore-rnaseq/bin:$PATH
