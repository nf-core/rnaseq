FROM continuumio/miniconda
MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
LABEL authors="phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for the NGI-RNAseq pipeline"

RUN conda update -n base conda
COPY environment.yml /
RUN conda env create -f /environment.yml
ENV PATH /opt/conda/envs/nfcore-methylseq/bin:$PATH
