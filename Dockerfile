FROM nfcore/base:dev
LABEL authors="Phil Ewels, Rickard Hammarén" \
      description="Docker image containing all requirements for nf-core/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnaseq-1.4/bin:$PATH
