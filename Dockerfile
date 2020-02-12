FROM nfcore/base:1.8
LABEL authors="phil.ewels@scilifelab.se" \
      description="Docker image containing all requirements for the nfcore/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnaseq-1.4.3dev/bin:$PATH
RUN conda env export --name nf-core-rnaseq-1.4.3dev > nf-core-rnaseq-1.4.3dev.yml
