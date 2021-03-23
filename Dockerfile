FROM nfcore/base:1.13.2
LABEL authors="Phil Ewels, Rickard Hammarén" \
      description="Docker image containing all software requirements for the nf-core/rnaseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnaseq-3.1dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnaseq-3.1dev > nf-core-rnaseq-3.1dev.yml
