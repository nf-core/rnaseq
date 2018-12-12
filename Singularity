From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
    DESCRIPTION Singularity image containing all requirements for the nf-core/rnaseq pipeline
    VERSION 1.2

%environment
    PATH=/opt/conda/envs/nf-core-rnaseq-1.2/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
