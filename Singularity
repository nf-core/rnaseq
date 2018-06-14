From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
    DESCRIPTION Container image containing all requirements for the nf-core/rnaseq pipeline
    VERSION 1.5dev

%environment
    PATH=/opt/conda/envs/nfcore-rnaseq-1.5dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
