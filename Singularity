From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
    DESCRIPTION Container image containing all requirements for the nf-core/rnaseq pipeline
    VERSION 1.5dev

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
    /opt/conda/bin/conda clean -a
