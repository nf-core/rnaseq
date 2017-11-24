FROM openjdk:8

LABEL authors="phil.ewels@scilifelab.se,rickard.hammaren@scilifelab.se,denis.moreno@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAseq pipeline"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        make \
        python-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py

RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o /opt/fastqc_v0.11.5.zip && \
    unzip /opt/fastqc_v0.11.5.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v0.11.5.zip

# Install bedops
RUN mkdir /opt/bedops && \
    curl -fsSL https://github.com/bedops/bedops/releases/download/v2.4.20/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 -o /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 && \
    tar xvjf /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 -C /opt/bedops && \
    ln -s /opt/bedops/bin/* /usr/local/bin/ && \
    rm /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2

# Install cutadapt
RUN pip install cutadapt

# Install TrimGalore
RUN mkdir /opt/TrimGalore && \
    curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.2.zip -o /opt/TrimGalore/trim_galore_v0.4.2.zip && \
    unzip /opt/TrimGalore/trim_galore_v0.4.2.zip -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/trim_galore_v0.4.2.zip

# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# Install RSeQC
RUN pip install RSeQC

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o /opt/samtools-1.3.1.tar.bz2 && \
    tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.3.1 && \
    make && \
    make install && \
    rm /opt/samtools-1.3.1.tar.bz2

# Install PreSeq
RUN curl -fsSL http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 -o /opt/preseq_linux_v2.0.tar.bz2 && \
    tar xvjf /opt/preseq_linux_v2.0.tar.bz2 -C /opt/ && \
    ln -s /opt/preseq_v2.0/preseq /usr/local/bin/preseq && \
    ln -s /opt/preseq_v2.0/bam2mr /usr/local/bin/bam2mr && \
    rm /opt/preseq_linux_v2.0.tar.bz2 && \
    # Make sure that libgsl.so.0 exists beause PreSeq links to that
    ln -s /usr/lib/x86_64-linux-gnu/libgsl.so /lib/x86_64-linux-gnu/libgsl.so.0

# Install PicardTools
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.0.1/picard-tools-2.0.1.zip -o /opt/picard-tools-2.0.1.zip && \
    unzip /opt/picard-tools-2.0.1.zip -d /opt/ && \
    rm /opt/picard-tools-2.0.1.zip
ENV PICARD_HOME /opt/picard-tools-2.0.1

# Install R
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/R-3.4.2.tar.gz -o /opt/R-3.4.2.tar.gz && \
    tar xvzf /opt/R-3.4.2.tar.gz -C /opt/ && \
    cd /opt/R-3.4.2 && \
    ./configure && \
    make && \
    make install && \
    rm /opt/R-3.4.2.tar.gz

# Install R Packages v2
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("Rsubread", "dupRadar", "limma", "lattice", "locfit", "edgeR", "chron", "data.table", "gtools", "gdata", "bitops", "caTools", "gplots", "markdown"))' >> /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir -p  /usr/local/lib/R/site-library

# Install featureCounts
RUN curl -fsSL http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz -o /opt/subread-1.5.1-Linux-x86_64.tar.gz && \
    tar xvzf /opt/subread-1.5.1-Linux-x86_64.tar.gz -C /opt/ && \
    ln -s /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts && \
    rm /opt/subread-1.5.1-Linux-x86_64.tar.gz

# Install StringTie
RUN curl -fsSL http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3.Linux_x86_64.tar.gz -o /opt/stringtie-1.3.3.Linux_x86_64.tar.gz && \
    tar xvzf /opt/stringtie-1.3.3.Linux_x86_64.tar.gz -C /opt/ && \
    ln -s /opt/stringtie-1.3.3.Linux_x86_64/stringtie /usr/local/bin/stringtie && \
    rm /opt/stringtie-1.3.3.Linux_x86_64.tar.gz

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git

# Install Hisat2
RUN git clone https://github.com/infphilo/hisat2.git /opt/hisat2 && \
    cd /opt/hisat2/ && \
    make && \
    cp /opt/hisat2/hisat2 /opt/hisat2/hisat2-align-s /opt/hisat2/hisat2-align-l /opt/hisat2/hisat2-build /opt/hisat2/hisat2-build-s /opt/hisat2/hisat2-build-l /opt/hisat2/hisat2-inspect /opt/hisat2/hisat2-inspect-s /opt/hisat2/hisat2-inspect-l /usr/local/bin/ && \
    cp /opt/hisat2/*.py /usr/local/bin

# Create UPPMAX root directories
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw
