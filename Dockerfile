FROM openjdk:8

LABEL authors="phil.ewels@scilifelab.se,rickard.hammaren@scilifelab.se,denis.moreno@scilifelab.se"

#Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update
RUN apt-get install -y libreadline-dev
RUN apt-get install -y gcc 
RUN apt-get install -y make
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y python-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y gfortran
RUN apt-get install -y g++
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y libpcre3-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libgsl0-dev
RUN wget -O /opt/get-pip.py --no-check-certificate https://bootstrap.pypa.io/get-pip.py
RUN python /opt/get-pip.py
RUN rm /opt/get-pip.py

RUN wget -O /opt/fastqc_v0.11.5.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip /opt/fastqc_v0.11.5.zip -d /opt/
RUN chmod 755 /opt/FastQC/fastqc
RUN ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc
RUN rm /opt/fastqc_v0.11.5.zip

#Install bedops
RUN mkdir /opt/bedops
RUN wget -q -O /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 https://github.com/bedops/bedops/releases/download/v2.4.20/bedops_linux_x86_64-v2.4.20.v2.tar.bz2
RUN tar xvjf /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 -C /opt/bedops
RUN ln -s /opt/bedops/bin/* /usr/local/bin/
RUN rm /opt/bedops_linux_x86_64-v2.4.20.v2.tar.bz2

#Install cutadapt
RUN pip install cutadapt

#Install TrimGalore
RUN mkdir /opt/TrimGalore
RUN wget -q -O /opt/TrimGalore/trim_galore_v0.4.2.zip http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.2.zip
RUN unzip /opt/TrimGalore/trim_galore_v0.4.2.zip -d /opt/TrimGalore
RUN ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore
RUN rm /opt/TrimGalore/trim_galore_v0.4.2.zip

#Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR
RUN ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR
RUN ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

#Install RSeQC
RUN pip install RSeQC

#Install SAMTools
RUN wget -q -O /opt/samtools-1.3.1.tar.bz2 https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/
RUN cd /opt/samtools-1.3.1;make;make install
RUN rm /opt/samtools-1.3.1.tar.bz2 

#Install PreSeq
RUN wget -q -O /opt/preseq_linux_v2.0.tar.bz2 http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2
RUN tar xvjf /opt/preseq_linux_v2.0.tar.bz2 -C /opt/
RUN ln -s /opt/preseq_v2.0/preseq /usr/local/bin/preseq
RUN ln -s /opt/preseq_v2.0/bam2mr /usr/local/bin/bam2mr
RUN rm /opt/preseq_linux_v2.0.tar.bz2

#Install PicardTools
RUN wget -q -O /opt/picard-tools-2.0.1.zip https://github.com/broadinstitute/picard/releases/download/2.0.1/picard-tools-2.0.1.zip
RUN unzip /opt/picard-tools-2.0.1.zip -d /opt/
RUN rm /opt/picard-tools-2.0.1.zip
ENV PICARD_HOME /opt/picard-tools-2.0.1

#Install R
RUN wget -q -O /opt/R-3.2.3.tar.gz https://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz
RUN tar xvzf /opt/R-3.2.3.tar.gz -C /opt/
RUN cd /opt/R-3.2.3;./configure;make;make install
RUN rm /opt/R-3.2.3.tar.gz 

#Install R Packages v2
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r
RUN echo 'biocLite()' >> /opt/packages.r
RUN echo 'biocLite(c("Rsubread", "dupRadar", "limma", "lattice", "locfit", "edgeR", "chron", "data.table", "gtools", "gdata", "bitops", "caTools", "gplots"))' >> /opt/packages.r
RUN Rscript /opt/packages.r
RUN mkdir /usr/local/lib/R/site-library

#Install featureCounts
RUN wget -q -O /opt/subread-1.5.1-Linux-x86_64.tar.gz http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz
RUN tar xvzf /opt/subread-1.5.1-Linux-x86_64.tar.gz -C /opt/
RUN ln -s /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts
RUN rm /opt/subread-1.5.1-Linux-x86_64.tar.gz

#Install StringTie
RUN wget -q -O /opt/stringtie-1.3.0.Linux_x86_64.tar.gz  http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
RUN tar xvzf /opt/stringtie-1.3.0.Linux_x86_64.tar.gz -C /opt/
RUN ln -s /opt/stringtie-1.3.0.Linux_x86_64/stringtie /usr/local/bin/stringtie
RUN rm /opt/stringtie-1.3.0.Linux_x86_64.tar.gz

#Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git
