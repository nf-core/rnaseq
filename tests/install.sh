#!/bin/bash

# usage: install.sh nst
#   letters represent tools:
#     n = nextflow
#     s = singularity
#     t = nf-core/tools

# Install Nextflow
if [[ $1 = *n* ]]; then
  cd $HOME
  curl -fsSL get.nextflow.io | bash
  chmod +x nextflow
  sudo mv nextflow /usr/local/bin/
fi

# Install Singularity
if [[ $1 = *s* ]]; then
  sudo apt-get install squashfs-tools
  cd $HOME
  wget https://github.com/singularityware/singularity/releases/download/$SGT_VER/singularity-$SGT_VER.tar.gz
  tar xvf singularity-$SGT_VER.tar.gz
  cd singularity-$SGT_VER
  ./configure --prefix=/usr/local
  make
  sudo make install
  cd ..
  rm -rf singularity-$SGT_VER*
fi

# Install nf-core/tools
if [[ $1 = *t* ]]; then
    pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git
fi
