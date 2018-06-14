#!/bin/bash
set -euo pipefail

CONTAINER_ENGINE="singularity"
SGT_VER="2.5.1"

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -e|--engine)
    CONTAINER_ENGINE="$2"
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

# Install Nextflow
mkdir /tmp/nextflow
cd /tmp/nextflow
wget -qO- get.nextflow.io | bash
sudo ln -s /tmp/nextflow/nextflow /usr/local/bin/nextflow

# Install nf-core/tools
git clone https://github.com/nf-core/tools.git /tmp/nf-core-tools
cd /tmp/nf-core-tools
pip install --user -e .

# Install Singularity
if [[ "$CONTAINER_ENGINE" = singularity ]]
then
  sudo apt-get install libarchive-dev squashfs-tools
  mkdir /tmp/singularity
  cd /tmp/singularity
  wget https://github.com/singularityware/singularity/releases/download/$SGT_VER/singularity-$SGT_VER.tar.gz
  tar xvf singularity-$SGT_VER.tar.gz
  cd singularity-$SGT_VER
  ./configure --prefix=/usr/local
  make
  sudo make install
  cd ${TRAVIS_BUILD_DIR}
  rm -rf /tmp/singularity
fi
