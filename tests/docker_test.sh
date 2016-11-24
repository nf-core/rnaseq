#!/usr/bin/env bash

script_path="../main.nf"
if [ -z $1]
then
    echo "No argument given, going to try to run ../main.nf"

else
    script_path=$1
fi

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed.  Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed.  Aborting."; exit 1; }
docker -v >/dev/null 2>&1 || { echo >&2 "I require docker, but it's not installed. Visit https://www.docker.com/products/overview#/install_the_platform  ."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

echo "Downloading test set..."
curl https://s3-eu-west-1.amazonaws.com/ngi-rnaseq/ngi-rna_test_set.tar.bz2 > /tmp/ngi-rna_test_set.tar.bz2
echo "Done"
if [ -d /tmp/ngi-rna_test_set ]
then
    echo "Removing previous test set..."
    rm -r /tmp/ngi-rna_test_set
    echo "Done"
fi
echo "Unpacking test set..."
tar xvjf /tmp/ngi-rna_test_set.tar.bz2 -C /tmp
echo "Done"
echo "Starting nextflow..."
nextflow run $script_path -profile docker_test --reads '/tmp/ngi-rna_test_set/*.fastq.gz'

