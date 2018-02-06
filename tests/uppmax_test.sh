#!/usr/bin/env bash

script_path="../main.nf"
if [ -z $1] 
then
    echo "No argument given, going to try to run ../main.nf"
else
    script_path=$1
fi

data_path="./test_data"

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

data_dir=${data_path}/ngi-rna_test_set
if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Unpacking test set..."
    mkdir -p ${data_path}
    tar xvjf ngi-rna_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

run_name="Test RNA Run: "$(date +%s)
cmd="nextflow run $script_path -resume -name \"$run_name\" -profile uppmax_devel --gtf ${data_dir}/genes.gtf --bed12 ${data_dir}/genes.bed --fasta ${data_dir}/genome.fa --singleEnd --reads \"${data_dir}/*.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-----"
eval $cmd
