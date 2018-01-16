#!/usr/bin/env bash

script_path="../main.nf"
PROFILE=docker

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    --path)
    script_path=$2
    shift # past argument
    shift # past value
    ;;
    -p|--profile)
    PROFILE=$2
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    shift # past argument
    ;;
  esac
done

data_path="/tmp"
if [ -d "./test_data" ]
then
    data_path="./test_data"
    echo "Found data directory in current working directory, using ./test_data/"
fi

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
docker -v >/dev/null 2>&1 || { echo >&2 "I require docker, but it's not installed. Visit https://www.docker.com/products/overview#/install_the_platform  ."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

data_dir=${data_path}/ngi-rna_test_set
if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Downloading test set..."
    curl https://export.uppmax.uu.se/b2013064/test-data/ngi-rna_test_set.tar.bz2 > ${data_path}/ngi-rna_test_set.tar.bz2
    echo "Unpacking test set..."
    tar xvjf ${data_path}/ngi-rna_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

run_name="Test RNA Run: "$(date +%s)

cmd="nextflow run $script_path -resume -name \"$run_name\" -profile $PROFILE --max_memory '7.GB' --max_cpus 2 --max_time '48.h' --gtf ${data_dir}/genes.gtf --bed12 ${data_dir}/genes.bed --hisat2_index ${data_dir}/r64/ --aligner hisat2 --singleEnd --reads \"${data_dir}/*.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-----"
eval $cmd
