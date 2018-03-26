#!/usr/bin/env bash

# print_usage()
function print_usage {
  echo -e  "\nUsage:\t$0\n" \
    "\t\t[-a (aligner to use)\n" \
    "\t\t[-b (build genome references)\n" \
    "\t\t[-r (run in RRBS mode)\n" \
    "\t\t[-n (run in notrim mode)\n" \
    "\t\t[-u (run UPPMAX test)\n" \
    "\t\t[-t <test data directory>]\n" \
    "\t\t[-d <docker image>]\n" \
    "\t\t[-s <singularity image>]\n" \
    "\t\t[-h (show this help message)]" >&2 ;
}

# Check that we have required commands
curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

# Detect Travis fork for dockerhub image if we can
containerfl="-with-docker"
if [ ! -z "$TRAVIS_REPO_SLUG" ]; then
    echo "Detected repo as '$TRAVIS_REPO_SLUG'"
    dockerimg=$(echo "$TRAVIS_REPO_SLUG" | awk '{print tolower($0)}' | sed s/^nf-core/nfcore/)
    if [ ! -z "$TRAVIS_BRANCH" ] && [ "$TRAVIS_BRANCH" != "master"]; then
        dockerimg="$dockerimg:$TRAVIS_BRANCH"
    fi
    echo "Using docker image '$dockerimg'"
    containerfl="-with-docker $dockerimg"
fi

# Look for an existing test data directory
data_path="./test_data"
data_dir=${data_path}/ngi-rna_test_set

if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Unpacking test set..."
    tar xvjf ${data_path}/ngi-rna_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

# command line options
pipelinescript="../main.nf"
aligner="star"
profile="--max_cpus 2 --max_memory '6.GB' --max_time '48.h'"
customrefs=""

while getopts ":bpusht:d:c:a:" opt; do
  case $opt in
    a)
      echo "Using aligner $OPTARG" >&2
      aligner=$OPTARG
      ;;
    b)
      echo "Building genome references" >&2
      customrefs="--saveReference --fasta ${data_dir}/genome.fa"
      buildrefs=1
      ;;
    u)
      echo "Running UPPMAX config" >&2
      profile="-profile uppmax_devel"
      ;;
    s)
      echo "Running in singularity mode" >&2
      containerfl=${containerfl/-with-docker/-with-singularity}
      ;;
    t)
      echo "Test data path specified" >&2
      data_path=$OPTARG
      ;;
    d)
      echo "Using docker image $OPTARG" >&2
      containerfl="-with-docker $OPTARG"
      ;;
    c)
      echo "Using singularity image $OPTARG" >&2
      containerfl="-with-singularity $OPTARG"
      ;;
    h)
      print_usage
      exit
      ;;
    :)
      echo -e "\nOption -$OPTARG requires an argument." >&2
      print_usage
      exit 1;
      ;;
    \?)
      print_usage
      echo -e "\nInvalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


if [ ! -z "$customrefs" ]; then
    refs=$customrefs
elif [ -d "${data_dir}/star/" ] && [ "$aligner" == "star" ]; then
    refs="--star_index ${data_dir}/star/"
elif [ ! -d "${data_dir}/star/" ] && [ "$aligner" == "star" ]; then
    refs="--star_index results/reference_genome/star/"
    echo "Attempting to use a reference genome from previous run in ./results/reference_genome"
elif [ -d "${data_dir}/r64/" ] && [ "$aligner" == "hisat2" ]; then
    refs="--hisat2_index ${data_dir}/r64/ --fasta ${data_dir}/references/WholeGenomeFasta/genome.fa"
elif [ ! -d "${data_dir}/r64/" ] && [ "$aligner" == "hisat2" ]; then
    refs="--hisat2_index results/reference_genome/r64/"
    echo "Attempting to use a reference genome from previous run in ./results/reference_genome/"
fi


cmd="nextflow run $pipelinescript -resume --aligner $aligner $profile $containerfl --gtf ${data_dir}/genes.gtf --bed12 ${data_dir}/genes.bed $refs --singleEnd --reads \"${data_dir}/*.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-----"
eval $cmd
