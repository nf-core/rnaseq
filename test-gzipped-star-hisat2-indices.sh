#!/bin/bash
set -ev
bundle exec rake:units
if [ "${SUITE}" = "test_gz" ]; then
  # Run without premade STAR indices - made within nextflow
  nextflow run ${TRAVIS_BUILD_DIR} -profile $SUITE,docker --star_index false --skipQC

  # Run with HiSat2 and no premade indices -- do-it-yourself reference
  nextflow run ${TRAVIS_BUILD_DIR} -profile $SUITE,docker --aligner hisat2 --hisat2_index false --skipQC
fi
