#!/bin/bash -l
#SBATCH -A b2013064
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J nfcore-rnaseq-standalone

######## DESCRIPTION ####
# This scipt was created in order to do some quick analysis with BAM files as input.
# It mimics the commands used in the RNA pipeline, but from aligned files onwards
# It was written to work on Uppmax's SLURM cluster Milou, but can easily be tweaked to work elsewere
# NB: You will need to change the project number at the top of this file ( -A [your-project] )

######## USAGE ##########
# sbatch RNA-pipeline-from-BAM.sh my_file.bam
bam=$1

################### Reference genome locations
bed12="/sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"
gtf="/sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
star="/sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/"

#################### RSeQC
echo "Starting RSeQC"
module load bioinfo-tools
module load rseqc
samtools index $bam
infer_experiment.py -i $bam -r $bed12 > $bam.infer_experiment.txt
RPKM_saturation.py -i $bam -r $bed12 -d $strandRule -o $bam.RPKM_saturation
junction_annotation.py -i $bam -o $bam.rseqc -r $bed12
bam_stat.py -i $bam 2> $bam.bam_stat.txt
junction_saturation.py -i $bam -o $bam.rseqc -r $bed12 2> $bam.junction_annotation_log.txt
inner_distance.py -i $bam -o $bam.rseqc -r $bed12
geneBody_coverage.py -i $bam -o $bam.rseqc -r $bed12
read_distribution.py -i $bam -r $bed12 > $bam.read_distribution.txt
read_duplication.py -i $bam -o $bam.read_duplication

################### preseq
echo Starting preseq
module unload rseqc
module load preseq
preseq lc_extrap -v -B $bam -o $bam.ccurve.txt

################## MarkDuplicates
module unload preseq
module load picard
java -Xmx2g -jar $PICARD_HOME/picard-1.118.jar MarkDuplicates \
    INPUT=${bam}_markduplicates \
    OUTPUT=${bam}_markduplicates.markDups.bam \
    METRICS_FILE=${bam}_markduplicates.markDups_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    PROGRAM_RECORD_ID=null \
    VALIDATION_STRINGENCY=LENIENT

################## Featurecounts
module unload picard
module load subread
featureCounts -a $gtf -g gene_id -o ${bam}_gene.featureCounts.txt -p -s 2 $bam
featureCounts -a $gtf -g gene_biotype -o ${bam}_biotype.featureCounts.txt -p -s 2 $bam
cut -f 1,7 ${bam}_biotype.featureCounts.txt > ${bam}_biotype_counts.txt

################# Stringtie
module unload subread
module load StringTie
stringtie $bam \
    -o ${bam}_transcripts.gtf \
    -v \
    -G $gtf \
    -A ${bam}.gene_abund.txt \
    -C ${bam}.cov_refs.gtf \
    -e \
    -b ${bam}_ballgown
