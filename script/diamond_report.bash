#!/bin/bash
# Grid Engine options
#$ -N diamond_report
#$ -cwd
#$ -o log/diamond_report.out
#$ -e log/diamond_report.err
#$ -pe sharedmem 1
#$ -l h_vmem=20G
#$ -l h_rt=20:00:00
module unload python
module add anaconda/5.3.1
module add blast+/2.9.0


# Organise variables
WORKING_DIR=$(pwd)
MAGPY_DIR=${magpy_dir}
SAMPLE_LIST=${sample_list}

CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
cd $MAGPY_DIR

CURRENT_SAMPLE=$(basename $CURRENT_SAMPLE)
IN_FAA=proteins/${CURRENT_SAMPLE}.faa
IN_TSV=diamond/${CURRENT_SAMPLE}.diamond.tsv
OUT_TSV=diamond_report/${CURRENT_SAMPLE}.tsv
mkdir -p diamond_report

if [[ ! -s $OUT_TSV ]]; then
source activate bioperl # please install from https://anaconda.org/bioconda/perl-bioperl

scripts/diamond_report.pl $IN_TSV $IN_FAA "diamond_report"


conda deactivate

fi

