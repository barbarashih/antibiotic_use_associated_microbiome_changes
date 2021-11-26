#!/bin/bash
# Grid Engine options
#$ -N prodigal
#$ -cwd
#$ -o log/prodigal.out
#$ -e log/prodigal.err
#$ -pe sharedmem 1
#$ -l h_vmem=10G
#$ -l h_rt=40:00:00
module unload python
module add anaconda/5.3.1


# Organise variables
WORKING_DIR=$(pwd)
MAGPY_DIR=${magpy_dir}
SAMPLE_LIST=${sample_list}

CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
cd $MAGPY_DIR

CURRENT_SAMPLE=$(basename $CURRENT_SAMPLE)
IN_FA=mags/$CURRENT_SAMPLE.fa
OUT_FAA=proteins/${CURRENT_SAMPLE}.faa
OUT_GFF=proteins/${CURRENT_SAMPLE}_prodigal.gff
mkdir -p proteins


source activate prodigal # install prodigal https://anaconda.org/bioconda/prodigal

prodigal -p meta -a $OUT_FAA -q -i $IN_FA -f gff -o $OUT_GFF

conda deactivate


