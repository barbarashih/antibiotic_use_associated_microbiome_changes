#!/bin/bash
# Grid Engine options
#$ -N checkm
#$ -cwd
#$ -o log/checkm.out
#$ -e log/checkm.err
#$ -pe sharedmem 2
#$ -l h_vmem=30G
#$ -l h_rt=100:00:00

module add blast+/2.9.0
module add anaconda/5.3.1


# Organise variables
WORKING_DIR=$(pwd)
MAGPY_DIR=${magpy_dir}
CHECKM_DB=checkm_data
MAGS_SPLIT=mags_split
BATCH_SIZE=${batch_size}
CHECKM_OUT_FILE=checkm_batch/checkm_${BATCH_SIZE}_${SGE_TASK_ID}.txt
CHECKM_OUT_DIR=checkm_batch/checkm_${BATCH_SIZE}_${SGE_TASK_ID}
CHECKM_OUT_TMP_DIR=checkm_batch/checkm_${BATCH_SIZE}_${SGE_TASK_ID}_tmp
CURRENT_BATCH_DIR="$MAGS_SPLIT/mags_${BATCH_SIZE}_${SGE_TASK_ID}"



cd $MAGPY_DIR

mkdir -p checkm_batch $CHECKM_OUT_TMP_DIR

source activate magpy_install # install from https://github.com/WatsonLab/MAGpy

checkm data setRoot ${CHECKM_DB}
checkm lineage_wf -f $CHECKM_OUT_FILE --reduced_tree -t 4 -x fa $CURRENT_BATCH_DIR ./$CHECKM_OUT_DIR --tmpdir $CHECKM_OUT_TMP_DIR
scripts/add_tax.py $CHECKM_OUT_FILE > ${CHECKM_OUT_FILE/.txt/_plus.txt}

source deactivate
