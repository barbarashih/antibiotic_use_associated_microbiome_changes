#!/bin/bash
#$ -N abricate
#$ -cwd
#$ -o log/abricate.out
#$ -e log/abricate.err
#$ -l h_rt=5:00:00
#$ -l h_vmem=40G
#$ -pe sharedmem 2
module load anaconda/5.0.1


WORKING_DIR=$(pwd)
SAMPLE_LIST=$WORKING_DIR/${passed_bins}
MAGPY_DIR=${magpy_dir}
OUT_DIR=${out_dir}
USE_DB=${use_db}

mkdir -p $OUT_DIR

CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
IN_FP=${MAGPY_DIR}/mags/${CURRENT_SAMPLE}.fa
mkdir -p $OUT_DIR/$USE_DB
OUT_FP=$OUT_DIR/$USE_DB/$CURRENT_SAMPLE.csv



if [[ ! -s $OUT_FP ]]; then
source activate abricate # install from https://anaconda.org/bioconda/abricate

abricate $IN_FP --db $USE_DB --threads 2 --csv > $OUT_FP


source deactivate
fi
