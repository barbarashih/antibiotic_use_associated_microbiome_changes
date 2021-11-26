#!/bin/bash
# Grid Engine options
#$ -N repair
#$ -cwd
#$ -o log/repair.out
#$ -e log/repair.err
#$ -pe sharedmem 2
#$ -l h_vmem=30G
#$ -l h_rt=40:00:00
module add samtools/1.10

BBREFORMAT=${private_module_dir}/BBMap_38.71/reformat.sh
BBREPAIR=${private_module_dir}/BBMap_38.71/repair.sh

#### Organise data
# Organise variables
WORKING_DIR=$(pwd)
SAMPLE_LIST=${sample_list}
OUT_DIR=${out_dir}
IN_DIR=${fastq_trimmed_dir}

CURRENT_SAMPLE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
OUT_DIR=${OUT_DIR}/$CURRENT_SAMPLE
mkdir -p $(dirname $OUT_DIR)
mkdir -p $OUT_DIR

FQ1_PAIRED=$IN_DIR/${CURRENT_SAMPLE}_R1_trimmed_paired.fastq.gz
FQ2_PAIRED=$IN_DIR/${CURRENT_SAMPLE}_R2_trimmed_paired.fastq.gz

FQ1_REFORMATTED=$OUT_DIR/reformatted_1.fastq
FQ2_REFORMATTED=$OUT_DIR/reformatted_2.fastq
FQ1_REPAIRED=$OUT_DIR/repaired_1.fastq
FQ2_REPAIRED=$OUT_DIR/repaired_2.fastq
SINGLETON_REPAIRED=$OUT_DIR/singleton.fastq

# reformat the reads into a temp dir before running bwa
$BBREFORMAT -Xmx50g in1=$FQ1_PAIRED in2=$FQ2_PAIRED out1=$FQ1_REFORMATTED out2=$FQ2_REFORMATTED overwrite=t trimreaddescriptions=t addslash=t slashspace=f
$BBREPAIR -Xmx50g in1=$FQ1_REFORMATTED in2=$FQ2_REFORMATTED out1=$FQ1_REPAIRED out2=$FQ2_REPAIRED outs=$SINGLETON_REPAIRED overwrite=t 
rm $FQ1_REFORMATTED $FQ2_REFORMATTED
rm ${FQ1_REPAIRED}.gz ${FQ2_REPAIRED}.gz
gzip $FQ1_REPAIRED
gzip $FQ2_REPAIRED
