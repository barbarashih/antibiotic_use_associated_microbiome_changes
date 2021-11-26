#!/bin/bash
# Grid Engine options
#$ -N megahit
#$ -cwd
#$ -o log/megahit.out
#$ -e log/megahit.err
#$ -pe sharedmem 4
#$ -l h_vmem=30G
#$ -l h_rt=60:00:00
module add bwa/0.7.17
module add samtools/1.10
module add anaconda/5.3.1
module add megahit/1.1.3

#source activate megahit_metabat

IN_DIR=${in_fastq_dir}
OUT_DIR=${out_dir}
SAMPLE_LIST=${sample_list}

IN_MEM=${in_h_vmem_G}
IN_CORE=${in_core}
#IN_MEM=$(( $IN_MEM * 1000000000 - 3000000000 )) # convert to byte

CURRENT_SAMPLE_LINE=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
CURRENT_SAMPLE_NAME=$CURRENT_SAMPLE_LINE
FQ1=${IN_DIR}/${CURRENT_SAMPLE_NAME}_R1_trimmed_paired.fastq.gz
FQ2=${IN_DIR}/${CURRENT_SAMPLE_NAME}_R2_trimmed_paired.fastq.gz 
OUT_DIR=${OUT_DIR}/${CURRENT_SAMPLE_NAME}
OUT_CONTIGS=$OUT_DIR

mkdir -p $(dirname $OUT_DIR)
#mkdir -p $OUT_DIR 
#rm -r $OUT_CONTIGS


# make contigs with megahit
# -m option is in bytes. 21474836480 is 20G
megahit -1 $FQ1 -2 $FQ2 --continue --kmin-1pass --presets meta-sensitive --min-contig-len 500 -t $IN_CORE -o $OUT_CONTIGS

#source deactivate