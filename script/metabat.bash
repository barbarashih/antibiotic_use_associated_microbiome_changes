#!/bin/bash
# Grid Engine options
#$ -N metabat
#$ -cwd
#$ -o log/metabat.out
#$ -e log/metabat.err
#$ -pe sharedmem 4
#$ -l h_vmem=64G
#$ -l h_rt=40:00:00
module unload python
module add bwa/0.7.17
module add samtools/1.10
module add anaconda/5.3.1
module add metabat/0.32.4


SAMPLE_LIST=${sample_list}
CURRENT_SAMPLE_NAME=$(awk "NR==$SGE_TASK_ID" $SAMPLE_LIST)
IN_DIR=${in_dir}
IN_FQ_DIR=${in_dir_fastq}
OUT_DIR=${out_dir}
IN_CONTIG_FILENAME=${in_contig_filename}
IN_DIR_BASE=$(basename $IN_DIR)

CURRENT_IN_DIR=$IN_DIR/$CURRENT_SAMPLE_NAME
CONTIG_FA=$CURRENT_IN_DIR/$IN_CONTIG_FILENAME

OUT_DIR=${OUT_DIR}/$IN_DIR_BASE/$CURRENT_SAMPLE_NAME
OUT_DIR_BAM=${OUT_DIR}/bam
OUT_DEPTH=$OUT_DIR/depth.txt
OUT_DIR_BIN=$OUT_DIR/bin/${CURRENT_SAMPLE_NAME}_bin

mkdir -p ${out_dir} $out_dir/$IN_DIR_BASE $OUT_DIR $(dirname $OUT_DIR_BIN)


# Run metabat
jgi_summarize_bam_contig_depths --minContigDepth 2 --minContigLength 2000 --outputDepth $OUT_DEPTH $OUT_DIR_BAM/*.bam
metabat2 -i $CONTIG_FA -a $OUT_DEPTH --verysensitive --minContig 2000 -t 4 --outFile $OUT_DIR_BIN


