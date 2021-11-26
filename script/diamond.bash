#!/bin/bash
# Grid Engine options
#$ -N diamond
#$ -cwd
#$ -o log/diamond.out
#$ -e log/diamond.err
#$ -pe sharedmem 2
#$ -l h_vmem=40G
#$ -l h_rt=5:00:00
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
OUT_TSV=diamond/${CURRENT_SAMPLE}.diamond.tsv
mkdir -p diamond


if [[ ! -f $OUT_TSV ]]; then

source activate magpy_install # please install magpy from https://github.com/WatsonLab/MAGpy

echo $OUT_TSV
UNIPORT_DB=uniprot_trembl.dmnd # please download the uniprot_trembl database here. Downloaded on 11/Dec/2020
diamond blastp --threads 2 --max-target-seqs 10 --db $UNIPORT_DB --query $IN_FAA --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out $OUT_TSV

conda deactivate

fi

