#!/bin/bash
# Grid Engine options
#$ -N graphlan
#$ -cwd
#$ -o log/graphlan.out
#$ -e log/graphlan.err
#$ -pe sharedmem 1
#$ -l h_vmem=20G
#$ -l h_rt=40:00:00
module unload python
module add blast+/2.9.0
module add anaconda/5.3.1

source activate graphlan # install from https://anaconda.org/bioconda/graphlan

CURRENT_GRAPHLAN_DIR=${current_graphlan_dir}
mkdir -p $CURRENT_GRAPHLAN_DIR
CURRENT_TRE=$CURRENT_GRAPHLAN_DIR/input.tre
CURRENT_ANNOTATION_BASE=$CURRENT_GRAPHLAN_DIR/graphlan_annotation_base.txt
CURRENT_ANNOTATION_AMR=$CURRENT_GRAPHLAN_DIR/graphlan_annotation_amr.txt
CURRENT_XML_BASE=$CURRENT_GRAPHLAN_DIR/graphlan_base.xml
CURRENT_XML_AMR=$CURRENT_GRAPHLAN_DIR/graphlan_amr.xml
CURRENT_OUT_BASE=$CURRENT_GRAPHLAN_DIR/annotated_tree_base.pdf
CURRENT_OUT_AMR=$CURRENT_GRAPHLAN_DIR/annotated_tree_amr.pdf



graphlan_annotate.py $CURRENT_TRE $CURRENT_XML_BASE --annot $CURRENT_ANNOTATION_BASE
graphlan.py $CURRENT_XML_BASE $CURRENT_OUT_BASE  --dpi 150 --size 8
graphlan_annotate.py $CURRENT_TRE $CURRENT_XML_AMR --annot $CURRENT_ANNOTATION_AMR
graphlan.py $CURRENT_XML_AMR $CURRENT_OUT_AMR  --dpi 150 --size 10


source deactivate

