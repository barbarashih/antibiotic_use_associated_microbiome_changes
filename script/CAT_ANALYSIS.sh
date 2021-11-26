#!/bin/bash
#$ -N PLASMID_TAXA_GENE_SCREEN3
# Hard runtime limit
#$ -l h_rt=99:00:00
#$ -hold_jid Stagein_05043
#$ -l h_vmem=16G
#$ -pe sharedmem 8


module load anaconda/5.0.1


source activate DIAMOND # please install cat in this enviroment https://github.com/dutilh/CAT
CATBAT="" # please download the cat database in this direcotry
TRIMMED_READ_DIR="" # please add in the directory for trimmed reads here
cd $TRIMMED_READ_DIR



for i in *.fasta;

do

CAT contigs -c $i -d $CATBAT/2019-01-08_CAT_database/ -t $CATBAT/2019-01-08_taxonomy/ -n 8 --out_prefix $i.CAT_run

CAT add_names -i $i.CAT_run.contig2classification.txt -o $i.CAT_run.contig2classification.official_names.txt -t $CATBAT/2019-01-08_taxonomy/

CAT summarise -c $i -i $i.CAT_run.contig2classification.official_names.txt -o $i_run.summary.txt

rm $i.CAT_run.predicted_proteins.gff
rm $i.CAT_run.predicted_proteins.faa
rm $i.CAT_run.ORF2LCA.txt
rm $i.CAT_run.alignment.diamond
done

source deactivate
