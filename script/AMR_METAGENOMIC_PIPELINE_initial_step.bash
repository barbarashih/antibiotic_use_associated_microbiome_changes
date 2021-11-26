#!/bin/bash

############### CREATED BY ADRIAN MUWONGE 10/08/2017###################################################################################################################################################################
#######################################################################################################################################################################################################################

#THIS PIPELINE WAS ASSEMBLED TO ACHIEVE THE FOLLOWING OUTPUTS
#TO QUALITY CHECK METAGENOMIC READS
#TO SCREEN FOR PRESENCE OF AMR GENES(CLASSES,MODE OF ACTION AND SITE OF ACTIVITY)
#TO SCREEN FOR SNP LEVEL CHANGES  THAT HAVE OCCURED IN THE AMR GENES-POSSIBLY USE THEM TO INFER PHYLOGENY
#TO OUT PUT GENE AND PATHWAY ABUNDANCES VIA THE HUMANN2 PIPELINE
#TO OUTPUT TAXONOMIC ABUNDANCES
#TO OUTPUT GENEFAMILY ABUNDANCES AND COVERAGE
#TO PUTPUT PATHWAY ABUNDANCES AND COVERAGE

#THESE OUTPUTS CAN THEN FEED INTO STATISTICAL, EPIDEMIOLOGICAL AND GRAPHICAL ANALYSIS


echo "*******CHECKING AND LOADING REQUIRED MODULES******"

module add bowtie2/2.2.6

READS_DIR="READS_YORK" # the reads are placed in here

### We start by removing the contaminants, Phix and Eukaryotic DNA by mapping to the pig an indexed pig and phix reference###

echo "*******CREATING BOWTIE PARAMETERS******"


### Make sure to only give the prefix name of the indexed file for option - - db
echo "--very-sensitive-local" >> ./bowtie2_config.txt

echo "*******RUNNING YOUR SCREENING FOR THE PIG GENOME NOW***********"

script/run_contaminant_filter.pl $READS_YORK/*/*.fastq -p 30 --db GCA_000003025 -o screened_reads/ -c ./bowtie2_config.txt

echo "*******RUNNING YOUR SCREENING FOR THE PHIX174 GENOME NOW***********"

mkdir Step_2

cd Step_2

script/run_contaminant_filter.pl ../screened_reads/*.fastq -p 30 --db hiX174 -o screened_reads_2/ -c ../bowtie2_config.txt


cd ../


echo "*******RENAME FILES TO REMOVE UNDERSCORES***********"

### THE UNDERSCORES INTERFERE WITH THE TRIMOMATIC CODE

for i in Step_2/screened_reads_2/*_screened_screened.fastq; do mv $i ${i%_screened_screened.fastq}.fastq; done

cd Step_2/screened_reads_2/

for i in *.fastq; do target=$(echo $i | sed -e "s/_//"); mv "$i" "$target"; done

for i in *.fastq; do target=$(echo $i | sed -e "s/_//"); mv "$i" "$target"; done

cd ../../

echo "*******CLEAN UP STEP ***********"

rm -r screened_reads/

echo "*******RUNNING TRIMOMMATIC***********"

for i in Step_2/screened_reads_2/*_R1.fastq;

do
# please download trimmomatic from http://www.usadellab.org/cms/?page=trimmomatic 
run_trimmomatic.pl -l 5 -t 15 -r 15 -w 4 -m 70 -j trimmomatic/0.36/trimmomatic-0.36.jar --thread 30 -o trimmomatic_filtered  Step_2/screened_reads_2/*.fastq;

done

echo "*******CLEAN UP STEP ***********"

rm -r Step_2/

echo "*******STARTING THE AMR GENE SCREENING  PIPELINE***********"




