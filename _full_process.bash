#### Pre-requisite
# QC and trimming on reads were carried out by script/AMR_METAGENOMIC_PIPELINE_initial_step.bash

#### SETTINGs
## directories (need to be set)
WORKING_DIR="" # please put in your working path
SCRATCH_DIR="tmp"
PRIVATE_MODULE_DIR="private_modules" # please download and unzip BBMap_38.71 in this directory
FASTQ_ENDING="fastq.gz"
FASTQ_TRIMMED_DIR="data/trimmed_fastq" # please put the trimmed fastq from running script/AMR_METAGENOMIC_PIPELINE_initial_step.bash in this directory
FASTQ_REPAIRED_DIR="data/fastq_trimmed_repair"
DATA_DIR="data"

# Generated file paths/dirs
ANALYSIS_DIR="analysis"
SCRIPT_DIR="script"
FILE_LIST_DIR="file_list"
MEGAHIT_DIR="${ANALYSIS_DIR}/megahit/"
METAPLASMID_DIR="${ANALYSIS_DIR}/metaplasmid"
METABAT_DIR="${ANALYSIS_DIR}/metabat"
MAGPY_DIR="${ANALYSIS_DIR}/MAGpy/" # Download MAGpy here
ABRICATE_DIR="${ANALYSIS_DIR}/abricate/"
PLASMIDFINDER_DIR="${ANALYSIS_DIR}/plasmidfinder/"
GRAPHLAN_DIR="${ANALYSIS_DIR}/graphlan"
MAGPY_MAGS="${MAGPY_DIR}/mags"
MAGPY_MAGS_SPLIT="${MAGPY_DIR}/mags_split"
MAGPY_PHYLOPHLAN="${MAGPY_DIR}/phylophlan"

SAMPLE_LIST=$FILE_LIST_DIR/sample_list.txt # list the samples in this file
CHECKM_PASSED_BINS=$FILE_LIST_DIR/checkm_passed_bins.txt
CHECKM_BIN_ANNOTATION=$FILE_LIST_DIR/checkm_bin_annotation.txt
ALL_BINS=$FILE_LIST_DIR/all_bins.txt
MEGARESDB="${DATA_DIR}/MegaresDB/megares_database_v1.01.fasta" # download the megares database here
MEGARESDB_ANNOTATION="${DATA_DIR}/MegaresDB/megares_annotations_v1.01.csv" # download the annotation for megares database here
SAMPLE_NUM=27

mkdir -p $WORKING_DIR
cd $WORKING_DIR

module add python/3.4.3


mkdir -p $FILE_LIST_DIR $ANALYSIS_DIR $MAGPY_DIR $ABRICATE_DIR $GRAPHLAN_DIR $MAGPY_MAGS $MAGPY_PHYLOPHLAN $MAGPY_MAGS_SPLIT \
$MEGAHIT_DIR $METABAT_DIR $METABAT_DIR $FASTQ_REPAIRED_DIR $PLASMIDFINDER_DIR


# generate sample list
cp data/YORK_METAGENOME_READS/IDS.txt $SAMPLE_LIST # a list of sample names

#### repair 
qsub -t 1-27 -v private_module_dir=$PRIVATE_MODULE_DIR,out_dir=$FASTQ_REPAIRED_DIR,fastq_trimmed_dir=$FASTQ_TRIMMED_DIR,sample_list=$SAMPLE_LIST ${SCRIPT_DIR}/bbmap_repair.bash


#### Run megahit
qsub -t 1-27 -l h_vmem=32G -pe sharedmem 4 -v in_h_vmem_G=100,in_core=4,out_dir=$MEGAHIT_DIR,in_fastq_dir=$FASTQ_TRIMMED_DIR,sample_list=$SAMPLE_LIST ${SCRIPT_DIR}/megahit.bash
# summarise contig length

#### Run contig analysis
# use the script script/CAT_ANALYSIS.sh to analyse contigs 

#### Bin the contigs into genomes ---------------------------------
#### Run metabat
qsub -t 1-27 -hold_jid "metabat_bwa" -v in_dir=$MEGAHIT_DIR,in_contig_filename="final.contigs.fa",out_dir=$METABAT_DIR,sample_list=$SAMPLE_LIST ${SCRIPT_DIR}/metabat.bash



#### MAGpy
# because I cannot run snakemake continuously on the HPC setup, the steps in MAGpy were extracted and run individually 
# after all the files are generated, MAGpy was run to ensure all output files are present and ready for analysis
# checkM
qsub -v magpy_dir=$MAGPY_DIR,metabat_dir=$METABAT_DIR/megahit ${SCRIPT_DIR}/copy_metabat_data.bash 
ls $MAGPY_DIR/mags > $ALL_BINS

# Run checkM in batches
qsub -hold_jid copy_metabat -v magpy_dir=$MAGPY_DIR,batch_size=250,mag_list=$ALL_BINS ${SCRIPT_DIR}/mags_batch.bash
qsub -hold_jid magpy_batch_mags -t 1-13 -v magpy_dir=$MAGPY_DIR,batch_size=250 ${SCRIPT_DIR}/checkM_array.bash


# prodigal
qsub -t 1-$(sed -n '$=' $CHECKM_PASSED_BINS) -v magpy_dir=$MAGPY_DIR,sample_list=$CHECKM_PASSED_BINS ${SCRIPT_DIR}/prodigal.bash

# manual diamond
qsub -hold_jid prodigal -t 1-$(sed -n '$=' $CHECKM_PASSED_BINS) -v magpy_dir=$MAGPY_DIR,sample_list=$CHECKM_PASSED_BINS ${SCRIPT_DIR}/diamond.bash
qsub -hold_jid diamond -t 1-$(sed -n '$=' $CHECKM_PASSED_BINS) -v magpy_dir=$MAGPY_DIR,sample_list=$CHECKM_PASSED_BINS ${SCRIPT_DIR}/diamond_report.bash


# combine diamond report
module unload python
module add anaconda/5.3.1
source activate magpy_install # install MAGpy in this conda enviroment
cd $MAGPY_DIR
rm diamond_bin_report.tsv "diamond_bin_report_plus.tsv"
echo -e 'name\tnprots\tnhits\tnfull\tgenus\tngenus\tspecies\tnspecies\tavgpid' >> diamond_bin_report.tsv
find diamond_report/ -name "bin*.tsv" | xargs -I {{}} cat {{}} >> diamond_bin_report.tsv
scripts/add_tax_diamond.py "diamond_bin_report.tsv" > "diamond_bin_report_plus.tsv" # this script is downloaded from MAGpy
source deactivate
cd -


# run this in wildwest to ensure all prodigal and diamond files have been generated correctly
#cd $MAGPY_DIR
#module unload python
#module add anaconda/5.3.1
#conda activate magpy_install
#snakemake -s MAGpy_2 --cores 2
#cd $WORKING_DIR


#### End of MAGpy	

# Run all samples together (generate the final tree)
# This is done by moving all the tmp files from all samples into a single folder
mkdir -p $MAGPY_DIR/phylophlan_all $MAGPY_DIR/phylophlan_all/input
cp $MAGPY_DIR/proteins/*.faa $MAGPY_DIR/phylophlan_all/input

#### Run abricate
qsub -t 1-$(sed -n '$=' $CHECKM_PASSED_BINS) -v passed_bins=$CHECKM_PASSED_BINS,magpy_dir=$MAGPY_DIR,out_dir=$ABRICATE_DIR,use_db=MegaresDB ${SCRIPT_DIR}/abricate.bash
qsub -v abricate_dir=$ABRICATE_DIR ${SCRIPT_DIR}/abricate.bash

# run abricate plasmid finder on all metaplasmid results
qsub -t 1-$(sed -n '$=' $SAMPLE_LIST) -v sample_list=$SAMPLE_LIST,in_dir=$MEGAHIT_DIR,out_dir=$ABRICATE_DIR,use_db=MegaresDB,in_filename="final.contigs.fa" ${SCRIPT_DIR}/abricate_prebinning.bash

# run metaplasmid results against abricate MEGARESDB
qsub -t 1-$(sed -n '$=' $SAMPLE_LIST) -v sample_list=$SAMPLE_LIST,in_dir=$METAPLASMID_DIR,out_dir=$ABRICATE_DIR,use_db=MegaresDB,in_filename="contigs.fasta" ${SCRIPT_DIR}/abricate_prebinning.bash

# run megahit results against abricate MEGARESDB
qsub -t 1-$(sed -n '$=' $SAMPLE_LIST) -v sample_list=$SAMPLE_LIST,in_dir=$MEGAHIT_DIR,out_dir=$ABRICATE_DIR,use_db=MegaresDB,in_filename="final.contigs.fa" ${SCRIPT_DIR}/abricate_prebinning.bash



#### plot trees
# for all bins 
python script/sample_annotation_reorganise.py
CURRENT_PHYLOPHLAN_DIR=phylophlan_all
CURRENT_CHECKM_THRESHOLD=checkm75_15
CURRENT_PARAMETER=${CURRENT_PHYLOPHLAN_DIR}_${CURRENT_CHECKM_THRESHOLD}
CURRENT_GRAPHLAN_DIR=$GRAPHLAN_DIR/$CURRENT_PARAMETER/
mkdir -p $CURRENT_GRAPHLAN_DIR
# copy phylophlan output to graphlan results folder for plotting trees
cp $MAGPY_DIR/$CURRENT_PHYLOPHLAN_DIR/output/*.tre $CURRENT_GRAPHLAN_DIR
# create annotation files for graphlan
python ${SCRIPT_DIR}/graphlan_annotation.py "${CHECKM_PASSED_BINS}" $MAGPY_DIR/diamond_bin_report_plus.tsv $ABRICATE_DIR/plasmidfinder $ABRICATE_DIR/MegaresDB $MEGARESDB_ANNOTATION no $CHECKM_BIN_ANNOTATION $MAGPY_DIR/$CURRENT_PHYLOPHLAN_DIR/output/input_concatenated.aln $CURRENT_GRAPHLAN_DIR/graphlan
# run graphlan
qsub -v current_graphlan_dir=$CURRENT_GRAPHLAN_DIR ${SCRIPT_DIR}/graphlan.bash


#### Other summaries

# get an overview for the classified bins
python script/diamond_plus_summarise.py file_list/checkm_passed_bins.txt ${MAGPY_DIR}/diamond_bin_report_plus.tsv ${MAGPY_DIR}

# Get an overview on the contigs
python script/megahit_contig_summary.py