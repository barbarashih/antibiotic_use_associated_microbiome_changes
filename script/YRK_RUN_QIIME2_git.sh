#Step 1 ### Read in the data

    qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
   --input-path YORK_MANIFEST_select2.csv\
   --output-path YORK_MANIFEST_select2.qza\
    --input-format PairedEndFastqManifestPhred33
    
    ### Visualise demultiplexed data
    
     
    qiime demux summarize --i-data YORK_MANIFEST_select2.qza --o-visualization YORK_MANIFEST_select.qzv
    
  #****************************************************  ***************************************************************************************************************************************************
   
#Step 2 DENOISING(supply the original out from step 1 dada does its own QC)


qiime dada2 denoise-paired\
  --i-demultiplexed-seqs YORK_MANIFEST_select.qza \
  --o-table YORK_SELECT_dadatable.qza\
  --output-dir denoiseS_select_output\
  --p-trim-left-f 10\
  --p-trim-left-r 10\
  --p-trunc-len-f 120\
  --p-trunc-len-r 100\
  --p-n-threads 30\
  --o-representative-sequences rep-seqs.qza\
  --o-denoising-stats denoisingS-stats.qza\
  --p-n-reads-learn 200000
  
  
  qiime metadata tabulate \
  --m-input-file denoise_select_output/denoisingS-stats.qza \
  --o-visualization stats-dada2.qzv
  
  #Step 3 add metadata to qiime artifact

time qiime feature-table summarize\
  --i-table denoise_select_output/YORK_SELECT_dadatable.qza\
  --o-visualization YORK_SELECT_dadatable.qzv\
  --m-sample-metadata-file MetadataYork_shtgun_samples.txt
  #******************************************************************************************************************************************************
  #Step 4 Alignment, phylogenetic trees and Rarefaction curves
  
  time qiime alignment mafft\
  --i-sequences denoise_select_output/rep-seqs.qza\
  --p-n-threads 20\
  --output-dir alignment_output\
  --o-alignment selectYRK-rep-seqs.qza
  
  time qiime alignment mask\
  --i-alignment selectYRK-rep-seqs.qza\
  --o-masked-alignment selectYRK_masked-alignedrep-seqs.qza 
  
  
  time qiime phylogeny fasttree\
  --i-alignment selectYRK_masked-alignedrep-seqs.qza \
  --o-tree selectYRK_unrooted-tree.qza
  
  Midpoint rooting
  time qiime phylogeny midpoint-root\
  --i-tree selectYRK_unrooted-tree.qza\
  --o-rooted-tree selectYRK_rooted-tree.qza
  
  qiime diversity alpha-rarefaction \
--i-table denoise_select_output/YORK_SELECT_dadatable.qza \
--i-phylogeny selectYRK_rooted-tree.qza \
--p-max-depth 10000 \
--m-metadata-file MetadataYork_shtgun_samples.txt \
--p-metrics observed_otus \
--p-metrics shannon \
--p-metrics faith_pd \
--o-visualization rarefaction_10000.qzv

#*******************************************************************************************************************************************************************************************************

	
  
## Diversity indicies (Alpha and Beta)  
  
qiime diversity core-metrics-phylogenetic --i-phylogeny selectYRK_rooted-tree.qza --i-table denoise_select_output/YORK_SELECT_dadatable.qza --p-sampling-depth 10000 --m-metadata-file MetadataYork_shtgun_samples.txt --output-dir metrics

qiime diversity alpha-group-significance --i-alpha-diversity metrics/faith_pd_vector.qza --m-metadata-file MetadataYork_shtgun_samples.txt --o-visualization metrics/faith-pd-group-significance.qzv

qiime diversity alpha-correlation --i-alpha-diversity metrics/evenness_vector.qza --m-metadata-file MetadataYork_shtgun_samples.txt --o-visualization metrics/evenness-alpha-correlation.qzv


qiime emperor plot --i-pcoa metrics/unweighted_unifrac_pcoa_results.qza --m-metadata-file MetadataYork_shtgun_samples.txt --o-visualization metrics/unweighted-unifrac-emperor.qzv

#******************************************************************************************************************************************************

# Taxonomic assignment
# Here we are using 99 and consensus all taxonomic levels
# please download the SILVA_TAXONOMIC_DB for the following commands
qiime tools import --type 'FeatureData[Sequence]' --input-path SILVA_TAXONOMIC_DB/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna   --output-path silva_132_99_16S.qza

qiime tools import --type 'FeatureData[Taxonomy]'  --input-format HeaderlessTSVTaxonomyFormat --input-path SILVA_TAXONOMIC_DB/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_all_levels.txt --output-path ref-99taxonomy.qza


qiime feature-classifier extract-reads --i-sequences silva_132_99_16S.qza --p-f-primer CCTACGGGAGGCAGCAG --p-r-primer ATTACCGCGGCTGCTGG --p-trunc-len 150 --o-reads ref-seqs.qza

#Train Classifier

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-99seqs.qza --i-reference-taxonomy ref-99taxonomy.qza --o-classifier classifier_99.qza


#USING 99

qiime feature-classifier extract-reads --i-sequences silva_132_99_16S.qza --p-f-primer CCTACGGGAGGCAGCAG --p-r-primer ATTACCGCGGCTGCTGG --p-trunc-len 220 --o-reads ref-99seqs.qza

Train Classifier

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-99seqs.qza --i-reference-taxonomy ref-99taxonomy.qza --o-classifier classifier99.qza

qiime feature-classifier classify-sklearn --i-classifier classifier_99.qza --i-reads denoise_select_output/rep-seqs.qza  --o-classification Classify_ref-99_taxonomy.qza

qiime metadata tabulate --m-input-file Classify_ref-99_taxonomy.qza --o-visualization Classify_ref-99_taxonomy.qzv
  
  
 qiime taxa barplot --i-table denoise_select_output/YORK_SELECT_dadatable.qza --i-taxonomy Classify_ref-99_taxonomy.qza --m-metadata-file MetadataYork_shtgun_samples.txt --o-visualization YRK_classify_taxa-bar-plots.qzv
#****************************************************************************************************************************************************************************************************

  
## REMOVE CONTAMINANTS(mitochondria and chloroplast)
  
  qiime tools export --output-path taxonomy-export \
  --input-path Classify_ref-99_taxonomy.qza

grep -v -i "mitochondia|chloroplast|Feature" taxonomy-export/taxonomy.tsv | cut  -f 1 > no-chloro-mito-ids.txt

#****************************************************************************************************************************************************************************************************


### Files for use in R for Phyloseq artifact
# Export data to biom format
qiime tools export --output-path dada2-table-export --input-path denoise_select_output/YORK_SELECT_dadatable.qza 
# Move into the directory
cd dada2-table-export
#****************************************************************************************************************************************************************************************************

 

# Convert the HDF5 biom file to a tsv biom file
biom subset-table \
  --input-hdf5-fp feature-table.biom \
  --axis observation \
  --ids ../no-chloro-mito-ids.txt \
  --output-fp feature-table-subset.biom
  
  biom convert -i table.biom -o table.from_biom.txt --to-tsv
  
#****************************************************************************************************************************************************************************************************
 
  

biom convert  --input-fp collapse.frequency/feature-table.biom --output-fp YRK_collapse.frequency.table.txt --sample-metadata-fp MetadataYork_shtgun_samples.txt --header-key “pig.type” --to-tsv
