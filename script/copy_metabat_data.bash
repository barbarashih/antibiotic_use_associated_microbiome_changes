#!/bin/bash
# Grid Engine options
#$ -N copy_metabat
#$ -cwd
#$ -o log/copy_metabat.out
#$ -e log/copy_metabat.err
#$ -m e
#$ -pe sharedmem 1
#$ -l h_vmem=10G
#$ -l h_rt=40:00:0
# Initialise the modules framework

MAGPY_DIR=${magpy_dir}
METABAT_DIR=${metabat_dir}

# do it for all samples
for FOLDER in $METABAT_DIR/*/; do
	sample=$(basename $FOLDER)
	mkdir -p $MAGPY_DIR/copy_bin_tmp
	rm $MAGPY_DIR/copy_bin_tmp/*
	cp -r $FOLDER/bin/*.fa $MAGPY_DIR/copy_bin_tmp
	for file in $MAGPY_DIR/copy_bin_tmp/*.fa; do
		file_basename=$(basename $file)
		file_basename=${file_basename/.fa/}
		file_basename=${file_basename/./}
		TMP_FILE=$MAGPY_DIR/mags/$file_basename.tmp
		awk '/^>/{print ">" ++i; next}{print}' < $file > $TMP_FILE
		mv $TMP_FILE $MAGPY_DIR/mags/${file_basename}.fa
		sed -i "s/>/>${file_basename}_contig/g" $MAGPY_DIR/mags/${file_basename}.fa
	done
	rm $MAGPY_DIR/copy_bin_tmp/*
done
rm -d $MAGPY_DIR/copy_bin_tmp
