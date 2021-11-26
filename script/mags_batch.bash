#!/bin/bash
# Grid Engine options
#$ -N magpy_batch_mags
#$ -cwd
#$ -o log/magpy_batch_mags.out
#$ -e log/magpy_batch_mags.err
#$ -pe sharedmem 1
#$ -l h_vmem=5G
#$ -l h_rt=40:00:00

# split mags into batches
MAGPY_DIR=${magpy_dir}
WORKING_DIR=$(pwd)
cd $MAGPY_DIR

# split it into batches of 1000 genomes per batch
MAGS_SPLIT_DIR=mags_split
MAG_LIST=$WORKING_DIR/${mag_list}
i=1
j=1
k=1
batch_size=${batch_size}
lineNum=$(wc -l < $MAG_LIST)
while [ "$lineNum" -gt "$i" ]; do
if [ "$k" -gt "$batch_size" ]; then
k=1
(( j++ ))
fi
current_sample=$(awk "NR==$i" $MAG_LIST)
current_batch="${MAGS_SPLIT_DIR}/mags_${batch_size}_${j}"
mkdir -p $current_batch
cp mags/$current_sample $current_batch
(( i++ ))
(( k++ ))
done

