#!/bin/bash
# Grid Engine options
#$ -N phylophlan_all
#$ -cwd
#$ -o log/phylophlan_all.out
#$ -e log/phylophlan_all.err
#$ -m e
#$ -pe sharedmem 2
#$ -l h_vmem=32G
#$ -l h_rt=60:00:00
. /etc/profile.d/modules.sh
module unload python
module add blast+/2.9.0
module add anaconda/5.3.1


# Organise variables
WORKING_DIR=$(pwd)
MAGPY_DIR=${magpy_dir}
PHYLOPHLAN_DIR=${phylophlan_dir}

cd $MAGPY_DIR
mkdir -p $PHYLOPHLAN_DIR $PHYLOPHLAN_DIR/input $PHYLOPHLAN_DIR/output
PHYLOPHLAN_IN_DIR=${PHYLOPHLAN_DIR}/input
PHYLOPHLAN_OUT_DIR=${PHYLOPHLAN_DIR}/output
PHYLOPHLAN_LOG_DIR=${PHYLOPHLAN_DIR}/logs
mkdir -p $(dirname $PHYLOPHLAN_OUT_DIR)
mkdir -p $PHYLOPHLAN_IN_DIR $PHYLOPHLAN_LOG_DIR $PHYLOPHLAN_OUT_DIR


source activate phylophlan # install phylophlan https://huttenhower.sph.harvard.edu/phylophlan/#:~:text=PhyloPhlAn%20is%20an%20integrated%20pipeline,at%20multiple%20levels%20of%20resolution.


phylophlan \
    -i $PHYLOPHLAN_IN_DIR \
    -o $PHYLOPHLAN_OUT_DIR \
    -d phylophlan \
    -t a \
    --accurate \
    --nproc 4 \
    --min_num_entries 4 \
    -f isolates_config.cfg \
    --diversity high \
    --verbose 2>&1 | tee $PHYLOPHLAN_LOG_DIR/all.log


source deactivate
