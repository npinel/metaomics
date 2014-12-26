#!/bin/bash

# load ESB-specific shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

# list and load modules
SCRIPTS='~/git/laoPipeline/'
MAPPER='bwa'
DATA=MT
QUAL=''

# Start workflow
#
# do in a while loo
# e.g.,
#while [  $COUNTER -lt 10 ]; do
#    let COUNTER=COUNTER+1
#done
YIELD=0
ITER=0
YIELD_CHANGE=100
TOTAL=100
MAPPED=50
while [[ $YIELD_CHANGE -ge 5 ]]; do
# 1. assemble
# ${SCRIPTS}04b.analysis.idba.sh -f <sample_list> -d <data_type>
#echo -e "${SCRIPTS}04b.analysis.idba.sh -f sample_list.iterative.txt -d MT -o -w"

# 2. merge new  assembly with pre-existing ones

# 3. map reads to assembly
#echo -e "${SCRIPTS}05a.mapping.${MAPPER}.sh"

# 4. count assembly yields (how many reads were mapped out of the total library)
# 5. extract unmapped reads from assembly, and prepare for re-assembly

# 6. repeat from 1.
let "ITER += 1"
NYIELD=$(( $TOTAL - $MAPPED ))
YIELD_CHANGE=$NYIELD
#YIELD_CHANGE=$(( $YIELD - $NYIELD ))

echo -e "iteration: $ITER\nmapped: $MAPPED\nyield: $YIELD"

let "MAPPED += 6"
#YIELD=$NYIELD
done
