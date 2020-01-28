#!/usr/bin/env bash

set -x
set -e

PROCESS=$1
SAVEDIR=$2
NORUNS=$3

### SET ENVIRONMENT VARIABLES HERE
RUNLABEL="pilotrun"
###

### Loop over number of runs
for ((i=0 ; i<$NORUNS ; i++)); do
	qsub -q hep.q -l h_rt=600 run_batch.sh ${PROCESS} tmp_yodafiles/${i}
done
qsub -q hep.q -l h_rt=60 -hold_jid run_batch.sh merge_yodafiles.sh ${SAVEDIR} ${NORUNS}


