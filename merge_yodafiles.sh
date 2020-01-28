#!/usr/bin/env bash
cd /home/hep/mdk16/Masters/EFT2Obs-Demo2

source env.sh

set -x
set -e

SAVEDIR=$1
NORUNS=$2

OUTPUT=${SAVEDIR}.yoda

INPUT=

### Loop over number of runs
for ((i=0 ; i<$NORUNS ; i++)); do
	INPUT=$INPUT$" "tmp_yodafiles/${i}.yoda
done

yodamerge -o ${OUTPUT} ${INPUT}
cp tmp_yodafiles/0.txt ${SAVEDIR}.txt
rm tmp_yodafiles/* 
