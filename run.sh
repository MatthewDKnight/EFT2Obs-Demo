#!/usr/bin/env bash
source env.sh

set -x
set -e

if [[ $# -lt 1 ]]; then
    echo "Insufficient number of arguments, usage is ./setup_process.sh [name]]"
    exit 1
fi

PROCESS=$1
SAVEDIR=$2

### SET ENVIRONMENT VARIABLES HERE
RUNLABEL="pilotrun"
TMPDIR=/home/hep/mdk16/Masters/EFT2Obs-Demo2
###

cp cards/${PROCESS}/{param,reweight,run,pythia8}_card.dat ${MG_DIR}/${PROCESS}/Cards/

# Need to harcode the path to the FIFO, madgraph won't expand environment variables
sed -i "s@XTMPDIRX@${TMPDIR}@g" ${MG_DIR}/${PROCESS}/Cards/pythia8_card.dat

pushd ${MG_DIR}/${PROCESS}
# Create MG config
{
  echo "shower=Pythia8"
  echo "reweight=ON"
  echo "done"
} > mgrunscript

if [ -d "${MG_DIR}/${PROCESS}/Events/${RUNLABEL}" ]; then rm -r ${MG_DIR}/${PROCESS}/Events/${RUNLABEL}; fi
./bin/generate_events pilotrun < mgrunscript
popd

#rivet -v --analysis=HiggsTemplateCrossSectionsStage1 "${TMPDIR}/fifo.hepmc" -o Rivet.yoda
#rivet -v --analysis=SimpleHiggs "${TMPDIR}/fifo.hepmc" -o "${SAVEDIR}.yoda"
rivet -v --analysis=VH_bb "${TMPDIR}/fifo.hepmc" -o "${SAVEDIR}.yoda"
python save_info.py ${SAVEDIR} ${PROCESS}
#yoda2root -t Rivet.yoda
