#!/usr/bin/env bash
cd /home/hep/mdk16/Masters/EFT2Obs-Demo2

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
###

### Make pipe
mkfifo ${TMPDIR}/fifo.hepmc

cp -r ${MG_DIR}/${PROCESS} ${TMPDIR}/${PROCESS}

cp cards/${PROCESS}/{param,reweight,run,pythia8}_card.dat ${TMPDIR}/${PROCESS}/Cards/

# Need to harcode the path to the FIFO, madgraph won't expand environment variables
sed -i "s@XTMPDIRX@${TMPDIR}@g" ${TMPDIR}/${PROCESS}/Cards/pythia8_card.dat

pushd ${TMPDIR}/${PROCESS}
# Create MG config
{
  echo "shower=Pythia8"
  echo "reweight=ON"
  echo "done"
} > mgrunscript

if [ -d "${TMPDIR}/${PROCESS}/Events/${RUNLABEL}" ]; then rm -r ${TMPDIR}/${PROCESS}/Events/${RUNLABEL}; fi
./bin/generate_events pilotrun < mgrunscript
popd

#rivet -v --analysis=HiggsTemplateCrossSectionsStage1 "${TMPDIR}/fifo.hepmc" -o Rivet.yoda
#rivet -v --analysis=SimpleHiggs "${TMPDIR}/fifo.hepmc" -o "${SAVEDIR}.yoda"
rivet -v --analysis=VH ${TMPDIR}/"fifo.hepmc" -o ${TMPDIR}/"results.yoda"
cp ${TMPDIR}/results.yoda ${SAVEDIR}.yoda
python save_info.py ${SAVEDIR} ${PROCESS}
#yoda2root -t Rivet.yoda
