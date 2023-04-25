#!/bin/bash

#
#SBATCH -D /lustre/home/user/v/vtroshin/output/MCTests/7bins
#SBATCH -J MCTest
#SBATCH -p mephi
#SBATCH -t 1-4
#SBATCH -a 1-2499
#SBATCH -x n02p002,n02p003,n02p004,n02p005
#
#SBATCH -o /lustre/home/user/v/vtroshin/output/MCTests/7bins/logs/slurm_%A_%a.out
#SBATCH -e /lustre/home/user/v/vtroshin/output/MCTests/7bins/logs/slurm_%A_%a.err
#

date
hostname

source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev
export MPDROOT=/lustre/home/user/n/nazarova/mpd
source $MPDROOT/config/env.sh

module add GVR/v1.0-1

echo $SLURM_ARRAY_TASK_ID
((SLURM_ARRAY_TASK_ID--))
echo $SLURM_ARRAY_TASK_ID
SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID + 0))
echo $SLURM_ARRAY_TASK_ID

export MACRO_DIR=/lustre/home/user/v/vtroshin/polarization-analysis-framework-main/train_macros
export OUT_DIR=/lustre/home/user/v/vtroshin/output/MCTests/7bins

if [ -f "${OUT_DIR}/Anal_bin${SLURM_ARRAY_TASK_ID}.root" ]; then
  exit 0
fi

mkdir d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
cd d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} 
cp ${MACRO_DIR}/MpdGlobalPolarization.cxx MpdGlobalPolarization.cxx
cp ${MACRO_DIR}/MpdGlobalPolarization.h MpdGlobalPolarization.h
cp ${MACRO_DIR}/MpdAnalysisTask2.cxx MpdAnalysisTask2.cxx
cp ${MACRO_DIR}/MpdAnalysisTask2.h MpdAnalysisTask2.h
cp ${MACRO_DIR}/mpdloadlibs.C mpdloadlibs.C
cp ${MACRO_DIR}/nTr_Centr_Req30-PHSD.root nTr_Centr_Req30-PHSD.root
cp ${MACRO_DIR}/pCentr.txt pCentr.txt
cp ${MACRO_DIR}/TrackRecEff.root TrackRecEff.root
cp ${MACRO_DIR}/pGlobalPol.txt pGlobalPol.txt
cp ${MACRO_DIR}/pEP.txt pEP.txt

DIRECTORY=/lustre/home/mpd/amoshkin/req30/data

printf -v SLURM_ARRAY_TASK_ID_COUNTER "%04d" $SLURM_ARRAY_TASK_ID

export INFILE0=${DIRECTORY}/phsd-BiBi-09.2GeV-mb-eos0-2k-${SLURM_ARRAY_TASK_ID_COUNTER}-0.reco.root
export INFILE1=${DIRECTORY}/phsd-BiBi-09.2GeV-mb-eos0-2k-${SLURM_ARRAY_TASK_ID_COUNTER}-1.reco.root
export INFILE2=${DIRECTORY}/phsd-BiBi-09.2GeV-mb-eos0-2k-${SLURM_ARRAY_TASK_ID_COUNTER}-2.reco.root
export INFILE3=${DIRECTORY}/phsd-BiBi-09.2GeV-mb-eos0-2k-${SLURM_ARRAY_TASK_ID_COUNTER}-3.reco.root

echo $INFILE0 > list.txt
echo $INFILE1 >> list.txt
echo $INFILE2 >> list.txt
echo $INFILE3 >> list.txt

echo "---------------"
echo "Start Analysis"
date
echo "---------------"

export OUTFILE=${OUT_DIR}/Anal_bin${SLURM_ARRAY_TASK_ID}

root -b -q MpdAnalysisTask2.cxx MpdGlobalPolarization.cxx -e "runit(\"$OUTFILE\");"

cd ..
rm -rf d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

echo "End"
echo "---------------"
date
