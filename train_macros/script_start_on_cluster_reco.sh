#!/bin/bash

#
#SBATCH -D /lustre/home/user/v/vtroshin/output/Reco/Analysis/Chi
#SBATCH --mem-per-cpu=8000M
#SBATCH -J Reco
#SBATCH -p mephi
#SBATCH -t 1-4
#SBATCH -a 1-2499
#SBATCH -x n02p002,n02p003,n02p004,n02p005
#
#SBATCH -o /lustre/home/user/v/vtroshin/output/Reco/Analysis/Chi/logs/slurm_%A_%a.out
#SBATCH -e /lustre/home/user/v/vtroshin/output/Reco/Analysis/Chi/logs/slurm_%A_%a.err
#

date

source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev
export MPDROOT=/lustre/home/user/n/nazarova/mpd
source $MPDROOT/config/env.sh

module add GVR/v1.0-1

echo $SLURM_ARRAY_TASK_ID
((SLURM_ARRAY_TASK_ID--))
echo $SLURM_ARRAY_TASK_ID
SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID + 0))
#7476
echo $SLURM_ARRAY_TASK_ID

export MACRO_DIR=/lustre/home/user/v/vtroshin/polarization-analysis-framework-main/train_macros
export OUT_DIR=/lustre/home/user/v/vtroshin/output/Reco/Analysis/Chi

if [ -f "${OUT_DIR}/Anal_bin${SLURM_ARRAY_TASK_ID}.root" ]; then
  exit 0
fi

mkdir d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
cd d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} 
cp ${MACRO_DIR}/MpdGlobalPolarization_RECO.cxx MpdGlobalPolarization_RECO.cxx
cp ${MACRO_DIR}/MpdGlobalPolarization_RECO.h MpdGlobalPolarization_RECO.h
cp ${MACRO_DIR}/MpdAnalysisTask2.cxx MpdAnalysisTask2.cxx
cp ${MACRO_DIR}/MpdAnalysisTask2.h MpdAnalysisTask2.h
cp ${MACRO_DIR}/mpdloadlibs.C mpdloadlibs.C
cp ${MACRO_DIR}/RunAnalyses.C RunAnalyses.C
cp ${MACRO_DIR}/nTr_Centr_Req30-PHSD.root nTr_Centr_Req30-PHSD.root
cp ${MACRO_DIR}/pCentr.txt pCentr.txt
cp ${MACRO_DIR}/pEP.txt pEP.txt
cp ${MACRO_DIR}/pGlobalPol_RECO.txt pGlobalPol_RECO.txt
cp ${MACRO_DIR}/TrackRecEff.root TrackRecEff.root
cp ${MACRO_DIR}/Struct_L0_v1.h Struct_L0_v1.h

cp ${MACRO_DIR}/ChiSelection_values_Cent_4.txt ChiSelection_values_Cent_4.txt
cp ${MACRO_DIR}/Omega2_values_Cent_4.txt Omega2_values_Cent_4.txt

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

rootcling -v4 -f mydict.cxx -rmf libmydict.rootmap -rml libmydict.so Struct_L0_v1.h
g++ -shared -o libmydict.so mydict.cxx `root-config --cflags --libs` -fPIC

root -b -q MpdAnalysisTask2.cxx MpdGlobalPolarization_RECO.cxx RunAnalyses.C\(\"${OUTFILE}\",\"analysis\",\"omega2\"\)

cd ..
rm -rf d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

echo "End"
echo "---------------"
date
