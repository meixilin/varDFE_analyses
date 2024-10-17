#!/bin/bash
#$ -l highp,h_data=4G,h_vmem=INFINITY,h_rt=72:00:00
#$ -pe shared 25

# @version      v2
# @usage        qsub -t <rowid> DFE1D_refspectra_final.sh
# @description  wrapper to run python DFE1D_refspectra.py by workflow reference table
# Author: Meixi Lin
# Date: 2023-03-23 18:22:28

################################################################################
## import packages
sleep $((RANDOM % 20))
set -eo pipefail


conda activate dfe

################################################################################
## def variables
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse HEAD)
WORKDIR="${HOMEDIR}/data_varDFE"

REFDICT="DFE1D_refspectra_final.txt"

# find values
pop=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $4}' ${REFDICT})
mask_singleton=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $5}' ${REFDICT})
demog_model=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $6}' ${REFDICT})
demog_params=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $7}' ${REFDICT})
ns=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $8}' ${REFDICT})
outprefix=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $9}' ${REFDICT})

outdir=$(echo ${outprefix} | cut -d '/' -f1-5)
logprefix=$(echo ${outprefix} | cut -d '/' -f6)
mkdir -p "${WORKDIR}/${outdir}"

TODAY=$(date "+%Y%m%d")
LOG="${WORKDIR}/${outdir}/${logprefix}_${TODAY}.log"

################################################################################
## def functions

################################################################################
## main

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}" > ${LOG}

cd ${WORKDIR}

python "${HOMEDIR}/scripts_varDFE/varDFE/workflow/DFE/DFE1D_refspectra.py" \
"$demog_model" "$demog_params" "$ns" "$outprefix" &>> "${LOG}"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}


