#!/bin/bash
#$ -l highp,h_data=8G,h_vmem=10G,h_rt=72:00:00

# @version      v2
# @usage        qsub -t <rowid> Demog1D_sizechangeFIM_final.sh
# @description  wrapper to run python Demog1D_sizechangeFIM.py by workflow reference table
# Author: Meixi Lin
# Date: 2023-03-21 11:20:50

################################################################################
## import packages
sleep $((RANDOM % 60))
set -eo pipefail

conda activate dfe

################################################################################
## def variables
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse HEAD)
WORKDIR="${HOMEDIR}/data_varDFE"

REFDICT="Demog1D_sizechangeFIM_final.txt"

# find values
pop=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $4}' ${REFDICT})
mu=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $5}' ${REFDICT})
Lcds=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $6}' ${REFDICT})
NS_S_scaling=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $7}' ${REFDICT})
initval=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $8}' ${REFDICT})
mask_singleton=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $9}' ${REFDICT})
sfs=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $10}' ${REFDICT})
demog_model=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $11}' ${REFDICT})
outdir=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $12}' ${REFDICT})

# assemble variables
TODAY=$(date "+%Y%m%d")
SFS="${WORKDIR}/SFS/${sfs}"
OUTDIR="${WORKDIR}/${outdir}"
LOG="${WORKDIR}/${outdir}/${pop}_demog_${demog_model}_${TODAY}.log"

# remove directory if it exists
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

mkdir -p "$OUTDIR"

################################################################################
## def functions
# requires all the model to finish within 4 hrs
DEFAULTRUN="python ${HOMEDIR}/scripts_varDFE/varDFE/workflow/Demography/Demog1D_sizechangeFIM.py --pop $pop --mu $mu --Lcds $Lcds --NS_S_scaling $NS_S_scaling --initval $initval --impatient 14400"

################################################################################
## main

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}" > ${LOG}

cd ${WORKDIR}

# add mask_singleton flag if required to mask the singletons
if [ $mask_singleton == 'TRUE' ]; then
    echo $DEFAULTRUN --mask_singleton "$SFS" "$demog_model" "$OUTDIR"
    $DEFAULTRUN --mask_singleton "$SFS" "$demog_model" "$OUTDIR" &>> "${LOG}"
elif [ $mask_singleton == 'FALSE' ]; then
    echo $DEFAULTRUN "$SFS" "$demog_model" "$OUTDIR"
    $DEFAULTRUN "$SFS" "$demog_model" "$OUTDIR" &>> "${LOG}"
else
    echo -e "mask_singleton=${mask_singleton} not found"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}


