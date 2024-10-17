#!/bin/bash
#$ -l highp,h_data=2G,h_vmem=INFINITY,h_rt=72:00:00
#$ -pe shared 20

# @version      v2
# @usage        qsub -t <rowid> DFE1D_inferenceFIM_final.sh
# @description  wrapper to run python DFE1D_inferenceFIM.py by workflow reference table
# Author: Meixi Lin
# Date: 2023-03-28 14:41:24

################################################################################
## import packages
sleep $((RANDOM % 20))
set -eo pipefail


conda activate dfe

################################################################################
## def variables
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse HEAD)
WORKDIR="${HOMEDIR}/data_varDFE"

REFDICT="DFE1D_inferenceFIM_final.txt"


# find values
pop=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $4}' ${REFDICT})
mu=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $5}' ${REFDICT})
Lcds=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $6}' ${REFDICT})
NS_S_scaling=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $7}' ${REFDICT})
mask_singleton=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $8}' ${REFDICT})
demog_model=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $9}' ${REFDICT})
sfs=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $10}' ${REFDICT})
ref_spectra=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $11}' ${REFDICT})
pdfname=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $12}' ${REFDICT})
theta_syn=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $13}' ${REFDICT})
outdir=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $14}' ${REFDICT})

mkdir -p "${WORKDIR}/${outdir}"

TODAY=$(date "+%Y%m%d")
SFS="${WORKDIR}/SFS/${sfs}"
LOG="${WORKDIR}/${outdir}/${pop}_${demog_model}_sel_mask_${mask_singleton}_${pdfname}_${TODAY}.log"

################################################################################
## def functions
DEFAULTRUN="python ${HOMEDIR}/scripts_varDFE/varDFE/workflow/DFE/DFE1D_inferenceFIM.py --pop $pop --mu $mu --Lcds $Lcds --NS_S_scaling $NS_S_scaling"

################################################################################
## main

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}" > ${LOG}

cd ${WORKDIR}

if [ $mask_singleton == 'TRUE' ]; then
    echo $DEFAULTRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_syn" "$outdir"
    $DEFAULTRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_syn" "$outdir" &>> "${LOG}"
elif [ $mask_singleton == 'FALSE' ]; then
    echo $DEFAULTRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_syn" "$outdir"
    $DEFAULTRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_syn" "$outdir" &>> "${LOG}"
else
    echo -e "mask_singleton=${mask_singleton} not found"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}


