#!/bin/bash
#$ -l h_data=2G,h_vmem=INFINITY,h_rt=23:59:00
#$ -pe shared 20

# @version      v2
# @usage        qsub -t <rowid> DFE1D_gridsearch_final.sh
# @description  wrapper to run python DFE1D_gridsearch.py by workflow reference table
# Author: Meixi Lin
# Date: 2023-05-10 15:40:03

################################################################################
## import packages
sleep $((RANDOM % 60))
set -eo pipefail


conda activate dfe

################################################################################
## def variables
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse HEAD)
WORKDIR="${HOMEDIR}/data_varDFE"

REFDICT="DFE1D_gridsearch_final.txt"

# find values
pop=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $4}' ${REFDICT})
max_bound=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $5}' ${REFDICT})
min_bound=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $6}' ${REFDICT})
dfe_scaling=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $7}' ${REFDICT})
Nanc=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $8}' ${REFDICT})
mask_singleton=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $9}' ${REFDICT})
# demog_model=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $10}' ${REFDICT})
sfs=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $11}' ${REFDICT})
ref_spectra=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $12}' ${REFDICT})
pdfname=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $13}' ${REFDICT})
theta_nonsyn=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $14}' ${REFDICT})
outprefix=$(awk -v sgeid=${SGE_TASK_ID} '$1 == sgeid {print $15}' ${REFDICT})

outdir=$(echo ${outprefix} | cut -d '/' -f1-4)
logsuffix=$(echo ${outprefix} | cut -d '/' -f5)
mkdir -p "${WORKDIR}/${outdir}"

TODAY=$(date "+%Y%m%d")
SFS="${WORKDIR}/SFS/${sfs}"
LOG="${WORKDIR}/${outdir}/${logsuffix}_${TODAY}.log"

################################################################################
## def functions
SCALEDRUN="python ${HOMEDIR}/scripts_varDFE/varDFE/workflow/DFE/DFE1D_gridsearch.py --max_bound $max_bound --min_bound $min_bound --Npts 2000"
UNSCALEDRUN="$SCALEDRUN --dfe_scaling --Nanc $Nanc"

################################################################################
## main

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}" > ${LOG}

cd ${WORKDIR}

if [ $mask_singleton == 'TRUE' ] && [ $dfe_scaling == 'TRUE' ]
then
    echo $UNSCALEDRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix"
    $UNSCALEDRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix" &>> "${LOG}"
elif [ $mask_singleton == 'FALSE' ] && [ $dfe_scaling == 'TRUE' ]
then
    echo $UNSCALEDRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix"
    $UNSCALEDRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix" &>> "${LOG}"
elif [ $mask_singleton == 'TRUE' ] && [ $dfe_scaling == 'FALSE' ]
then
    echo $SCALEDRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix"
    $SCALEDRUN --mask_singleton "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix" &>> "${LOG}"
elif [ $mask_singleton == 'FALSE' ] && [ $dfe_scaling == 'FALSE' ]
then
    echo $SCALEDRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix"
    $SCALEDRUN "$SFS" "$ref_spectra" "$pdfname" "$theta_nonsyn" "$outprefix" &>> "${LOG}"
else
    echo -e "[$(date "+%Y-%m-%d %T")] Wrong dfe_scaling = ${dfe_scaling} or mask_singleton = ${mask_singleton}"
    exit 3
fi


echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}


