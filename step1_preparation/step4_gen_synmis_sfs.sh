#!/bin/bash

################################################################################
### Run projection previews using a modified easySFS module
python easySFS.py \
-i "${VCFN}_SYN.vcf.gz" \
-maxHetFilter 0.75 \
-p ${POPFILE} \
--preview -a -v \
> "${VCFN}_SYN_prj_preview.txt"

python easySFS.py \
-i "${VCFN}_MIS.vcf.gz" \
-maxHetFilter 0.75 \
-p ${POPFILE} \
--preview -a -v \
> "${VCFN}_MIS_prj_preview.txt"

################################################################################
### Generate final SYN-SFS and MIS-SFS given the selected NSITES
python direct1DSFS.py \
--proj ${NSITES} \
--pop ${VCFN} \
"${VCFN}_SYN.vcf.gz" \
${POPFILE} \
"${VCFN}_SYN.sfs"

python direct1DSFS.py \
--proj ${NSITES} \
--pop ${VCFN} \
"${VCFN}_MIS.vcf.gz" \
${POPFILE} \
"${VCFN}_MIS.sfs"

