#!/bin/bash

################################################################################
### Download VCF files

################################################################################
### Add gene annotations if needed
# Build snpEff database
java -Xmx10G -jar snpEff.jar build -debug "${DATABASE}"

# Run snpEff on VCF files
java -Xmx15G -jar snpEff.jar -nodownload -v \
-stats "${VCFN}_ann.html" \
"${DATABASE}" \
"${VCFN}.vcf.gz" > \
"${VCFN}_ann.vcf"

bgzip "${VCFN}_ann.vcf"
tabix -p vcf "${VCFN}_ann.vcf.gz"

################################################################################
### Annotate VCF files with standard annotations
gatk3 -Xmx20G -Djava.io.tmpdir="./temp" \
-T VariantAnnotator \
-R "${REFERENCE}" \
-G StandardAnnotation \
-A VariantType \
-A AlleleBalance \
-A ClippingRankSumTest \
-A GenotypeSummaries \
-L "${AUTOSOMES}" \
-V "${VCFN}_ann.vcf.gz" \
-o "${VCFN}_ann_gatk.vcf.gz"

################################################################################
### Mask the repeat and CpG regions
gatk3 -Xmx20G -Djava.io.tmpdir="./temp" \
-T VariantFiltration \
-R "${REFERENCE}" \
--logging_level ERROR \
--mask "${RMCPGMASK}" --maskName "FAIL_CpGRep" \
-L "${AUTOSOMES}" \
-V "${VCFN}_ann_gatk.vcf.gz" \
-o "${VCFN}_ann_gatk_filtered.vcf.gz"

################################################################################
### Apply site-level and genotype-level filters if needed
# gatk standard filters
gatk3 -Xmx20G -Djava.io.tmpdir="./temp" \
-T VariantFiltration \
-R "${REFERENCE}" \
--logging_level ERROR \
--mask "${RMCPGMASK}"--maskName "FAIL_CpGRep" \
-filter "QD < 2.0" --filterName "FAIL_QD" \
-filter "FS > 60.0" --filterName "FAIL_FS" \
-filter "MQ < 40.0" --filterName "FAIL_MQ" \
-filter "MQRankSum < -12.5" --filterName "FAIL_MQRS" \
-filter "ReadPosRankSum < -8.0" --filterName "FAIL_RPRS" \
-filter "SOR > 3.0" --filterName "FAIL_SOR" \
-filter "DP < [minDP]" --filterName "FAIL_minDP" \
-filter "DP > [maxDP]" --filterName "FAIL_maxDP" \
-filter "VariantType != 'NO_VARIATION' && VariantType != 'SNP'" --filterName "FAIL_mutType" \
-L "${AUTOSOMES}" \
-V "${VCFN}_ann_gatk.vcf.gz" \
-o "${VCFN}_ann_gatk_f1.vcf.gz"

# custom genotype-level filters
python customVCFfilter.py \
"${VCFN}_ann_gatk_f1.vcf.gz" \
"${VCFN}_ann_gatk_filtered.vcf" \
"geno_maxD.csv"

bgzip "${VCFN}_ann_gatk_filtered.vcf"
tabix -p vcf "${VCFN}_ann_gatk_filtered.vcf.gz"

# use bcftools to get coverage and percent missing statistics
bcftools stats -f "PASS,." -s- "${VCFN}_ann_gatk_filtered.vcf.gz"

