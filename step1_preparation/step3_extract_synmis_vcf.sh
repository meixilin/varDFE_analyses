#!/bin/bash

################################################################################
### Select individuals to keep for the SYN-VCF and MIS-VCF
gatk3 -Xmx15G -T SelectVariants \
-R "${REFERENCE}" \
--sample_name "sample1" \
--sample_name "sample2" \
--sample_name "sample3" \
--sample_name "sampleN" \
--setFilteredGtToNocall \
--excludeFiltered \
--removeUnusedAlternates \
-L "${CDSREGION}" \
-V "${VCFN}_ann_gatk_filtered.vcf.gz" \
-o "${VCFN}_ann_gatk_PASS_samp.vcf.gz"

bcftools +counts "${VCFN}_ann_gatk_PASS_samp.vcf.gz"

################################################################################
### Filter for excessive heterozygosity or missingness

bcftools +fill-tags "${VCFN}_ann_gatk_PASS_samp.vcf.gz" -- --tags 'NS' | \
bcftools filter -i 'COUNT(GT=="het")/NS < 0.75' -s 'WARNING_excHet75' -m + | \
bcftools filter -i 'F_MISSING < 0.2' -s 'WARNING_missing20' -m + | \
bcftools view -f '.,PASS' -O z -o "${VCFN}_allsites_samp.vcf.gz"

bcftools stats -f 'PASS,.' "${VCFN}_allsites_samp.vcf.gz"

################################################################################
### Keep only biallelic SNPs to reduce file size

gatk3 -Xmx15G -T SelectVariants \
-R "${REFERENCE}" \
--selectTypeToInclude "SNP" \
--excludeFiltered \
--preserveAlleles \
-V "${VCFN}_allsites_samp.vcf.gz" \
-o "${VCFN}_SNPs_samp.vcf.gz"

################################################################################
### Generate SYN-VCF and MIS-VCF
query_ann() {
    bcftools query -f '%INFO/ANN\n' ${1} | cut -d '|' -f 2 | sort | uniq -c
}

# SYN-VCF
java -jar -Xmx15G SnpSift.jar filter \
"ANN[0].EFFECT == 'synonymous_variant'" \
"${VCFN}_SNPs_samp.vcf.gz" > \
"${VCFN}_SYN.vcf"

bgzip "${VCFN}_SYN.vcf"
tabix -p vcf "${VCFN}_SYN.vcf.gz"

# check that all selected are synonymous_variant
query_ann "${VCFN}_SYN.vcf.gz"

# MIS-VCF
java -jar -Xmx15G SnpSift.jar filter \
"ANN[0].EFFECT == 'missense_variant'" \
"${VCFN}_SNPs_samp.vcf.gz" > \
"${VCFN}_MIS.vcf"

bgzip "${VCFN}_MIS.vcf"
tabix -p vcf "${VCFN}_MIS.vcf.gz"

# check that all selected are missense_variant
query_ann "${VCFN}_MIS.vcf.gz"
