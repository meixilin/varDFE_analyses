# Title: Population structure QC
# 1. Generate gdsfile from all autosomes
# 2. Select only SNPs that passed filters
# 3. LD pruning
# 4. Run PCA
# 5. Run kinship
# 6. Select individuals passing filters

# Author: Meixi Lin
# Date: 2023-02-20 14:57:33

options(echo = TRUE, stringsAsFactors = FALSE)
library(gdsfmt)
library(SNPRelate)

vcffn <- paste0(VCFN, "_ann_gatk_filtered")
gdsname <- paste0(vcffn, ".gds")
gdsnamef <- paste0(vcffn, "_SNPs_PASS.gds")
gdsnamefl <- paste0(vcffn, "_SNPs_PASS_pruned.gds")

################################################################################
### Convert VCF to GDS file
# use all the data from the annotated VCF
snpgdsVCF2GDS(vcf.fn = paste0(vcffn, ".vcf.gz"),
              out.fn = gdsname,
              method = "biallelic.only")

# write summary about this file
snpgdsSummary(paste0(vcffn, ".gds"))

# open this file
gf <- snpgdsOpen(paste0(vcffn, ".gds"), readonly = TRUE)

# use only pass sites
g <- (read.gdsn(index.gdsn(gf, "snp.annot/filter")) == "PASS")
snpsetf <- read.gdsn(index.gdsn(gf, "snp.id"))[g]
snpgdsCreateGenoSet(gdsname, gdsnamef, snp.id = snpsetf)
snpgdsClose(gf)
gff = snpgdsOpen(gdsnamef, readonly = TRUE)

# run ld pruning (use seed)
set.seed(7)
snpsetfl <- unname(unlist(snpgdsLDpruning(gff, autosome.only = FALSE)))
SNPRelate::snpgdsCreateGenoSet(gdsnamef, gdsnamefl, snp.id = snpsetfl)
snpgdsClose(gff)

################################################################################
### Run PCA analysis
# open the paste0(vcffn, "_SNPs_PASS_pruned.gds") file
genofile <- snpgdsOpen(gdsnamefl, readonly = TRUE)
pca <- snpgdsPCA(genofile, autosome.only = FALSE, verbose = TRUE,
                 sample.id = sampleset)
# examine subpopulations by plotting

################################################################################
### Run kinship analysis
ibddt <- snpgdsIBDKING(gdsobj = genofile,
                       sample.id = sampleset,
                       verbose = TRUE,
                       autosome.only = FALSE)


################################################################################
### Mark individuals that are in high kinship or outliers in the PCA
# returns character(0) if no outliers
pca_outliers <- function(pca) {
    indid <- apply(pca$eigenvect[, 1:3], 2, function(x) {
        which(abs(x - median(x)) > (6 * sd(x)))
    })
    indid <- unlist(indid)
    outliers <- pca$sample.id[indid]
    return(outliers)
}

sum_kinship <- function(kin) {
    samples = unique(c(kin$ID1,kin$ID2))
    avekin = unlist(lapply(samples, function(xx){
        nxx = which(kin$ID1 == xx | kin$ID2 == xx)
        if (length(nxx) != (length(samples)-1)) {
            stop('Wrong kinship shape.')
        }
        out = sum(kin[nxx, 'kinship'])
    }))
    names(avekin) = samples
    avekin = data.frame(avekin, stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column(var = 'sample.id') %>%
        dplyr::mutate(sample.id = as.character(sample.id))
    return(avekin)
}

get_highkin <- function(kin, cutoff) {
    removedt = kin %>%
        dplyr::filter(kinship > cutoff) %>%
        dplyr::arrange(ID1,ID2) %>%
        dplyr::mutate(Pair = paste(ID1,ID2,sep='-')) %>%
        dplyr::arrange(Pair)
    avekin = sum_kinship(kin = kin)
    removedt = dplyr::left_join(removedt, avekin, by = c('ID1'='sample.id'))
    removedt = dplyr::left_join(removedt, avekin, by = c('ID2'='sample.id'), suffix = c('.ID1','.ID2'))
    return(removedt)
}

pca_outliers(pca)
get_highkin(kin = kinship_all, cutoff = 0.15)
