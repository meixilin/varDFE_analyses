# Title: utilities for dfe plots
# Author: Meixi Lin
# Date: Fri Mar 31 10:35:32 2023
# Last modified: Tue Sep 10 18:43:28 2024

require(dplyr)
require(RColorBrewer)
require(ggplot2)
require(mmutil) # loadSomeRData, SourceSome and write_table
options(stringsAsFactors = FALSE, warn = 1)

# working directory --------
workdir = "workdir"

# variables --------
npop = 11
pops = c("AC136","DM100","LC18","FH18","FA100","CL24","CA26","PS24","BP44","MM16","HS100")
popinsects = c('AC136', 'DM100', 'LC18')
popbirds = c('FH18', 'FA100')
popmammals = c("CL24","CA26","PS24","BP44","MM16","HS100")
# common name labels
common_names = c("Mosquito","Drosophila","Halictid bee","Pied flycatcher","Collared flycatcher",
                 "Gray wolf","Arctic wolf","Vaquita","Fin whale","Mouse","Human")
# pals::pal.bands(popcols)
# cutoff for dfe integration check
pintverr = 1e-3

# best model selected for main text presentation
bestmodels = c("AC136 True three_epoch","DM100 False two_epoch","LC18 True two_epoch",
               "FH18 False two_epoch","FA100 False three_epoch","CL24 False three_epoch",
               "CA26 False two_epoch","PS24 False two_epoch","BP44 False two_epoch",
               "MM16 False two_epoch","HS100 False two_epoch")

## for robustness comparisons
bestmodelsdf = reshape2::colsplit(bestmodels, ' ', c('pop', 'mask_singleton', 'demog_model'))

# categories for s or 2ns ========
# unscaled
cats_us = 10^(seq(-5,-2))
# technically the last one is [1e-2,Inf] but > 0.5 is not possible so clumped
catlabs_us = c('[0, 1e-5)', '[1e-5, 1e-4)', '[1e-4, 1e-3)', '[1e-3, 1e-2)', '[1e-2, 0.5]')

# scaled
cats_s = 10^(seq(0,3))
catlabs_s = c('[0, 1)', '[1, 10)', '[10, 100)', '[100, 1000)', '[1000, Inf]')

# reading data --------
load_dadi_sfs <- function(dadisfs, id, singleton_mask = FALSE) {
    sfs <- read.table(dadisfs,skip = 1,header = FALSE, colClasses = 'numeric')
    rownames(sfs) = c('count', 'mask')
    nsites = dim(sfs)[2]-1
    sfs = t(sfs) %>%
        as.data.frame() %>%
        dplyr::mutate(site_freq = seq(0,nsites),
                      sfsid = id)
    if (singleton_mask == TRUE) {
        sfs[sfs$site_freq == 1, 'mask'] = 1
    }
    sfs = sfs %>%
        dplyr::filter(mask == 0)

    sfssum = sum(sfs$count)
    sfs = sfs %>%
        dplyr::mutate(prop = count/sfssum)
    return(sfs)
}

# filtercols: columns to filter on
# filterstrings: the strings to filter when filtercols are pasted together
# .filterid: should the pasted filtercols be kept in the data frame
load_DFEdf <- function(dfe_func, dftype = c('table', 'prop_us', 'prop'),
                       filtercols = "modelstring", filterstrings = NULL, .filterid = FALSE) {
    # check file inputs
    myfile = paste0('TabSX_DFE/', dfe_func, '/11pops_', dfe_func, '_', dftype, '.rds')
    stopifnot(file.exists(myfile))
    df = readRDS(file = myfile)
    # filter by the standards given
    if (!is.null(filterstrings)) {
        if (length(filtercols) > 1) {
            stopifnot(all(filtercols %in% colnames(df)))
            # generate new "modelstring"
            df$filterid = apply(df[,filtercols], 1, paste, collapse = ' ')
        } else {
            if (filtercols == "modelstring") {
                stopifnot(all(c('pop', 'mask_singleton', 'demog_model') %in% colnames(df)))
                # define filterid regardless if modelstring column exists
                df = df %>% dplyr::mutate(filterid = paste(pop, mask_singleton, demog_model))
            } else {
                stopifnot(all(filtercols %in% colnames(df)))
                df$filterid = df[,filtercols]
            }
        }
    }
    # print(df$filterid)
    df = df %>%
        dplyr::filter(filterid %in% filterstrings)
    # should be selecting for more than the length of modelstrings
    stopifnot(nrow(df) >= length(filterstrings))
    # remove filterid column by default
    if(!.filterid) df = df %>% dplyr::select(-filterid)
    return(df)
}

# pdf and stats ---------
# pgamma derived methods of obtaining pneugamma
# pneugamma ========
# pneugamma based on pgamma (here mins is not used)
pneugamma <- function(pneu, shape, scale, mins,
                      cats = cats_us,
                      catlabs = catlabs_us) {
    if(max(length(pneu), length(scale), length(shape)) > 1)
        stop("parameters must be of length 1")
    if (any(is.na(c(pneu, shape, scale)))) {
        pintv = rep(NA, times = length(catlabs))
    } else {
        # cdf for each sections (slightly different from pneugamma implementation)
        pcats = sapply(c(cats, Inf) , function(xx) {
            stats::pgamma(q = xx, shape = shape, scale = scale)
        })
        # get intervals
        pintv = c(pcats[1],pcats[-1]-pcats[-length(pcats)])
        # scale for pneu
        pintv = (1-pneu) * pintv
        pintv[1] = pintv[1] + pneu
        # sanity check
        if (abs(sum(pintv)-1) > pintverr) {
            stop('Given pdf does not integrate to 1')
        }
    }
    # format data
    outdt = data.frame(catlabs,pintv)
    colnames(outdt) = c('cat','prop')
    return(outdt)
}

pneugamma_s = pneugamma
formals(pneugamma_s)$cats = cats_s
formals(pneugamma_s)$catlabs = catlabs_s

# gamma ========
# pgamma with categories
pgamma_dfe <- function(shape, scale,
                       cats = cats_us,
                       catlabs = catlabs_us) {
    if(max(length(scale), length(shape)) > 1)
        stop("parameters must be of length 1")
    if (any(is.na(c(shape, scale)))) {
        pintv = rep(NA, times = length(catlabs))
    } else {
        # cdf for each sections (slightly different from pneugamma implementation)
        pcats = sapply(c(cats, Inf) , function(xx) {
            stats::pgamma(q = xx, shape = shape, scale = scale)
        })
        # get intervals
        pintv = c(pcats[1],pcats[-1]-pcats[-length(pcats)])
        # sanity check
        if (abs(sum(pintv)-1) > pintverr) {
            stop('Given pdf does not integrate to 1')
        }
    }
    # format data
    outdt = data.frame(catlabs,pintv)
    colnames(outdt) = c('cat','prop')
    return(outdt)
}

pgamma_dfe_s = pgamma_dfe
formals(pgamma_dfe_s)$cats = cats_s
formals(pgamma_dfe_s)$catlabs = catlabs_s

# gammalet ========
# define gammalet based on gamma distribution
pgammalet <- function(plet, shape, scale,
                      cats = cats_us,
                      catlabs = catlabs_us) {
    if(max(length(plet), length(scale), length(shape)) > 1)
        stop("parameters must be of length 1")
    if (any(is.na(c(plet, shape, scale)))) {
        pintv = rep(NA, times = length(catlabs))
    } else {
        # cdf for each sections, scale by the percent lethal so it sums up to 1-plet
        pcats = sapply(c(cats, Inf) , function(xx) {
            stats::pgamma(q = xx, shape = shape, scale = scale)*(1-plet)
        })
        # get intervals
        pintv = c(pcats[1],pcats[-1]-pcats[-length(pcats)])
        # for gammalet only, adds in the last plet
        pintv[length(pintv)] = pintv[length(pintv)] + plet
        # sanity check
        if (abs(sum(pintv)-1) > pintverr) {
            stop('Given pdf does not integrate to 1')
        }
    }
    # format data
    outdt = data.frame(catlabs,pintv)
    colnames(outdt) = c('cat','prop')
    return(outdt)
}

# scaled pgammalet
pgammalet_s = pgammalet
formals(pgammalet_s)$cats = cats_s
formals(pgammalet_s)$catlabs = catlabs_s

# lognormal ========
plognormal <- function(meanlog, sdlog,
                       cats = cats_us,
                       catlabs = catlabs_us) {
    if(max(length(meanlog), length(sdlog)) > 1)
        stop("parameters must be of length 1")
    if (any(is.na(c(meanlog, sdlog)))) {
        pintv = rep(NA, times = length(catlabs))
    } else {
        # cdf for each sections
        pcats = sapply(c(cats, Inf) , function(xx) {
            stats::plnorm(q= xx, meanlog = meanlog,sdlog = sdlog)
        })
        # get intervals
        pintv = c(pcats[1],pcats[-1]-pcats[-length(pcats)])
        # sanity check
        if (abs(sum(pintv)-1) > pintverr) {
            stop('Given pdf does not integrate to 1')
        }
    }
    # format data
    outdt = data.frame(catlabs,pintv)
    colnames(outdt) = c('cat','prop')
    return(outdt)
}

plognormal_s = plognormal
formals(plognormal_s)$cats = cats_s
formals(plognormal_s)$catlabs = catlabs_s

# utility --------
get_min0 <- function(x,sd,xname) {
    # the allowed parameters for get_min
    xnamelist = c("shape", "scale", "pneu", "plet", "mus", "sigma", "m", "Ne")
    # mus can be negative by definitions
    # shape, scale, sigma (in both lognormal and lourenco_eq), m, Ne > 0
    xname_nozero = c("shape", "scale", "sigma", "m", "Ne")
    # pneu, plet >= 0
    xname_zero = c("pneu", "plet")
    stopifnot(xname %in% xnamelist)
    if (is.na(sd)) {
        xmin = NA
    } else {
        xmin = x-sd
        # convert less than 0 values to NA if the given parameter has no definition at zero
        if (xmin < 0) {
            if (xname %in% xname_nozero) xmin = NA
            if (xname %in% xname_zero) xmin = 0
            # keep mus the same as it can be negative
            # output warnings
            warning(paste0(xname, ' x - sd is less than zero. x = ', x, ' sd = ', sd, ' xmin = ', xmin))
        }
    }
    return(xmin)
}

get_min <- Vectorize(get_min0)

get_max0 <- function(x,sd) {
    if (is.na(sd)) {
        xmax = NA
    } else {
        xmax = x+sd
    }
    return(xmax)
}

get_max <- Vectorize(get_max0)

# append taxonomic groups for populations
append_group <- function(df) {
    outdf = df %>%
        dplyr::mutate(group = case_when(pop %in% popinsects ~ 'Insecta',
                                        pop %in% popbirds ~ 'Aves',
                                        TRUE ~ 'Mammalia')) %>%
        dplyr::mutate(group = factor(group, levels = c('Insecta', 'Aves', 'Mammalia')))
    return(outdf)
}

# plotting utilities and standardization --------
# new color scheme
popcols = ggthemes::tableau_color_pal(palette = "Hue Circle")(19)[c(19,1,3,6,7,9,10,12,13,15,17)]
names(popcols) = pops

# plot theme
theme_set(cowplot::theme_cowplot(font_size = 8))
# label size
labelsize = 9 # A,B,C,D's sizes
# geom_point size default
pointsize = 3
# line width default
linewidth = 0.5
# line width for columns
linewidthcol = 0.1
# alpha default
alphavalue = 0.2
# default text size within plots
textsize = 2.5
# error bar color
linerangecol = 'darkgray'

# wrapper for ggsave so all plots' output width is the same
ggsaver <- function(myname, mypath, myplot, mywidth, myheight) {
    ggsave(filename = myname, path = mypath, plot = myplot, width = mywidth, height = myheight, units = 'cm', device = cairo_pdf, bg = "transparent")
}


# SCRATCH --------
# old color scheme
# # use ggplot default colors with slight twist (add more divergence btw insect - bird, bird - mammal)
# popcols = scales::hue_pal(h = c(0,300), l = 50)(13)[-c(4,7)]
# set up a legend for every species






