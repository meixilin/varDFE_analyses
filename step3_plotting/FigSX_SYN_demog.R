# Title: Figure SX: SYN SFS demographic fit
# Author: Meixi Lin
# Date: Wed Apr 27 23:14:02 2022
# Modification: Update to 11 species
# Date: Mon Mar 18 15:18:06 2024


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

theme_set(theme_bw(base_size = 7))

# def functions --------

# def variables --------
prefix = 'FigSX_SYN_demog'
dir.create(prefix)

demog_models = c('one_epoch','two_epoch','three_epoch', 'four_epoch')
mask_singletons = c(FALSE,TRUE)

# load data --------
syn_sfs = vector('list',length=npop*length(mask_singletons))
synsfsid=1
for (ii in seq_along(pops)) {
    for (kk in seq_along(mask_singletons)) {
        pop <- pops[ii]
        sfs <- load_dadi_sfs(dadisfs=paste0('../SFS/SYN-',pop,'.sfs'),id=pop,
                             singleton_mask = mask_singletons[kk]) %>%
            dplyr::mutate(model='data',
                          mask=mask_singletons[kk])
        syn_sfs[[synsfsid]] <- sfs
        synsfsid = synsfsid + 1
    }
}

# load best models ========
fit_sfs = vector('list',length=npop*length(demog_models)*length(mask_singletons))
fitsfsid = 1
for (ii in seq_along(pops)) {
    pop <- pops[ii]
    for (jj in seq_along(demog_models)) {
        demog_model <- demog_models[jj]
        for (kk in seq_along(mask_singletons)) {
            path = paste0('../Demography/final/mask_',as.character(mask_singletons[kk]),'/',
                          pop,'/',demog_model, '/bestrun/')
            modelsfsf = list.files(
                path = path,
                pattern = "_folded.expSFS")
            sfs <- load_dadi_sfs(dadisfs=paste0(path,modelsfsf),id=pop,
                                 singleton_mask = mask_singletons[kk]) %>%
                dplyr::mutate(model=demog_model,
                              mask=mask_singletons[kk])
            fit_sfs[[fitsfsid]] <- sfs
            fitsfsid = fitsfsid+1
        }
    }
}

# plotdt ========
plotdt = dplyr::bind_rows(syn_sfs, fit_sfs)
plotdt$sfsid = factor(plotdt$sfsid,levels=pops,labels = common_names)
plotdt$model = factor(plotdt$model,levels=c('data',demog_models))
plotdt$mask = factor(plotdt$mask,labels=c('singleton not masked', 'singleton masked'))

# main --------
rdgy = rev(brewer.pal(n = 8, 'RdGy')[c(1:4,8)])
pp <- ggplot(data = plotdt, aes(x = site_freq, y = prop, fill = model)) +
    geom_col(position = 'dodge') +
    labs(x = 'Minor allele frequency', y = 'Proportion', fill = 'Demographic model') +
    scale_fill_manual(values = rdgy) +
    facet_grid(sfsid ~ mask, scales = 'free_y') +
    scale_x_continuous(breaks = c(seq(0,60,5),68)) +
    theme(legend.position = 'top',
          legend.key.size = unit(0.5,"line"))

# output files --------
ggsaver(paste0(prefix, '.pdf'), prefix, pp, 22, 22)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()
