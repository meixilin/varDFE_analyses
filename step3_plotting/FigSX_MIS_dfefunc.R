# Title: Figure SX: MIS SFS fit with each given DFE functions
# Author: Meixi Lin
# Date: Wed Apr 27 23:14:02 2022
# Modification: Update to 11 species
# Date: Tue Mar 19 09:27:01 2024


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

theme_set(theme_bw(base_size = 7))

# def functions --------
load_fitsfs <- function(dfe_func) {
    fit_sfs = vector('list',length=npop)
    for (ii in seq_along(pops)) {
        pop <- pops[ii]
        mypath = paste0('../DFE_inferenceFIM/final/mask_',as.character(mask_singletons[ii]),
                        '/', pop,'/', demog_models[ii],'/',dfe_func,'/bestrun')
        modelsfsf = list.files(path = mypath, pattern = "_folded.expSFS", full.names = TRUE)
        message(paste0('Loading: ', modelsfsf))
        sfs <- load_dadi_sfs(dadisfs = modelsfsf, id = pop,
                             singleton_mask = mask_singletons[ii]) %>%
            dplyr::mutate(model=dfe_func,
                          mask=mask_singletons[ii])
        fit_sfs[[ii]] <- sfs
    }
    return(fit_sfs)
}

plot_sfs <- function(plotdt) {
    rdgy = rev(RColorBrewer::brewer.pal(n = 8, 'RdGy')[c(1,8)])
    pp <- ggplot(data = plotdt, aes(x = site_freq, y = count, fill = model)) +
        geom_col(position = 'dodge') +
        labs(x = 'Minor allele frequency', y = 'Count', fill = 'DFE function') +
        scale_fill_manual(values = rdgy) +
        facet_grid(sfsid ~ ., scales = 'free_y') +
        scale_x_continuous(breaks = c(seq(0,60,5),68)) +
        theme(legend.position = 'top',
              legend.key.size = unit(0.5,"line"))
    return(pp)
}

# def variables --------
prefix = 'FigSX_MIS_dfefunc'
dir.create(prefix)

dfe_funcs = c('gamma', 'neugamma', 'gammalet', 'lognormal', 'lourenco_eq')

# best models splits
bestmodels
mask_singletons = unname(sapply(bestmodels, function(xx) strsplit(xx, ' ')[[1]][2] == 'True'))
demog_models = unname(sapply(bestmodels, function(xx) strsplit(xx, ' ')[[1]][3]))

# load data --------
mis_sfsf = paste0('../SFS/MIS-',pops,'.sfs')

mis_sfs = vector('list',length=npop)

for (kk in seq_along(mask_singletons)) {
    pop <- pops[kk]
    sfs <- load_dadi_sfs(dadisfs=mis_sfsf[kk],id=pop,
                         singleton_mask = mask_singletons[kk]) %>%
        dplyr::mutate(model='data',
                      mask=mask_singletons[kk])
    mis_sfs[[kk]] <- sfs
}

# load best models and plot --------
for (dfe_func in dfe_funcs) {
    fit_sfs = load_fitsfs(dfe_func)
    plotdt = dplyr::bind_rows(mis_sfs, fit_sfs)
    plotdt$sfsid = factor(plotdt$sfsid,levels = pops,labels = common_names)
    plotdt$model = factor(plotdt$model,levels = c('data',dfe_funcs))
    pp = plot_sfs(plotdt)
    ggsaver(paste0('FigSX_MIS_', dfe_func, '.pdf'), prefix, pp, 16,22)
    save.image(file = paste0(prefix,'/FigSX_MIS_',dfe_func,'.RData'))
}

# cleanup --------
date()
closeAllConnections()


