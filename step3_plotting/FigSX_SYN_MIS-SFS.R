# Title: Figure SX: SYN and MIS SFS obtained for each species
# Author: Meixi Lin
# Date: Wed Apr 27 23:14:02 2022
# Modification: Update to 11 species
# Date: Mon Mar 18 14:19:11 2024

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
indir = '../SFS/'
prefix = 'FigSX_SYN_MIS-SFS'
dir.create(prefix)

# load data --------
syn_sfs = vector('list',length=npop)

for (ii in seq_along(pops)) {
    pop <- pops[ii]
    sfs <- load_dadi_sfs(dadisfs=paste0(indir,'SYN-',pop,'.sfs'),id=pop) %>%
        dplyr::mutate(type='SYN')
    syn_sfs[[ii]] <- sfs
}

mis_sfs = vector('list',length=npop)

for (ii in seq_along(pops)) {
    pop <- pops[ii]
    sfs <- load_dadi_sfs(dadisfs=paste0(indir,'MIS-',pop,'.sfs'),id=pop) %>%
        dplyr::mutate(type='MIS')
    mis_sfs[[ii]] <- sfs
}

plotdt = dplyr::bind_rows(syn_sfs, mis_sfs)
plotdt$sfsid = factor(plotdt$sfsid,levels=pops, labels = common_names)

# main --------
pp <- ggplot(data = plotdt, aes(x = site_freq, y = count, fill = type)) +
    geom_col(position = 'dodge') +
    labs(x = 'Minor allele frequency', y = 'Count', fill = 'Mutation type') +
    facet_grid(sfsid ~ ., scales = 'free_y') +
    scale_x_continuous(breaks = c(seq(0,60,5),68)) +
    theme(legend.position = 'top',
          legend.key.size = unit(0.5,"line"))

# output files --------
ggsaver(paste0(prefix, '.pdf'), prefix, pp, 16,22)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()


