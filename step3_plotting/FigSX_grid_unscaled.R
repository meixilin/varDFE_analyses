# Title: Plot full gridsearch outputs following Fig_Gridsearch
# Author: Meixi Lin
# Date: Wed May  8 11:57:20 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(cowplot)
library(data.table)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

# def functions --------


# def variables --------
prefix = 'FigSX_grid_unscaled'
dir.create(prefix)

# load data --------
grid_df = readRDS(file = 'TabSX_gridsearch/unscaled/11pops_gamma_unscaled_FREEgrid.rds')
load('Fig_Gridsearch/Fig_Gridsearch.RData')

# main --------
allplots <- lapply(pops, plot_grids_pop)
leg <-  ggpubr::get_legend(allplots[[1]], position = 'right')

pp <- cowplot::plot_grid(plotlist = allplots, align = 'hv', axis = 'tblr', ncol = 3, labels = 'AUTO', label_size = labelsize)

# add legend while plotting
outpp <- pp +
    draw_plot(leg, x = 0.7, y = 0, width = 0.3, height = 0.3)

# output files --------
ggsaver('FigSX_grid_unscaled.pdf', prefix, outpp, 18, 22)

# save outputs
rm(grid_df)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()

