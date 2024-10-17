# Title: Plot full gridsearch outputs for scaled DFE following Fig_Gridsearch and TabSX_Gridsearch
# Author: Meixi Lin
# Date: Thu May  9 15:19:52 2024

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
prefix = 'FigSX_grid_scaled'
dir.create(prefix)

myscale = 'scaled'

# load data --------
grid_df = readRDS(file = 'TabSX_gridsearch/scaled/11pops_gamma_scaled_FREEgrid.rds')
loadSomeRData(c('scale_vars','check_grids','add_grids','pairlrtcast', 'mledf', 'lrt_11pops'),
              'TabSX_gridsearch/scaled/TabSX_gridsearch_scaled_nogrid.RData')
loadSomeRData(c('generate_grid0', 'aggregate_grid0', 'plot_grids', 'plot_grids_pop', 'append_group2', 'common_names'),
              'Fig_Gridsearch/Fig_Gridsearch.RData')

# main --------
# plot grids for constrained model ========
grid0agg_11pops = aggregate_grid0(generate_grid0(grid_df, pops, lrt_11pops))

pp_Cons <- plot_grids(mygrid = grid0agg_11pops,
                      mledf = mledf,
                      ll_grid0 = lrt_11pops, grid0col = "black",
                      mytitle = 'Constrained')

# plot grids for individual populations ========
allplots <- lapply(pops, plot_grids_pop)
allleg <-  ggpubr::get_legend(allplots[[1]], position = 'top')

# plot pretty LRT clustering heatmap ========
pairlrtcast1 = pairlrtcast
rownames(pairlrtcast1) = unname(common_names)
colnames(pairlrtcast1) = unname(common_names)
raw_pairlrt <- heatmaply::heatmaply(pairlrtcast1,
                                    symm = TRUE,
                                    colors = heat.colors(n = 256),
                                    distfun = function(c) as.dist(c),
                                    hclustfun = function(d) hclust(d, method = "ward.D"),
                                    heatmap_layers = list(cowplot::theme_cowplot(font_size = 8),
                                                          theme(axis.title = element_blank(),
                                                                axis.text.x = element_text(angle = 45, hjust = 1))),
                                    return_ppxpy = TRUE,
                                    branches_lwd = 0.2,
                                    key.title = 'LRT\nstatistics (Λ)')

pp_pairlrt <- heatmaply:::arrange_plots(raw_pairlrt)

# add the LRT comparisons (MAIN TEXT) ========
pairlrtcast2 = pairlrtcast
pairlrtcast2[upper.tri(pairlrtcast2, diag = TRUE)] <- NA
pairlrtlong = append_group2(pairlrtcast2)

pairlrtlong %>%
    dplyr::group_by(samegroup0) %>%
    dplyr::summarise(meanlr = mean(value),
                     sdlr = sd(value))

pairlrtlong %>%
    dplyr::group_by(samegroup1) %>%
    dplyr::summarise(meanlr = mean(value),
                     sdlr = sd(value),
                     nobs = n())

pp_lrt <- ggplot(data = pairlrtlong, aes(x = samegroup1, y = value)) +
    geom_boxplot() +
    scale_y_continuous(breaks = 1:4, labels = c('10','100','1000','1e+4')) +
    ggpubr::stat_compare_means(size = textsize) +
    labs(x = 'Same Class', y = 'LRT statistics (Λ)')

pairlrtlong %>% dplyr::slice_min(value)
pairlrtlong %>% dplyr::slice_max(value)


# excluding CA26 and CL24 ========
# (pattern becomes less significant when excluding the values reaching upper boundaries)
# not included in the main text
pairlrtlong_nowolves = pairlrtlong %>%
    dplyr::filter(!(Var1 %in% c('CA26', 'CL24')), !(Var2 %in% c('CA26', 'CL24')))
pairlrtlong_nowolves %>%
    dplyr::group_by(samegroup1) %>%
    dplyr::summarise(meanlr = mean(value),
                     sdlr = sd(value),
                     nobs = n())

pp_lrt_nowolves = ggplot(data = pairlrtlong_nowolves, aes(x = samegroup1, y = value)) +
    geom_boxplot() +
    ggpubr::stat_compare_means()

# arrange allplots with the constrained grid plot ========
pp_all <- cowplot::plot_grid(plotlist = c(list(pp_Cons), allplots),
                             align = 'hv', axis = 'tblr', ncol = 3, labels = 'AUTO', label_size = labelsize) +
    draw_plot(allleg, x = 0.8, y = 0.18, width = 0.3, height = 0.1, scale = 0.8)

# arrange the pairlrt and wilcoxon plots ========
pp_pair <- cowplot::plot_grid(pp_pairlrt, pp_lrt, rel_widths = c(2,1), ncol = 2, labels = 'AUTO', label_size = labelsize)

# output files --------
ggsaver('FigSX_grid_scaled.pdf', prefix, pp_all, 18, 22)
ggsaver('FigSX_gridlr_scaled.pdf', prefix, pp_pair, 18, 10)

# save outputs
rm(grid_df)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()

