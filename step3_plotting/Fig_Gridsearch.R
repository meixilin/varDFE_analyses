# Title: Grid search output
# Author: Meixi Lin
# Date: 2024-03-05 14:45:36

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(data.table)
library(raster)
library(cowplot)
library(heatmaply)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

# def functions --------
# generate grid0
generate_grid0 <- function(grid_df, mypops, lrtdf) {
    # message(myscale)
    grid0 = add_grids(grid_df[mypops], myscale)
    # append delta_ll for grid0
    grid0$delta_ll = lrtdf$ll_data - grid0$ll_model
    return(grid0)
}

# aggregate grid0 for plotting
aggregate_grid0 <- function(grid0, myfac = 10) {
    # message(myscale)
    myvars = c(rev(scale_vars(myscale)), 'delta_ll')
    rs = raster::rasterFromXYZ(grid0[,..myvars])
    outrs = raster::aggregate(rs, fact = myfac, fun = mean)
    outgrid0 = raster::rasterToPoints(outrs) %>%
        data.table::data.table(.)
    setnames(outgrid0, myvars)
    return(outgrid0)
}

plot_grids <- function(mygrid, mledf, ll_grid0, grid0col, mytitle) {
    mledf$pop = factor(mledf$pop, levels = pops)
    stopifnot(all(rev(scale_vars(myscale)) == colnames(mygrid)[1:2]))
    xname = colnames(mygrid)[1]
    yname = colnames(mygrid)[2]
    message(paste0('x = ', xname, ' y = ', yname))
    message('grid log10(delta_ll) summary: ')
    print(summary(log10(mygrid[['delta_ll']])))
    pp <- ggplot() +
        geom_tile(data = mygrid, aes_string(x = xname, y = yname, fill = 'delta_ll')) +
        geom_point(data = mledf, shape = 19, size = pointsize, # stroke = 2,
                   aes_string(x = xname, y = yname, color = 'pop')) +
        geom_point(data = ll_grid0, shape = 8, size = pointsize, color = grid0col,
                   aes_string(x = xname, y = yname)) +
        scale_color_manual(values = popcols) +
        guides(color = 'none') +
        labs(x = 'Scale', y = 'Shape',fill = '∆LL', subtitle = mytitle) +
        scale_fill_gradientn(colors = heat.colors(n = 256),
                             breaks = 10^(c(0:6)),
                             limits = 10^c(0,6.252),
                             labels = c("1","10","100", "", "1e+4", "", "1e+6"),
                             trans = 'pseudo_log') +
        theme(aspect.ratio = 1,
              legend.key.width = unit(0.5, 'cm'),
              legend.position = 'none')
    return(pp)
}

# note that the mledf here is the results from gridsearch not the DFE inference
plot_grids_pop <- function(mypop) {
    stopifnot(!is.na(common_names[mypop]))
    print(mledf[pop == mypop]) # cannot use `pop` as function input as it confuses data table
    pp <- plot_grids(mygrid = aggregate_grid0(grid_df[[mypop]]),
                     mledf = mledf,
                     ll_grid0 = mledf[pop == mypop], grid0col = "black",
                     mytitle = unname(common_names[mypop]))
    pp$layers[[2]] <- NULL
    return(pp)
}

append_group2 <- function(df) {
    outdf = reshape2::melt(df) %>%
        dplyr::mutate(group1 = case_when(Var1 %in% popinsects ~ 'Insecta',
                                         Var1 %in% popbirds ~ 'Aves',
                                         TRUE ~ 'Mammalia'),
                      group2 = case_when(Var2 %in% popinsects ~ 'Insecta',
                                         Var2 %in% popbirds ~ 'Aves',
                                         TRUE ~ 'Mammalia')) %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::mutate(samegroup0 = paste0(group1, group2),
                      samegroup1 = group1 == group2)
    return(outdf)
}

# def variables --------
outdir = './Fig_Gridsearch/'
dir.create(outdir)
myscale = 'unscaled'

# load data --------
grid_df = readRDS(file = 'TabSX_gridsearch/unscaled/11pops_gamma_unscaled_FREEgrid.rds')
loadSomeRData(c('scale_vars','check_grids','add_grids','pairlrtcast', 'mledf', 'lrt_11pops'),
                           'TabSX_gridsearch/unscaled/TabSX_gridsearch_unscaled_nogrid.RData')

# name common_names
names(common_names) = pops
common_names

# main --------
grid0agg_11pops = aggregate_grid0(generate_grid0(grid_df, pops, lrt_11pops))

pp_Cons <- plot_grids(mygrid = grid0agg_11pops,
                      mledf = mledf,
                      ll_grid0 = lrt_11pops, grid0col = "black",
                      mytitle = 'Constrained')

pp_BP44 <- plot_grids_pop(mypop = 'BP44')

pp_DM100 <- plot_grids_pop(mypop = 'DM100')

ABCleg <- ggpubr::get_legend(pp_Cons, position = 'top')

# add a pretty heatmap ========
# https://stats.stackexchange.com/questions/195446/choosing-the-right-linkage-method-for-hierarchical-clustering/217742#217742
# picked the ward.D method so that the higher-order break points are not "singular"
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
                                    branches_lwd = 0.2)

Dleg <- ggpubr::get_legend(raw_pairlrt$p, position = 'right')
pp_pairlrt <- heatmaply:::arrange_plots(raw_pairlrt, hide_colorbar = TRUE)

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

# excluding CA26 and CL24 (pattern becomes less significant when excluding the values reaching upper boundaries)
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

pairlrtlong %>% dplyr::slice_min(value)
pairlrtlong %>% dplyr::slice_max(value)

# combine plots ========
topp <- cowplot::plot_grid(pp_Cons, pp_BP44, pp_DM100, align = 'hv', axis = 'tblr', nrow = 1)

pp3 <- ggdraw() +
    draw_plot(topp, x = 0, y = 0.6, width = 1, height = 0.4) +
    draw_plot(pp_pairlrt, x = 0.03, y = 0, width = 0.5, height = 0.5) +
    draw_plot(pp_lrt, x = 0.67, y = 0, width = 0.33, height = 0.5) +
    draw_plot_label(label = c("A", "B", "C", "D","E"),
                    x = c(0, 0.34, 0.67, 0, 0.67), y = c(1,1,1,0.5,0.5),
                    size = labelsize)

# output files --------
ggsaver('Fig_Gridsearch_main.pdf', outdir, pp3, 18, 16)
ggsaver('Fig_Gridsearch_ABCleg.pdf', outdir, ABCleg, 6, 2)
ggsaver('Fig_Gridsearch_Dleg.pdf', outdir, Dleg, 2, 4)
# remove huge files
rm(grid_df)
rm(grid0agg_11pops)

save.image(file = paste0(outdir, 'Fig_Gridsearch.RData'))

# cleanup --------
date()
closeAllConnections()
