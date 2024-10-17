# Title: Fig1 gamma unscaled DFE
# Author: Meixi Lin
# Date: 2024-03-05 10:23:58

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(ggtree)
library(ape)
library(cowplot)
library(deeptime)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

# def functions --------

# def variables --------
outdir = './Fig_GammaDFE/'
dir.create(outdir)

# load data --------
# the main table
gammadf = load_DFEdf("gamma", "table", filterstrings = bestmodels) %>%
    dplyr::mutate(Es_rank = rank(Es, ties.method = "first"),
           shape_rank = rank(shape_us, ties.method = "first"))

# neugamma table
neugammadf = load_DFEdf("neugamma", "table", filterstrings = bestmodels) %>%
    dplyr::mutate(pneu_rank = rank(pneu, ties.method = "first"))

stopifnot(all(gammadf$pop == neugammadf$pop))
neugammadf$AIC_comp = neugammadf$aic < gammadf$aic

# the phylogenetic tree
tree11 = read.tree('../Timetree/species11.nwk')

# the proportion table
plotdf = load_DFEdf("gamma", "prop_us", filterstrings = bestmodels)

# main --------
# ggtree ========
ppA0 <- ggplot(tree11) +
    geom_tree()
ppA0 <- ppA0 %>%
    ggtree::flip(.,9,21) %>%
    ggtree::flip(.,10,11) %>%
    ggtree::flip(.,1,2) %>%
    ggtree::flip(.,15,18)
ppA0 <- revts(ppA0)
ppA <- ppA0 +
    # coord_geo(pos = "bottom", ylim = c(-0.001,11), size = 2, abbrv = TRUE, neg = TRUE, color = 'white', alpha = 0.5) +
    scale_x_continuous(breaks = c(-707,-600,-400,-200,0),
                       labels = c(707,600,400,200,0)) +
    labs(x = 'Million years ago (MYA)') +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

# plot geotime
Dcoords <- ggplot() +
    coord_geo(pos = "bottom", xlim = c(-706, 0), size = 4, abbrv = TRUE, neg = TRUE, color = 'white', alpha = 0.5) +
    theme(axis.line = element_blank())

# Es gamma ========
ppB <- ggplot(data = gammadf, aes(x = pop, y = Es, ymin = Es_min, ymax = Es_max, color = pop, label = Es_rank)) +
    scale_x_discrete(labels = common_names) +
    scale_y_log10(limits = c(0.0001,1), breaks = 10^seq(-4,0), labels = c('0.0001', '0.001', '0.01', '0.1', '1')) +
    geom_point(size = pointsize) +
    geom_linerange(color = linerangecol) +
    geom_text(size = textsize, nudge_x = 0.4) +
    scale_color_manual(values= popcols) +
    labs(x = '', y = 'E[|s|]') +
    coord_flip() +
    theme(legend.position = 'none')

# shape gamma ========
ppC <- ggplot(data = gammadf, aes(x = pop, y = shape, ymin = shape_min, ymax = shape_max, color = pop, label = shape_rank)) +
    geom_point(size = pointsize) +
    geom_linerange(color = linerangecol) +
    geom_text(size = textsize, nudge_x = 0.4) +
    scale_color_manual(values= popcols) +
    labs(x = '', y = 'Shape') +
    coord_flip() +
    theme(legend.position = 'none',
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

# proportions ========
ppD <- ggplot(data = plotdf, aes(x = cat, y = prop, ymin = prop_min, ymax=prop_max, fill = pop)) +
    geom_col(position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
    geom_linerange(position = position_dodge(width = 0.9), color = linerangecol) +
    labs(x = '|s|', y = 'Probability mass', fill = 'Pop') +
    scale_x_discrete(labels = c("-10⁻⁵ < s ≤ 0", "-0.0001 < s ≤ -10⁻⁵", "-0.001 < s ≤ -0.0001", "-0.01 < s ≤ -0.001", "s ≤ -0.01")) +
    scale_fill_manual(values= popcols) +
    theme(legend.position = 'none',
          axis.title.x = element_blank())

# pneu gamma ========
ppE <- ggplot(data = neugammadf, aes(x = 1, y = pneu, ymin = pneu_min, ymax = pneu_max,
                                     fill = pop, group = pop, alpha = AIC_comp)) +
    geom_col(position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
    geom_linerange(position = position_dodge(width = 0.9), color = linerangecol) +
    scale_fill_manual(values= popcols) +
    scale_alpha_manual(values = c(alphavalue,1), labels = c('gamma','neugamma')) +
    labs(x = '|s| = 0', y = 'Probability point mass', fill = 'Pop', alpha = "AIC") +
    guides(fill = 'none') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none')
Eleg <- ggpubr::get_legend(ppE + theme(legend.position = 'right'))

# combine plots AB ========
upp <- align_plots(ppA, ppB, ppC, align = 'h', axis = 'tb', greedy = FALSE)
lowp <- align_plots(ppD, ppE, align = 'hv', axis = 'tb')
pp1 <- ggdraw() +
    draw_plot(upp[[1]], x = 0, y = 0.45, width = .3, height = .55) +
    draw_plot(upp[[2]], x = 0.3, y = 0.45, width = .42, height = .55) +
    draw_plot(upp[[3]], x = 0.72, y = 0.45, width = .28, height = .55) +
    draw_plot(lowp[[1]], x = 0, y = 0, width = .78, height = .45) +
    draw_plot(lowp[[2]], x = 0.80, y = 0, width = .2, height = .45) +
    draw_plot_label(label = c("A", "B", "C", "D", "E"), size = labelsize,
                    x = c(0, 0.35, 0.71, 0, 0.80), y = c(1, 1, 1, 0.45, 0.45))

# output files --------
ggsaver('Fig_GammaDFE_main.pdf',outdir,pp1,18,18)
ggsaver('phylo_tree_coord.pdf',outdir,Dcoords,10,1)
ggsaver('Fig_GammaDFE_Eleg.pdf',outdir,Eleg,2,2)
save.image(file = paste0(outdir, 'Fig_GammaDFE.RData'))

# for pictures --------
# collection link
# # https://www.phylopic.org/permalinks/fee10159973987a9db2ac400b8d8d265d823ede212bce0141f51a1e55d37ca1c
# wget https://images.phylopic.org/images/036b96de-4bca-408e-adb6-2154fcd724ef/vector.svg -O human.svg
# wget https://images.phylopic.org/images/426c4d27-9898-4dbb-ab08-d1f242038f5d/vector.svg -O flycatcher.svg
# wget https://images.phylopic.org/images/1bbd9d5c-99bb-44f0-ab91-b04d523d4d34/vector.svg -O bee.svg
# wget https://images.phylopic.org/images/c8f71c27-71db-4b34-ac2d-e97fea8762cf/vector.svg -O mouse.svg
# wget https://images.phylopic.org/images/f538aa99-5c08-4f96-97d9-2e094ef5d84f/vector.svg -O mosquito.svg
# wget https://images.phylopic.org/images/ea8fa530-d856-4423-a81d-a74342cd1875/vector.svg -O drosophila.svg
# wget https://images.phylopic.org/images/5036f260-5a0d-42c5-a0bd-9eb0729e54e0/vector.svg -O wolf.svg
# wget https://images.phylopic.org/images/a23658b6-e872-45a8-b6ce-7eaba9e12a32/vector.svg -O vaquita.svg
# wget https://images.phylopic.org/images/8366482f-d41c-4ab2-8403-3d70fb04ad5f/vector.svg -O finwhale.svg

# cleanup --------
date()
closeAllConnections()
