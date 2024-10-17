# Title: Robustness summary
# Author: Meixi Lin
# Date: 2024-03-05 14:44:34

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(cowplot)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

# def functions --------
# compare the values
plot_compare <- function(df, comparevar, mytitle) {
    pos <- position_dodge2(width = 0.6, preserve = 'single')
    pp <- ggplot(df, aes(x = group, y = Es, color = pop, alpha = dll200, ymin = Es_min, ymax = Es_max, group = pop)) +
        geom_point(size = pointsize, position = pos) +
        geom_linerange(position = pos, color = linerangecol) +
        scale_color_manual(values = popcols, guide = "none") +
        # scale_shape_discrete(labels = c('Best', 'Medium', 'Worst')) +
        facet_wrap(as.formula(paste('~', comparevar))) +
        scale_y_log10(limits = c(4.9e-06,1.62), breaks = 10^seq(-4,0), labels = c('0.0001', '0.001', '0.01', '0.1', '1')) +
        scale_alpha_manual(values = c(alphavalue,1)) +
        labs(title = mytitle, x = '', y = 'E[|s|]', color = 'Pop', alpha ="\u0394LL < 200") +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 6)) # smaller text so it fits
    return(pp)
}

get_range <- function(dflist) {
    xx = unlist(lapply(dflist, function(df) {unlist(df[,c('Es', 'Es_min', 'Es_max')])}))
    xx = xx[xx != 0]
    return(range(xx, na.rm = TRUE))
}

rank_aic_dll <- function(df) {
    outdf = df %>%
        dplyr::group_by(pop) %>%
        dplyr::mutate(aic_rank = rank(aic)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dll200 = delta_ll < 200)
    return(outdf)
}

# def variables --------
outdir = './Fig_Robustness/'
dir.create(outdir)

# DFE compared
dfe_funcs = c('gamma', 'neugamma', 'lognormal')

# best demogs and best demogs
best_masks = apply(bestmodelsdf[,c("pop", "mask_singleton")], 1, paste, collapse=" ")
best_demogs = apply(bestmodelsdf[,c("pop", "demog_model")], 1, paste, collapse=" ")

# load data --------
# by dfe functions
dfedf = lapply(dfe_funcs, load_DFEdf, dftype = "table", filterstrings = bestmodels) %>%
    dplyr::bind_rows() %>%
    dplyr::select(pop, pdf_func, aic, Es, Es_min, Es_max, ll_data, ll_model, delta_ll) %>%
    dplyr::mutate(pdf_func = factor(pdf_func, levels = dfe_funcs))

# by demographic models (gamma dfe and best masking schemes)
gammademogdf = load_DFEdf("gamma", "table", filtercols = c("pop", "mask_singleton"),
                          filterstrings = best_masks)

# by mask singleton treatments (gamma dfe and best demographic models)
gammamaskdf = load_DFEdf("gamma", "table", filtercols = c("pop", "demog_model"),
                          filterstrings = best_demogs)

# rank aic and log-likelihood
dfedf = rank_aic_dll(dfedf)
gammademogdf = rank_aic_dll(gammademogdf)
gammamaskdf = rank_aic_dll(gammamaskdf)

# main --------
# plotting ========
dfedf = append_group(dfedf)
gammademogdf = append_group(gammademogdf)
gammamaskdf = append_group(gammamaskdf)
get_range(list(dfedf, gammademogdf, gammamaskdf))

# plot comparisons
ppA <- plot_compare(dfedf, 'pdf_func', 'Functional forms of DFE')
ppB <- plot_compare(gammademogdf, 'demog_model', 'Demographic model')
ppC <- plot_compare(gammamaskdf, 'mask_singleton', 'Mask singleton')

Cleg <- ggpubr::get_legend(ppC + theme(legend.position = 'right'))

# MAIN TEXT ========
table(dfedf[,c('pdf_func', 'aic_rank')])
dfedf$pdf_func = factor(dfedf$pdf_func, levels = dfe_funcs)

dfedf %>%
    dplyr::group_by(pdf_func) %>%
    dplyr::summarise(meanEs = mean(Es))

dfedf %>%
    dplyr::group_by(pop) %>%
    dplyr::summarise(rangeEs = max(Es)/min(Es)) %>%
    dplyr::arrange(rangeEs)

dfedf %>%
    dplyr::group_by(pdf_func, group) %>%
    dplyr::summarise(meanEs = mean(Es))

# assemble figures ========
pp <- align_plots(ppA, ppB, ppC, align = 'hv', axis = 'lr')

pp2 <- ggdraw() +
    draw_plot(pp[[1]], x = 0, y = 0.67, width = 1, height = .33) +
    draw_plot(pp[[2]], x = 0, y = 0.34, width = 1, height = .33) +
    draw_plot(pp[[3]], x = 0, y = 0.01, width = 0.72, height = .33) +
    draw_plot_label(label = c("A", "B", "C"), size = labelsize,
                    x = c(0,0,0), y = c(1,0.67,0.34))

# output files --------
ggsaver('Fig_Robustness_main.pdf', outdir, pp2, 9, 21)
ggsaver('Fig_Robustness_leg.pdf', outdir, Cleg, 2, 2)

save.image(file = paste0(outdir, 'Fig_Robustness.RData'))

# cleanup --------
date()
closeAllConnections()
