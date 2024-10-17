# Title: AIC comparisons for DFE functions
# Author: Meixi Lin
# Date: Tue Apr 30 18:49:41 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

require(cowplot)
sessionInfo()

# def functions --------
format_plotdt <- function(plotdt, pdflevels) {
    plotdt$common_name = factor(plotdt$pop, levels = pops, labels = common_name_labs)
    plotdt$pdf_func = factor(plotdt$pdf_func, levels = pdflevels)
    return(plotdt)
}

plot_aic <- function(plotdt) {
    pp = ggplot(plotdt, aes(x = pdf_func, y = aic, color = pdf_func)) +
        geom_point(size = pointsize) +
        facet_wrap(. ~ common_name, scales = 'free_y', nrow = 1) +
        labs(y = 'AIC', x = 'DFE function', color = 'DFE function') +
        theme(strip.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6),
              legend.position = 'top',
              axis.text.x = element_text(angle = 90))
    return(pp)
}

plot_ll <- function(plotdt) {
    pp = ggplot(plotdt, aes(x = pdf_func, y = ll_model, color = pdf_func)) +
        geom_hline(aes(yintercept = ll_data), color = linerangecol, linetype = 'dashed') +
        geom_point(size = pointsize) +
        facet_wrap(. ~ common_name, scales = 'free_y', nrow = 1) +
        labs(y = 'Log-likelihood', x = 'DFE function', color = 'DFE function') +
        theme(strip.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6),
              legend.position = 'top',
              axis.text.x = element_text(angle = 90))
    return(pp)
}

load_DFEdf_2 <- function(dfe_func) {
    df = load_DFEdf(dfe_func, "table", filterstrings = bestmodels) %>%
        dplyr::select(pop, pdf_func, aic, ll_model, ll_data)
    return(df)
}

# def variables --------
prefix = 'FigSX_AICll_dfefunc'
dir.create(prefix)

common_name_labs = c("Mosquito","Drosophila","Halictid\nbee","Pied\nflycatcher","Collared\nflycatcher",
                     "Gray\nwolf","Arctic\nwolf","Vaquita","Fin\nwhale","Mouse","Human")

# load data --------
# load TabSX_DFE outputs
gammadf = load_DFEdf_2("gamma")
neugammadf = load_DFEdf_2("neugamma")
gammaletdf = load_DFEdf_2("gammalet")
lognormaldf = load_DFEdf_2("lognormal")

# plotdts
plotdt_gam = format_plotdt(rbind(gammadf, neugammadf, gammaletdf), c('gamma', 'neugamma', 'gammalet'))
plotdt_rob = format_plotdt(rbind(gammadf, neugammadf, lognormaldf), c('gamma', 'neugamma', 'lognormal'))

# main --------
ppaic_gam <- plot_aic(plotdt_gam)
ppaic_rob <- plot_aic(plotdt_rob)

ppll_gam <- plot_ll(plotdt_gam)
ppll_rob <- plot_ll(plotdt_rob)

# output figure
pp_gam <- cowplot::plot_grid(ppaic_gam, ppll_gam, align = 'hv', axis = 'tblr', nrow = 2, labels = "AUTO", label_size = labelsize)
pp_rob <- cowplot::plot_grid(ppaic_rob, ppll_rob, align = 'hv', axis = 'tblr', nrow = 2, labels = "AUTO", label_size = labelsize)

# output files --------
ggsaver('FigSX_AICll_gammafam.pdf', prefix, pp_gam, 18, 18)
ggsaver('FigSX_AICll_gammalog.pdf', prefix, pp_rob, 18, 18)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()
