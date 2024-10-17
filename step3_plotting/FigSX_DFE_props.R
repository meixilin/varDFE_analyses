# Title: Plot mutation proportions for other treatments and DFE functions
# Author: Meixi Lin
# Date: Wed May  1 15:49:51 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

require(cowplot)
sessionInfo()

# def functions --------
plot_props <- function(plotdf, wrapvar, mytitle) {
    pp <- ggplot() +
        geom_col(data = plotdf, mapping = aes(x = cat, y = prop, fill = pop),
                 position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
        geom_linerange(data = plotdf, mapping = aes(x = cat, y = prop, ymin = prop_min, ymax=prop_max, group = pop),
                       position = position_dodge(width = 0.9), color = linerangecol) +
        facet_wrap(as.formula(paste('~', wrapvar)), ncol = 1) +
        labs(x = '|s|', y = 'Probability mass', fill = 'Pop', title = mytitle) +
        scale_fill_manual(values= popcols, labels = common_names) +
        theme(legend.position = 'right')
    return(pp)
}

add_pneu_plet <- function(pp, pneudf, pletdf) {
    outpp <- pp +
        geom_col(data = pneudf, mapping = aes(x = 1, y = pneu_us, group = pop),
                 position = position_dodge(width = 0.9), color = 'black', size = linewidthcol, alpha = 0.6) +
        geom_col(data = pletdf, mapping = aes(x = 5, y = plet_us, group = pop),
                 position = position_dodge(width = 0.9), color = 'black', size = linewidthcol, alpha = 0.6)
    return(outpp)
}

# def variables --------
prefix = 'FigSX_DFE_props'
dir.create(prefix)

# load data --------
# by dfe functions
gammadf = load_DFEdf("gamma", "prop_us", filterstrings = bestmodels)
neugammadf = load_DFEdf("neugamma", "prop_us", filterstrings = bestmodels)
gammaletdf = load_DFEdf("gammalet", "prop_us", filterstrings = bestmodels)
lognormaldf = load_DFEdf("lognormal", "prop_us", filterstrings = bestmodels)

# add pneu and plet proportions
pneudf = load_DFEdf("neugamma", "table", filterstrings = bestmodels) %>%
    dplyr::select(pop, pdf_func, starts_with('pneu_us'))
pletdf = load_DFEdf("gammalet", "table", filterstrings = bestmodels) %>%
    dplyr::select(pop, pdf_func, starts_with('plet_us'))

# by demographic models (gamma dfe and best masking schemes)
gammademogdf = load_DFEdf("gamma", "prop_us", filtercols = c("pop", "mask_singleton"),
                          filterstrings = apply(bestmodelsdf[,c("pop", "mask_singleton")], 1, paste, collapse=" "))

# by mask singleton treatments (gamma dfe and best demographic models)
gammamaskdf = load_DFEdf("gamma", "prop_us", filtercols = c("pop", "demog_model"),
                          filterstrings = apply(bestmodelsdf[,c("pop", "demog_model")], 1, paste, collapse=" "))

gammafamdf = rbind(gammadf, neugammadf, gammaletdf) %>%
    dplyr::mutate(pdf_func = factor(pdf_func, levels = c('gamma', 'neugamma', 'gammalet')))

gammalogdf = rbind(gammadf, neugammadf, lognormaldf) %>%
    dplyr::mutate(pdf_func = factor(pdf_func, levels = c('gamma', 'neugamma', 'lognormal')))

# main --------
pp_gam0 <- plot_props(gammafamdf, 'pdf_func', 'Gamma DFE variations')
pp_gammafam <- add_pneu_plet(pp_gam0, pneudf, pletdf)
pp_gammalog <- plot_props(gammalogdf, 'pdf_func', 'Functional forms of DFE')
pp_gammademog <- plot_props(gammademogdf, 'demog_model', 'Demographic model')
pp_gammamask <- plot_props(gammamaskdf, 'mask_singleton', 'Mask singleton')

# output files --------
ggsaver('FigSX_prop_gammafam.pdf', prefix, pp_gammafam, 16, 18)
ggsaver('FigSX_prop_gammalog.pdf', prefix, pp_gammalog, 16, 18)
ggsaver('FigSX_prop_gammademog.pdf', prefix, pp_gammademog, 16, 18)
ggsaver('FigSX_prop_gammamask.pdf', prefix, pp_gammamask, 16, 12)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()


