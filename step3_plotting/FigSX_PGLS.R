# Title: Plot additional PGLS results
# Author: Meixi Lin
# Date: Wed May  8 21:32:23 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(cowplot)

source("dfe_plot_util.R")
setwd(workdir)

# def functions --------
# plot correlations with the shape_us
plot_corr_shapeus <- function(PGLSl_summarydt, myvar, myvarname) {
    my.formula = y ~ x
    mysummary = PGLSl_summarydt %>%
        dplyr::filter(X_name == paste0('log10(', myvar,')'), Y_name == 'shape_us', N == 11)
    mytext = format_eq(mysummary[,c('Intercept_Value', 'Slope_Value', 'Lambda', 'Slope_pvalue')])
    pp <- ggplot(data = gammadf, aes(x = eval(as.name(myvar)), y = shape_us, ymin = shape_us_min, ymax=shape_us_max, color = pop)) +
        geom_abline(slope = mysummary$Slope_Value, intercept = mysummary$Intercept_Value,
                    linetype = 'dashed', color = 'black', size = linewidth) +
        geom_smooth(method = 'lm', formula = my.formula, se = FALSE,
                    linetype = 'dotted', color = 'darkgray', size = linewidth) +
        geom_point(size = pointsize) +
        geom_linerange(color = linerangecol) +
        scale_color_manual(values = popcols) +
        scale_x_log10() +
        # adding the labels manually so the plot is not too crowded
        # ggpmisc::stat_poly_eq(mapping = aes(x = eval(as.name(myvar)), y = Es, label = ..eq.label..),
        #                       formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'gray', vjust = -0.1) +
        # ggpmisc::stat_poly_eq(mapping = aes(x = eval(as.name(myvar)), y = Es, label = paste(stat(adj.rr.label), stat(p.value.label), sep = "~~~")),
        #                       formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'gray', vjust = 1.1) +
        # annotate('text', label = mytext, fontface = 'italic',
        #          x=min(gammadf[,myvar]), y=max(gammadf[,'Es_max']), hjust=0, vjust=1,
        #          size = textsize, color = 'darkgray') +
        labs(x = myvarname, y = 'Shape', color = 'Pop') +
        theme(legend.position = 'none')
    print(mytext)
    print(format_lm(pp))
    return(pp)
}

# def variables --------
prefix = 'FigSX_PGLS'
dir.create(prefix)

# load data --------
# load previous data
load('./Fig_TraitCor/Fig_TraitCor.RData')

# variables tested
mytraits = c("Nanc", "bodymass_gram", "gen_year", "agematurity_days", "maxlongevity_year", "mu")
mytraitnames = c("Long-term population size", "Body mass (g)", "Generation time (year)", "Age at maturity (day)", "Maximum longevity (year)", "Exon mutation rate (/bp/gen)")

names(mytraitnames) = mytraits

# main --------
# plot correlations with generation time and max longevity
pp_gen <- plot_corr(PGLSl_summarydt, myvar = 'gen_year', myvarname = mytraitnames['gen_year'])
pp_max <- plot_corr(PGLSl_summarydt, myvar = 'maxlongevity_year', myvarname = mytraitnames['maxlongevity_year'])

# plot correlations with shape_us
pp_shape <- lapply(mytraits, function(xx) {
    message(xx)
    plot_corr_shapeus(PGLSl_summarydt, myvar = xx, myvarname = mytraitnames[xx])
})


# output files --------
outpp <- cowplot::plot_grid(pp_gen, pp_max, align = 'hv', axis = 'tblr', ncol = 2, labels = 'AUTO', label_size = labelsize)
ggsaver('FigSX_PGLSEs_main.pdf', prefix, outpp, 15.5, 8)
outpp_shape <- cowplot::plot_grid(plotlist = pp_shape, align = 'hv', axis = 'tblr', ncol = 2, labels = 'AUTO', label_size = labelsize)
ggsaver('FigSX_PGLSshape_main.pdf', prefix, outpp_shape, 15.5, 21)

save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()
