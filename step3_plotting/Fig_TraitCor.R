# Title: Plot the correlation of Es with other traits
# Author: Meixi Lin
# Date: 2024-03-05 14:44:43

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

library(ggpubr)
library(cowplot)

sessionInfo()

# def functions --------
myformatter <- function(a,mydigits=2) {
    if(abs(a) < 1e-3) {
        out = format(a,digits=mydigits,scientific = T)
    } else {
        out = round(a, digits=mydigits)
    }
    return(out)
}

format_eq <- function(ablp) {
    aa = myformatter(ablp[1])
    bb = myformatter(ablp[2])
    ll = myformatter(ablp[3])
    pp = myformatter(ablp[4], mydigits=3)
    if (bb > 0) {bb = paste0('+',bb)}
    out = paste0('PGLS: y = ', aa, bb, 'x;', ' λ = ', ll,', P = ', pp)
    return(out)
}

format_lm <- function(myplot) {
    data <- ggplot_build(myplot)$data[[3]]
    lmres <- summary(lm(formula = y ~ x, data = data))
    aa = myformatter(lmres$coefficients[1,1])
    bb = myformatter(lmres$coefficients[2,1])
    r2 = myformatter(lmres$adj.r.squared)
    pp = myformatter(lmres$coefficients[2,4], mydigits=3)
    if (bb > 0) {bb = paste0('+',bb)}
    out = paste0('LM: y = ', aa, bb, 'x;', ' R2adj = ', r2,', P = ', pp)
    return(out)
}

plot_corr <- function(PGLSl_summarydt, myvar, myvarname) {
    my.formula = y ~ x
    mysummary = PGLSl_summarydt %>%
        dplyr::filter(X_name == paste0('log10(', myvar,')'), Y_name == 'log10(Es)', N == 11)
    mytext = format_eq(mysummary[,c('Intercept_Value', 'Slope_Value', 'Lambda', 'Slope_pvalue')])
    pp <- ggplot(data = gammadf, aes(x = eval(as.name(myvar)), y = Es, ymin = Es_min, ymax=Es_max, color = pop)) +
        geom_abline(slope = mysummary$Slope_Value, intercept = mysummary$Intercept_Value,
                    linetype = 'dashed', color = 'black', size = linewidth) +
        geom_smooth(method = 'lm', formula = my.formula, se = FALSE,
                    linetype = 'dotted', color = 'darkgray', size = linewidth) +
        geom_point(size = pointsize) +
        geom_linerange(color = linerangecol) +
        scale_color_manual(values = popcols) +
        scale_x_log10() +
        scale_y_log10(limits = c(0.0001,1), breaks = 10^seq(-4,0), labels = c('0.0001', '0.001', '0.01', '0.1', '1')) +
        # adding the labels manually so the plot is not too crowded
        # ggpmisc::stat_poly_eq(mapping = aes(x = eval(as.name(myvar)), y = Es, label = ..eq.label..),
        #                       formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'gray', vjust = -0.1) +
        # ggpmisc::stat_poly_eq(mapping = aes(x = eval(as.name(myvar)), y = Es, label = paste(stat(adj.rr.label), stat(p.value.label), sep = "~~~")),
        #                       formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'gray', vjust = 1.1) +
        # annotate('text', label = mytext, fontface = 'italic',
        #          x=min(gammadf[,myvar]), y=max(gammadf[,'Es_max']), hjust=0, vjust=1,
        #          size = textsize, color = 'darkgray') +
        labs(x = myvarname, y = 'E[|s|]', color = 'Pop') +
        theme(legend.position = 'none')
    message(mytext)
    message(format_lm(pp))
    return(pp)
}

# get legend for two models
model_leg <- function() {
    df = data.frame(a = 1:4, b = 1:4, Model = factor(c('PGLS', 'LM', 'PGLS', 'LM'), levels = c('PGLS', 'LM')))
    pp = ggplot(data = df, aes(x = a, y = b, linetype = Model, color = Model)) +
        geom_line() +
        scale_linetype_manual(values = c('dashed', 'dotted')) +
        scale_color_manual(values = c('black', 'darkgray')) +
        theme(legend.key.width = unit(1, 'cm'))
    leg = ggpubr::get_legend(pp)
    return(leg)
}

# def variables --------
outdir = './Fig_TraitCor/'
dir.create(outdir)

# load data --------
loadSomeRData(vars = c('gammadf', 'PGLSl_summarydt'), file = './TabSX_PGLS/TabSX_PGLS.RData')

# gammadf with info (a lot of edits added before, should not load from scratch)
gammadf$pop = factor(gammadf$pop, levels = gammadf$pop)
# get range
range(unlist(gammadf[,c('Es', 'Es_min', 'Es_max')]), na.rm = TRUE)

gammadf = gammadf %>%
    dplyr::mutate(E2Nas_rank = rank(E2Nas, ties.method = "first"))

# the proportion table
plotdf = load_DFEdf("gamma", "prop", filterstrings = bestmodels)

# main --------
# correlation with Nanc
ppA <- plot_corr(PGLSl_summarydt, myvar = 'Nanc', myvarname = 'Long-term population size')
# with body mass
ppB <- plot_corr(PGLSl_summarydt, myvar = 'bodymass_gram', myvarname = 'Body mass (g)')
# with age at maturity
ppC <- plot_corr(PGLSl_summarydt, myvar = 'agematurity_days', myvarname = 'Age at maturity (day)')
# with mutation rates
ppD <- plot_corr(PGLSl_summarydt, myvar = 'mu', myvarname = 'Exon mutation rate (/bp/gen)')

# population size based E2Nas models ========
ppE <- ggplot(data = gammadf, aes(x = pop, y = E2Nas, ymin = E2Nas_min, ymax = E2Nas_max, color = pop, label = E2Nas_rank)) +
    # scale_x_discrete(labels = common_names) +
    scale_y_log10() +
    geom_point(size = pointsize) +
    geom_linerange(color = linerangecol) +
    geom_text(size = textsize, nudge_y = 0.3) +
    scale_color_manual(values= popcols) +
    # scale_x_discrete(labels = common_names) +
    labs(x = 'Populations', y = expression("E[|2"*N[a]*"s|]")) +
    theme(legend.position = 'none',
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
          # axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# and proportions ========
ppF <- ggplot(data = plotdf, aes(x = cat, y = prop, ymin = prop_min, ymax=prop_max, fill = pop)) +
    geom_col(position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
    geom_linerange(position = position_dodge(width = 0.9), color = linerangecol) +
    scale_x_discrete(labels = c("-1 < 2Nₐs ≤ 0", "-10 < 2Nₐs ≤ -1", "-100 < 2Nₐs ≤ -10", "-1000 < 2Nₐs ≤ -100", "2Nₐs ≤ -1000")) +
    labs(x = expression("|2"*N[a]*"s|"), y = 'Probability mass', fill = 'Pop') +
    scale_fill_manual(values= popcols) +
    theme(legend.position = 'none',
          axis.title.x = element_blank())

# plot a legend ========
leg <- model_leg()

# arrange figures --------
alpp <- align_plots(ppA, ppB, ppC, ppD, ppE, ppF, align = 'hv', axis = 'tblr')

pp4 <- ggdraw() +
    draw_plot(alpp[[1]], x = 0, y = 0.7, width = .4, height = .3) +
    draw_plot(alpp[[2]], x = 0.45, y = 0.7, width = .4, height = .3) +
    draw_plot(alpp[[3]], x = 0, y = 0.4, width = .4, height = .3) +
    draw_plot(alpp[[4]], x = 0.45, y = 0.4, width = .4, height = .3) +
    draw_plot(alpp[[5]], x = 0, y = 0, width = .3, height = .4) +
    draw_plot(alpp[[6]], x = 0.3, y = 0, width = .7, height = .4) +
    draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = labelsize,
                    x = c(0,0.45,0,0.45,0,0.35), y = c(1,1,0.7,0.7,0.4,0.4))

# output -------
ggsaver('Fig_TraitCor_main.pdf', outdir, pp4, 18, 18)
ggsaver('Fig_TraitCor_leg.pdf', outdir, leg, 2, 2)
save.image(file = paste0(outdir, 'Fig_TraitCor.RData'))

# cleanup --------
date()
closeAllConnections()
