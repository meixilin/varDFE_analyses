# Title: Fig lourenco DFE
# Author: Meixi Lin
# Date: 2024-03-05 14:45:43

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(cowplot)

source("dfe_plot_util.R")
setwd(workdir)

# def functions --------
summarize_bygroup <- function(df,myvar) {
    outdf = df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(mean = mean(eval(as.name(myvar))),
                         sdlr = sd(eval(as.name(myvar))),
                         nobs = n())
    return(outdf)
}

reshape_tabledf <- function(tabledf) {
    sigmadf = tabledf %>%
        dplyr::select(pop, group, starts_with('sigma_us')) %>%
        dplyr::mutate(variable = 'scale')
    mdf = tabledf %>%
        dplyr::select(pop, group, starts_with('m_us')) %>%
        dplyr::mutate(variable = 'pleiotropy')
    colnames(sigmadf)[3:5] = c('value', 'value_min', 'value_max')
    colnames(mdf)[3:5] = c('value', 'value_min', 'value_max')
    outdf = rbind(sigmadf,mdf)
    return(outdf)
}

# def variables --------
outdir = './Fig_lourenco_DFE/'
dir.create(outdir)

# load data --------
# the main table
tabledf = load_DFEdf("lourenco_eq", "table", filterstrings = bestmodels)

# the proportions 2Ns
propdf = read.delim(file = './TabSX_DFE/lourenco_eq/11pops_lourenco_eq_prop_scaled.tsv')[,-1]

# the proportions s
propdf_us = read.delim(file = './TabSX_DFE/lourenco_eq/11pops_lourenco_eq_prop_unscaled.tsv')[,-1]

# summarize the parameters by group --------
tabledf = append_group(tabledf)
summarize_bygroup(tabledf,'m')
summarize_bygroup(tabledf,'sigma')

# main --------
my.formula = y ~ x
# plot correlation of Ne_dadi with Ne ========
summary(lm(log10(Ne) ~ log10(Ne_dadi), data = tabledf))

mylims = c(min(c(tabledf$Ne_dadi, tabledf$Ne)), max(c(tabledf$Ne_dadi, tabledf$Ne)))

ppA <- ggplot(data = tabledf, aes(x = Ne_dadi, y = Ne)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted', size = linewidth, color = 'darkgray') +
    scale_x_log10(limits = mylims, breaks = c(1e+4,1e+6,1e+8)) +
    scale_y_log10(limits = mylims, breaks = c(1e+4,1e+6,1e+8)) +
    geom_smooth(method = 'lm', formula = my.formula, se = FALSE, linetype = 'dashed',
                color = 'black', size = linewidth) +
    ggpmisc::stat_poly_eq(mapping = aes(x = Ne_dadi, y = Ne, label = ..eq.label..),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'black', vjust = -1.1,
                          label.x = "right", label.y = "bottom") +
    ggpmisc::stat_poly_eq(mapping = aes(x = Ne_dadi, y = Ne, label = paste(stat(adj.rr.label), stat(p.value.label), sep = "~~~")),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'black', vjust = 0.5,
                          label.x = "right", label.y = "bottom") +
    geom_point(size = pointsize, aes(color = pop)) +
    scale_color_manual(values= popcols) +
    # guides(color = guide_legend(nrow=2,byrow=TRUE)) +
    labs(x = 'Long-term population size (Demog)', y = 'Long-term population size (DFE)', color = "") +
    theme(legend.position = 'none')

# plot pleitropy and sigma as Fig4B ========
# didn't report CI: 1. very tight CI in pleiotropy parameter 2. no CI for scale

tabledfb = reshape_tabledf(tabledf)
tabledfb$variable = factor(tabledfb$variable, labels = c('mutation pleiotropy (m)', 'mutation effect scale (σ)'))
pos <- position_dodge2(width = 0.6, preserve = 'single')

ppB <- ggplot(data = tabledfb, aes(x = group, y = value, ymin = value_min, ymax = value_max, color = pop, group = pop)) +
    geom_point(size = 1.5, position = pos) +
    # geom_linerange(position = pos, color = 'black') +
    facet_wrap(. ~ variable, nrow = 2, scales = 'free_y') +
    scale_y_log10() +
    labs(y = 'Lourenco parameters') +
    scale_color_manual(values = popcols) +
    theme(axis.title.x = element_blank(),
          legend.position = 'none')

# proportion of beneficial mutations with long term Ne ========
propdfc = propdf %>%
    dplyr::filter(cat == '[0, inf]') %>%
    dplyr::select(pop,cat, prop) %>%
    dplyr::left_join(y = tabledf[,c('pop','Ne_dadi')], by = 'pop')
propdfc$pop = factor(propdfc$pop, levels = pops)

summary(lm(prop ~ log10(Ne_dadi), propdfc))

ppC <- ggplot(data = propdfc, aes(x = Ne_dadi, y = prop)) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_log10() +
    geom_smooth(method = 'lm', formula = my.formula, se = FALSE, linetype = 'dashed',
                 color = 'black', size = linewidth) +
    geom_point(size = pointsize, aes(color = pop)) +
    ggpmisc::stat_poly_eq(mapping = aes(x = Ne_dadi, y = prop, label = ..eq.label..),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'black', vjust = -1.1, label.y = "bottom") +
    ggpmisc::stat_poly_eq(mapping = aes(x = Ne_dadi, y = prop, label = paste(stat(adj.rr.label), stat(p.value.label), sep = "~~~")),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 3, color = 'black', vjust = 0.5, label.y = "bottom") +
    scale_color_manual(values = popcols) +
    # guides(color = guide_legend(nrow=2,byrow=TRUE)) +
    labs(x = 'Long-term population size', y = 'Proportion of beneficial mutations', color = "") +
    theme(legend.position = 'none')

# plot proportions like the gamma PDF ========
propdf_us = propdf_us %>%
    dplyr::mutate(pop = factor(pop, levels = pops),
                  cat = factor(cat, levels = rev(unique(cat))))

# sanity check
propdf_us %>%
    dplyr::group_by(pop) %>%
    dplyr::summarise(ptotal = sum(prop))

ppD <- ggplot(data = propdf_us, aes(x = cat, y = prop, fill = pop)) +
    geom_col(position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
    geom_vline(xintercept = 3.5, linetype = 'dashed', size = linewidth, color = 'black') +
    labs(x = 'selection coefficient (s)', y = 'Probability mass', fill = 'Pop') +
    scale_x_discrete(labels = c("0 < s ≤ 10⁻⁵", "10⁻⁵ < s ≤ 0.0001", "s > 0.0001",
                                "-10⁻⁵ < s ≤ 0", "-0.0001 < s ≤ -10⁻⁵", "-0.001 < s ≤ -0.0001", "-0.01 < s ≤ -0.001", "s ≤ -0.01")) +
    scale_fill_manual(values= popcols) +
    theme(legend.position = 'none',
          axis.title.x = element_blank())

# combine plots --------
alpp <- cowplot::align_plots(ppA, ppB, ppC, ppD, align = 'hv', axis = 'tblr')
pp5 <- ggdraw() +
    draw_plot(alpp[[1]], x = 0, y = 0.5, width = .38, height = .5) +
    draw_plot(alpp[[2]], x = 0.38, y = 0.5, width = .24, height = .5) +
    draw_plot(alpp[[3]], x = 0.62, y = 0.5, width = .38, height = .5) +
    draw_plot(alpp[[4]], x = 0, y = 0, width = 1, height = .5) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = labelsize,
                    x = c(0, 0.38, 0.62, 0), y = c(1,1,1,0.5))

# output files --------
ggsaver('Fig_lourenco_DFE.pdf', outdir, pp5, 18, 14)

save.image(file = paste0(outdir, 'Fig_lourenco_DFE.RData'))

# cleanup --------
date()
closeAllConnections()
