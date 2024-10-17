# Title: Plot correlations of the m and sigma parameters
# Author: Meixi Lin
# Date: Fri May 10 20:42:50 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()

# def functions --------

# def variables --------
prefix = 'FigSX_lourenco_msig'
dir.create(prefix)

# load data --------
tabledf = load_DFEdf("lourenco_eq", "table", filterstrings = bestmodels)

# main --------
my.formula = y ~ x

# correlation MAIN TEXT ========
summary(lm(m ~ sigma, data = tabledf)) # not very linear
stats::cor(tabledf[,c("m", "sigma")], method = "spearman") # but the rank-based correlation is very clear
summary(lm(log10(m) ~ log10(sigma), data = tabledf)) # linear correlation when both are on the log10 scale

# plot the values ========
pp <- ggplot(data = tabledf, aes(x = sigma, y = m)) +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = 'lm', formula = my.formula, se = FALSE, linetype = 'dashed', color = 'black', size = linewidth) +
    ggpmisc::stat_poly_eq(mapping = aes(x = sigma, y = m, label = ..eq.label..),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 2.5, color = 'black', vjust = 0,
                          label.x = "right", label.y = "top") +
    ggpmisc::stat_poly_eq(mapping = aes(x = sigma, y = m, label = paste(stat(adj.rr.label), stat(p.value.label), sep = "~~~")),
                          formula = my.formula, parse = TRUE, p.digits = 2, size = 2.5, color = 'black', vjust = 1,
                          label.x = "right", label.y = "top") +
    geom_point(size = pointsize, aes(color = pop)) +
    scale_color_manual(values= popcols, labels = common_names) +
    labs(x = 'Mutation effect scale (sigma)', y = 'Mutation pleiotropy (m)', color = "Population") +
    theme(legend.position = 'right')

# output files --------
ggsaver('FigSX_lourenco_msig.pdf', prefix, pp, 9, 6)

# save outputs
save.image(file = paste0(prefix, '/', prefix, '.RData'))


# cleanup --------
date()
closeAllConnections()
