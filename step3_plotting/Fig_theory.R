# Title: Theory charts on DFE
# Author: Meixi Lin
# Date: Thu Jul 11 12:30:50 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

library(cowplot)

# def functions --------
gen_dummy <- function(Ne, maxNe) {
    myend = Ne/maxNe
    out = expand.grid(seq(0, myend, by = 0.05), seq(0, myend, by = 0.05))
    return(out)
}

plot_waffle <- function(df, mycolor) {
    pp <- ggplot(df, aes(x = Var1, y = Var2)) +
        geom_point(color = mycolor) +
        lims(x = c(0,1), y = c(0,1)) +
        theme_void()
    return(pp)
}

plot_gammadf <- function(shape, scale_us, myfill) {
    # get the data frame
    plotdf = pgamma_dfe(shape = shape, scale = scale_us)
    print(plotdf)
    plotdf$cat = factor(plotdf$cat, levels = plotdf$cat)
    pp <- ggplot(data = plotdf, aes(x = cat, y = prop)) +
        geom_col(width = 0.8, fill = myfill, color = 'black', size = linewidthcol) +
        labs(x = '|s|', y = 'Probability mass') +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank()) +
        lims(y = c(0,0.6))
    return(pp)
}

# def variables --------
outdir = './Fig_theory/'

gammadf = load_DFEdf("gamma", "table", filterstrings = bestmodels)
summary(gammadf[,c('shape', 'scale', 'Nanc')])

shape = 0.20 # median gammadf shape (0.1933)
scale_s = 5000 # median gammadf shape (4614)
Ne_true = c(1e+6, 1e+4) # typical insect / mammal popsizes
# Ne_dummy = sqrt(Ne_true)/10 # take sqrt for the length
Ne_dummy = c(10,3)
scale_us = scale_s/Ne_true

# load data --------
# dummy data
nedf = lapply(Ne_dummy, gen_dummy, maxNe = max(Ne_dummy))
lapply(nedf, nrow)

# main --------
# waffle plots
ppA1 <- plot_waffle(nedf[[1]], popcols['DM100'])
ppB1 <- plot_waffle(nedf[[2]], popcols['HS100'])

# gamma DFE
ppA2 <- plot_gammadf(shape, scale_us[1], popcols['DM100'])
ppB2 <- plot_gammadf(shape, scale_us[2], popcols['HS100'])

# output files --------
alpp <- align_plots(ppA1, ppA2, ppB1, ppB2, align = 'hv', axis = 'tblr')

pp6 <- ggdraw() +
    draw_plot(alpp[[1]], x = 0, y = 0.5, width = .3, height = .4) +
    draw_plot(alpp[[2]], x = 0.65, y = 0.5, width = .35, height = .4) +
    draw_plot(alpp[[3]], x = 0, y = 0., width = .3, height = .4) +
    draw_plot(alpp[[4]], x = 0.65, y = 0., width = .35, height = .4)
    # draw_plot_label(label = c("A", "B"), size = labelsize,
                    # x = c(0, 0), y = c(1, 0.5))

ggsaver('Fig_theory_main.pdf', outdir, pp6, 18, 11)
save.image(file = paste0(outdir, 'Fig_theory.RData'))


# cleanup --------
date()
closeAllConnections()
