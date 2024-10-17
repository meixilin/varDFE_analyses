# Title: Plot PGLS log-likelihood surface
# Author: Meixi Lin
# Date: Tue Oct  8 22:42:45 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(cowplot)

source("dfe_plot_util.R")
setwd(workdir)

# def functions --------
build_plotdf <- function(ii, PGLSl_my_full, PGLSldt, PGLSldt_my) {
    plotdf = data.frame(lambda = seq(0,1,0.001),
                        logL = PGLSl_my_full[[ii]]$ll_lambdas)
    pt1 = data.frame(lambda = PGLSldt[ii, 'lambda'], logL = PGLSldt[ii, 'logL'])
    pt2 = data.frame(lambda = PGLSldt_my[ii, 'lambda'], logL = PGLSldt_my[ii, 'logL'])
    pts = rbind(PGLSldt[ii,c('lambda', 'logL', 'method')],
                PGLSldt_my[ii,c('lambda', 'logL', 'method')])
    return(list(plotdf, pts))
}

# def variables --------
prefix = 'FigSX_PGLSllsurface'
dir.create(prefix)

# load data --------
# load previous data
mmutil::loadSomeRData(c('PGLSl_my_full', 'PGLSldt', 'PGLSldt_my'), './TabSX_PGLS/TabSX_PGLS.RData')

# variables tested
mytraits = c("Nanc", "bodymass_gram", "gen_year", "agematurity_days", "maxlongevity_year", "mu")
mytraitnames = c("Long-term population size", "Body mass (g)", "Generation time (year)", "Age at maturity (day)", "Maximum longevity (year)", "Exon mutation rate (/bp/gen)")

names(mytraitnames) = mytraits

# main --------
# plot the log-likelihood for only the Nanc model
mydf <- build_plotdf(1, PGLSl_my_full, PGLSldt, PGLSldt_my)
# plot log-likelihood surface with outputs
pp <- ggplot(data = mydf[[1]], mapping = aes(x = lambda, y = logL)) +
    geom_line() +
    geom_point(data = mydf[[2]], aes(color = method)) +
    scale_color_manual(values = c('red', 'gray'),
                         labels = c('Global optimum', 'Local optimum'),
                         name = '') +
    labs(x = 'Phylogenetic signal in PGLS (λ)', y = 'Log-likelihood',
         title = 'Log-likelihood surface for estimating λ in PGLSλ',
         subtitle = ' model = log10(E[|s|]) ~ log10(Nₐ), correlation = corPagel(λ)')

# output files --------
ggsaver('FigSX_PGLSllsurface.pdf', prefix, pp, 9, 9)
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()
