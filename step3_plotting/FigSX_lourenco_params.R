# Title: Plot the lourenco exploration results on parameters
# Author: Meixi Lin
# Date: 2024-05-13 14:15:43

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

source("dfe_plot_util.R")
setwd(workdir)

sessionInfo()


# def functions --------
read_pdfdt <- function(filename, m, sigma) {
    # generate subset string
    df = read.csv(filename, row.names = 1)
    # get all combinations
    mycols = apply(expand.grid(m, sigma), 1, function(xx) {
        mysubset = paste0('X', xx[1], '_', xx[2])
        colids = grep(mysubset, colnames(df))
        return(colids)
    }) %>% as.vector

    # reformat data
    outdf = df[,c(mycols, ncol(df))] %>%
        reshape2::melt(id.vars = 'sr') %>%
        dplyr::mutate(variable = sub('X', '', variable))
    # m, sigma, Ne info
    info = reshape2::colsplit(outdf$variable, pattern = '_', names = c('m', 'sigma', 'Ne')) %>%
        dplyr::mutate_all(as.numeric)
    outdf = cbind(outdf, info)
    return(outdf)
}

# def variables --------
prefix = 'FigSX_lourenco_params'
dir.create(prefix)

pintverr

# the m and sigma values picked
mym = 0.5
mysigma = 0.01

# load data --------
# & Ne > 1e+2
pdfdt = read_pdfdt("Data_lourenco_explore/lourenco_params_pdf.csv", mym, mysigma) %>%
    dplyr::filter(Ne < 1e+9)
catdt = read.csv('Data_lourenco_explore/lourenco_params_summary.csv') %>%
    dplyr::filter(m == mym, sigma == mysigma, Ne < 1e+9)

# main --------
# select for expected mutation effects (did not plot in the end)
esdt = catdt %>%
    dplyr::select(m, sigma, Ne, e_all, e_neg, e_pos)

# fill in p_del2 where the integration did not work well (integration could be malfunctioning on both sides, use pintverr as tolerance)
# convert negative p to 0
newcatdt = catdt %>%
    dplyr::mutate(p_del2 = ifelse(p_sum > 1 + pintverr, p_del2, p_neg - p_del3 - p_del4 - p_del5 - p_neu)) %>%
    dplyr::select(m,sigma,Ne, p_del2,p_del3,p_del4,p_del5, p_neu, p_pos) %>%
    dplyr::mutate(Ne = formatC(Ne, format = 'e', digits = 0),
                  across(starts_with('p_'), ~ if_else(. < 0, 0, .)),
                  sump = rowSums(select(., starts_with("p"))))

print(summary(newcatdt$sump))

plotdt = newcatdt %>%
    dplyr::select(-sump) %>%
    # dplyr::mutate(m = paste0('m = ', m),
    #               sigma = paste0('sigma = ', sigma)) %>%
    reshape2::melt(id.vars = c('m','sigma','Ne'))

plotdt$variable = factor(plotdt$variable,
                         levels = rev(c('p_del2','p_del3','p_del4','p_del5','p_neu','p_pos')),
                         labels = c('[0,+inf)','[-1e-5,0)','[-1e-4,-1e-5)','[-1e-3,-1e-4)','[-1e-2,-1e-3)','(-inf,-1e-2)'))
# mycolors = c('lightblue', RColorBrewer::brewer.pal(n = 5, 'Reds'))

# plot the pdf
ppA <- ggplot(pdfdt, aes(x = sr, y = value, color = as.factor(Ne))) +
    geom_line(size = 0.3) +
    coord_cartesian(xlim = c(-0.01,0.01), ylim = c(0,0.05)) +
    scale_color_brewer(palette = "Reds", direction = -1) +
    labs(y = 'Probability density', x = 'selection coefficient', color = 'Ne') +
    theme(legend.position = 'none')

# plot only the representative proportions
ppB <- ggplot(plotdt, aes(x = variable, y = value, fill = as.factor(Ne))) +
    geom_col(position = position_dodge(width = 0.9), color = 'black', size = linewidthcol) +
    geom_vline(xintercept = 1.5, linetype = 'dashed', size = linewidth, color = 'black') +
    labs(x = 'selection coefficient (s)', y = 'Probability mass', fill = 'Ne') +
    scale_fill_brewer(palette = "Reds", direction = -1)

# assemble the output figure
pp <- cowplot::plot_grid(ppA, ppB, align = 'h', rel_widths = c(1,1.5), labels = 'AUTO', label_size = labelsize)


# output files --------
ggsaver('FigSX_lourenco_params.pdf', prefix ,pp, 18, 9)
# save outputs
save.image(file = paste0(prefix, '/', prefix, '.RData'))

# cleanup --------
date()
closeAllConnections()
