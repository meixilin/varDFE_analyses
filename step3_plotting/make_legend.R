rm(list = ls())
source("dfe_plot_util.R")
setwd(workdir)

# horizontal version

mypops = c(pops[1:5], "", pops[6:11])
mycols = c(popcols[1:5], "#ffffff", popcols[6:11])
names(mycols)[6] = ""
mylabs = c(common_names[1:5], "", common_names[6:11])
pp <- ggplot(data = data.frame(x = factor(mypops, levels = mypops), y = 1:12),
             mapping = aes(x = x, y = y, fill = x)) +
    geom_col() +
    scale_fill_manual(values = mycols,
                      labels = mylabs) +
    theme(legend.position = 'top',
          legend.text = element_text(size = 6.5),
          legend.spacing.x = unit(0.05, "cm"),
          legend.spacing.y = unit(0.01, "cm")) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE,title = NULL))

ppleg = ggpubr::as_ggplot(ggpubr::get_legend(pp))
ggsaver('poplegH.pdf', "./", ppleg, 11, 2)

# vertical version

pp2 <- ggplot(data = data.frame(x = factor(pops, levels = pops), y = 1:11),
             mapping = aes(x = x, y = y, fill = x)) +
    geom_col() +
    scale_fill_manual(values = popcols,
                      labels = common_names) +
    guides(fill=guide_legend(title = NULL, byrow = TRUE)) +
    theme(legend.position = 'right',
          legend.text = element_text(size = 6.5),
          legend.spacing.y = unit(0.05, "cm"))


pp2leg = ggpubr::as_ggplot(ggpubr::get_legend(pp2))
ggsaver('poplegV.pdf', "./", pp2leg, 3, 8)
