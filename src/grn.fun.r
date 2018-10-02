#{{{ load required libraries, define common variables
require(tidyverse)
require(igraph)
require(grid)
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
#}}}

dirg = file.path(Sys.getenv("genome"), "Zmays_v4/TF")
dirp = '~/projects/maize.grn'
dird = file.path(dirp, 'data')
dirr = file.path(dirp, 'Rmd')


radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

auc_barplot <- function(ti, fo, wd, ht) { 
    #{{{ ti w. 5 cols: nid tag ctag aupr auroc
    tp = ti %>% #filter(ctag %in% ctags, nid %in% nids) %>%
        rename(AUPR = aupr, AUROC = auroc) %>%
        gather(type, auc, -nid, -tag, -ctag) %>%
        mutate(auc = as.numeric(auc)) %>% filter(!is.na(auc)) %>%
        mutate(type = factor(type, levels = c("AUROC", "AUPR"))) %>%
        mutate(lab = str_remove(sprintf("%.03f", auc), '^0+'))
    tps = tp %>% distinct(nid, tag) %>% rename(lab = tag)
    tpl = tp %>% distinct(type, ctag) %>% filter(type == 'AUROC') %>%
        mutate(itc.y = .5)
    p1 = ggplot(tp) +
        geom_bar(mapping = aes(x = nid, y = auc, fill = type), stat = 'identity', width = .75, alpha = .7) +
        geom_text(mapping = aes(x = nid, y = auc , label = lab), hjust = 1, size = 2) +
        geom_hline(data = tpl, aes(yintercept = itc.y), size = .3, alpha = .5) +
        scale_x_discrete(breaks = tps$nid, labels = tps$lab, expand = c(0,0)) +
        scale_y_continuous(expand = expand_scale(mult=c(0,.05))) +
        scale_fill_npg() +
        coord_flip() +
        facet_wrap(type~ctag, scale = 'free_x', nrow = 2) +
        theme_bw() +
        theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,.2,0,'lines'))) + 
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) +
        theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
        theme(axis.title = element_blank()) +
        theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
    ggarrange(p1, nrow = 1, ncol = 1, labels = '', heights = c(2,2)) %>% 
        ggexport(filename = fo, width = wd, height = ht)
    #}}}
}


