#{{{ load required libraries, define common variables
require(tidyverse)
require(igraph)
require(grid)
require(tidyverse)
require(rlang)
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

dirg = '~/data/genome'
dirp = '~/projects/maize.grn'
dird = file.path(dirp, 'data')
dirr = file.path(dirp, 'Rmd')
#
f_cfg = file.path(dird, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
studies = t_cfg %>% distinct(study) %>% pull(study)
th = t_cfg %>% select(-tag, -mid) %>%
    filter(net_type %in% c("genotype","tissue",'tissue*genotype')) %>%
    filter(!nid %in% c("n99b_2","n99b_3","nc01","nc02")) %>%
    arrange(net_type, sample_size) %>%
    mutate(txt = sprintf("%s %s [%d]", study,note,sample_size))
thc = th %>% distinct(net_type) %>% mutate(col = pal_aaas()(3))
th = th %>% inner_join(thc, by='net_type')

load_maize_dataset <- function(id = 'me99b', opt = "exp") {
    #{{{
    if(opt == 'exp') {
        diri = '~/projects/maize.expression/data/15_output'
        fi = sprintf("%s/%s.rda", diri, id)
        stopifnot(file.exists(fi))
        x = load(fi)
        env(th = th, tm = tm)
    } else if (opt == 'grn') {
        diri = '~/projects/maize.grn/data/12_output'
        fi = sprintf("%s/%s.rda", diri, id)
        stopifnot(file.exists(fi))
        x = load(fi)
        env(rids = rids, tids = tids, reg.mat = reg.mat, tn = tn)
    } else {
        stop(sprintf("unknown dataset opt: %s", opt))
    }
    #}}}
}
read_gs <- function() {
    #{{{ read known TF targets & GO
    fi = file.path('~/projects/maize.grn/data', '09.gs.rda')
    x = load(fi)
    t_gs = t_gs %>%
        filter(! ctag %in% c("KN1_any","KN1_ear","KN1_tassel","KN1_leaf")) %>%
        mutate(ctag = ifelse(ctag=='KN1_all', 'KN1', ctag)) %>%
        mutate(binding = 1)
    t_gs %>% dplyr::count(ctag)
    t_grp_f %>% dplyr::count(ctag)
    ctags = unique(t_gs$ctag)
    t_gss = t_gs %>% distinct(ctag, reg.gid)
    env(tf = t_gs, tfs = t_gss, grp = t_grp, grp_f = t_grp_f, tf_ids = tf_ids)
    #}}}
}
read_briggs <- function() {
    #{{{ read briggs data
    diri = file.path("~/projects/briggs/data", "49_coop")
    fi = file.path(diri, "01.master.rda")
    x = load(fi)
    br = env()
    br$tissues = unique(tm$Tissue)
    br$de = tm %>% filter(silent == 0) %>%
        mutate(DE= ifelse(pDE == 'non_DE', 'non_DE',
            ifelse(abs(log2mb) < 1, 'DE1-2',
            ifelse(abs(log2mb) < 2, 'DE2-4',
            ifelse(SPE == 'non_SPE', 'DE4+', 'SPE'))))) %>%
        mutate(DEdir = ifelse(log2mb<0, 'B>M', 'B<M')) %>%
        select(Tissue, gid, DE, DEdir)
    br$de %>% count(Tissue, DE, DEdir)
    br$des = br$de %>% group_by(Tissue) %>%
        summarise(propDE = sum(DE!='non_DE')/n()) %>% ungroup()
    br
    #}}}
}
read_biomap <- function(opt = 'all') {
    #{{{ read biomap data
    bm = load_maize_dataset("me99c")
    if (opt == 'inbred') {
        sids = bm$th %>%
            separate(Genotype, c('pa1','pa2'), sep='x', fill='right', extra='merge') %>%
            filter(is.na(pa2)) %>% pull(SampleID)
        bm$th = bm$th %>% filter(SampleID %in% sids)
        bm$tm = bm$tm %>% filter(SampleID %in% sids)
        bm = bm$tm %>% inner_join(bm$th, by = 'SampleID') %>%
            select(Tissue,Genotype,gid,CPM)
    }
    bm
    #}}}
}
radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}



