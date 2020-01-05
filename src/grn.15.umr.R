source("functions.R")
diri = file.path(dirg, 'Zmays_B73', 'chromatin')
opts = c('promoter','gene_body','distal','none')
tc = read_tsv(file.path(diri, 'chromatin.tsv')) %>%
    mutate(acr = factor(acr, levels=opts)) %>%
    mutate(umr = factor(umr, levels=opts))
fn = sprintf("%s/raw/rf.500k.rds", dird)
tn = readRDS(fn)
fi = file.path(dird, '17_degree', 'degree.rds')
deg = readRDS(fi)
nids= c("nc01","n13e","n14a","n15a","n19a","n17a","n18a",'n18e','n18e_1','n99a','n18g')
dirw = file.path(dird, '17_degree')

#{{{ support of ACR in leaf GRNs
tn0 = tn %>% filter(nid %in% nids) %>%
    select(nid, tn) %>% unnest() %>% select(nid, reg.gid, tgt.gid, score) %>%
    group_by(nid, tgt.gid) %>% summarize(score=sum(score)) %>% ungroup()
tn1 = tn %>% filter(nid %in% nids) %>%
    select(nid, tids) %>% unnest() %>% rename(tgt.gid = tids) %>%
    left_join(tn0, by=c('nid','tgt.gid')) %>%
    replace_na(list(score = 0))

tp = tc %>% inner_join(tn1, by=c('gid'='tgt.gid')) %>%
    inner_join(t_cfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=t_cfg$lgd))
tps = tp %>% group_by(lgd, acr) %>% summarise(scores=list(score)) %>% ungroup() %>%
    spread(acr, scores) %>%
    mutate(res = map2(promoter, none, wilcox.test, alternative='greater')) %>%
    mutate(pval = map_dbl(res, 'p.value')) %>%
    mutate(pval = str_c("p = ", scientific(pval, digits=2), sep=''))
p = ggviolin(tp, x="acr", y="score", fill="acr",
        palette = pal_npg()(9),
        add = "boxplot", add.params = list(fill = "white")) +
    geom_text(data=tps, aes(x=1.5, y=1, label=pval), hjust=.5, vjust=1, size=3) +
    scale_y_continuous(name = 'Regulation strength') +
    facet_wrap(~lgd, ncol = 3) +
    otheme(ytitle=T, xtext=T,ytext=T,xtick=T,ytick=T, legend.pos='none')
fo = file.path(dirw, 't.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

umrs = c('aUMR', 'iUMR', 'no_UMR')
tc1 = tc %>%
    mutate(umr = as.character(umr)) %>%
    mutate(umr = ifelse(umr=='promoter', ifelse(acr=='promoter', 'aUMR', 'iUMR'), 'no_UMR')) %>%
    mutate(umr = factor(umr, levels=umrs))
tp = tc1 %>% inner_join(tn1, by=c('gid'='tgt.gid')) %>%
    inner_join(t_cfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=t_cfg$lgd))
tps = tp %>%
    group_by(lgd, umr) %>% summarise(scores=list(score)) %>% ungroup() %>%
    spread(umr, scores) %>%
    mutate(res = map2(aUMR, no_UMR, wilcox.test, alternative='greater')) %>%
    mutate(pval = map_dbl(res, 'p.value')) %>%
    mutate(pval = str_c("p = ", scientific(pval, digits=2), sep=''))
p <- ggviolin(tp, x="umr", y="score", fill="umr",
        palette = pal_npg()(9),
        add = "boxplot", add.params = list(fill="white",width=.1)) +
    geom_text(data=tps, aes(x=1.5, y=1, label=pval), hjust=.5, vjust=1, size=3) +
    scale_x_discrete(name = 'Type of promoter UMR') +
    scale_y_continuous(name = 'Regulation strength') +
    facet_wrap(~lgd, ncol = 3) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T, legend.pos='none')
fo = file.path(dirw, 't2.pdf')
ggsave(p, file=fo, width=8, height=8)

tp = tc %>% inner_join(deg$tgt, by=c('gid'='tgt.gid')) %>%
    mutate(bin = ifelse(deg < 100, ifelse(deg < 50, ifelse(deg < 20,
        ifelse(deg < 10, ifelse(deg < 5, ifelse(deg == 0, 0, 1), 2), 3), 4), 5), 6)) %>%
    mutate(bin = as.character(bin)) %>%
    mutate(bin = factor(bin, levels=0:6))
#

