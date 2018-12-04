source("functions.R")
dirw = file.path(dird, '14_eval_sum')

gs = read_gs()
#br = read_briggs()
#bm = read_biomap()

eval_pcc_briggs <- function(rids, tids, tn, tm3) {
    #{{{ eval functions
    to = tn %>%
        filter(row_number() <= net_size) %>%
        filter(reg.gid %in% gids, tgt.gid %in% gids) %>%
        mutate(pcc = map2_dbl(reg.gid, tgt.gid, function(i,j,mat) mat[i,j], mat = tm3))
    to$pcc
    #}}}
}


#{{{ # evaluate PCC
if(FALSE) {
ds = load_maize_dataset(id = 'nc03', opt = 'grn')

tissues = unique(br$th$Tissue)
net_size = 10000
tissue = tissues[16]
tm1 = br$tm %>%
    select(-FPKM) %>%
    inner_join(br$th, by = 'SampleID') %>%
    filter(Tissue == tissue) %>%
    select(Genotype, gid, CPM)
gids = tm1 %>% group_by(gid) %>%
    summarise(n.exp = sum(CPM >= 1)) %>% ungroup() %>%
    filter(n.exp >= 2) %>% pull(gid)
length(gids)
tm2 = tm1 %>% filter(gid %in% gids) %>% spread(Genotype, CPM)
em = t(as.matrix(tm2[,-1]))
colnames(em) = tm2$gid
tm3 = cor(em, method = 'pearson')

nids = c('nc01','nc03')
y = map(nids1, eval_pcc_briggs, tm3 = tm3, rids = rids, tids = tids)
to = tibble(nid = rep(nids, map_int(y, length)), pcc = flattern_dbl(y))
}
#}}}

#ev_tf = th %>% mutate(fi = sprintf("%s/13_eval/%s_tf.rds", dird, nid)) %>%
    #mutate(res = map(fi, readRDS)) %>%
    #mutate(roc = map(res, "roc"), pr = map(res, "pr"),
           #auroc = map(res, "auroc"), aupr = map(res, "aupr")) %>%
    #select(nid,roc,pr,auroc,aupr)
#fo = file.path(dirw, '01.tf.rds')
#saveRDS(ev_tf, file = fo)

#ev_go = th %>% mutate(fi = sprintf("%s/13_eval/%s_go.rds", dird, nid)) %>%
    #mutate(res = map(fi, readRDS))
#fo = file.path(dirw, '01.go.rds')
#saveRDS(ev_go, file = fo)

ev_br = th %>% mutate(fi = sprintf("%s/13_eval/%s_br.rds", dird, nid)) %>%
    mutate(res = map(fi, readRDS)) %>%
    select(nid,res) %>% unnest()
fo = file.path(dirw, '01.br.rds')
saveRDS(ev_br, file = fo)

ev_bm = th %>% mutate(fi = sprintf("%s/13_eval/%s_bm.rds", dird, nid)) %>%
    mutate(res = map(fi, readRDS)) %>%
    select(nid,res) %>% unnest()
fo = file.path(dirw, '01.bm.rds')
saveRDS(ev_bm, file = fo)
