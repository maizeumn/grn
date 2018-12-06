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


