source("functions.R")
dirw = file.path(dird, '71_share')
tm = v3_to_v4()


fi = file.path(dirw, '11.inari.tsv')
ti = read_tsv(fi, col_names = c('gid.v3')) %>%
    left_join(tm, by=c('gid.v3'='ogid'))
ti %>% print(n=40)
tl = ti %>% mutate(gid = ifelse(is.na(gid), gid.v3, gid)) %>% select(gid)

fn = file.path(dird, 'raw/01.tf.rds')
res = readRDS(fn)

ncfg = th %>% select(nid, lgd)
tn = res %>% select(nid, tn) %>% unnest()

to = tn %>% inner_join(tl, by=c('reg.gid'='gid')) %>%
    select(nid, reg.gid,tgt.gid, pcc) %>%
    inner_join(ncfg, by='nid') %>%
    select(-nid) %>%
    mutate(lgd = factor(lgd, levels=th$lgd)) %>%
    spread(lgd, pcc) %>%
    mutate(n_network = rowSums(!is.na(.[3:49]))) %>%
    mutate(tgt.is.TF = tgt.gid %in% tl$gid) %>%
    select(reg.gid, tgt.gid, tgt.is.TF, n_network, everything()) %>%
    arrange(desc(n_network))

fo = file.path(dirw, '15.supported.supported.tsv')
write_tsv(to, fo)
