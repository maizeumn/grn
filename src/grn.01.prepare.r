source("grn.fun.r")

#{{{ read TFs and v3tov4 mapping
fg = file.path(dirg, "../51.gtb")
tg = read_tsv(fg, col_types = "ccciiccccccccccccc") %>%
    transmute(tid = id, gid = par, chrom = chr)
#
fm = file.path(dirg, "../gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tmr = read_tsv(fm, col_names = F)
colnames(tmr) = c('ogid', 'gid', 'change', 'method', 'type')
tmr = tmr %>% select(-method)
tm = tmr %>% filter(type == '1-to-1')
#
ff = file.path(dirg, "10.tsv")
tf = read_tsv(ff)
#}}}

#{{{ Walley2016 and Huang2018 GRNs
dirw = file.path(dird, '05.previous.grns')
tr = tibble()
ctags = c(rep("huang", 4), rep("walley", 3))
tags = c('leaf', 'root', 'sam', 'seed',
         'rna', 'protein', 'all')
for (i in 1:length(ctags)) {
    ctag = ctags[i]; tag = tags[i]
    fi = sprintf("%s/%s_%s.txt", dirw, ctag, tag)
    if(ctag == 'huang')
        ti = read_tsv(fi)[,1:3]
    else
        ti = read_tsv(fi, col_names = F)[,1:3]
    colnames(ti) = c('regulator', 'target', 'score')
    ti = ti %>%
        mutate(ctag = ctag, tag = tag) %>%
        select(ctag, tag, everything()) %>%
        #top_n(200000, score)
        top_n(1000000, score)
    tr = rbind(tr, ti)
}
tr = tr %>% mutate(ctag = sprintf("%s_%s", ctag, tag)) %>% select(-tag)

tr1 = tr %>% 
    filter(regulator %in% tm$ogid, target %in% tm$ogid) %>%
    inner_join(tm[,c('ogid','gid')], by = c('regulator' = 'ogid')) %>%
    inner_join(tm[,c('ogid','gid')], by = c('target' = 'ogid')) %>%
    mutate(regulator = gid.x, target = gid.y) %>%
    select(-gid.x, -gid.y) %>%
    inner_join(tf, by = c("regulator"="gid"))
tr1 %>% count(ctag)
tr1 %>% distinct(ctag, regulator) %>% count(ctag)
tr1 %>% distinct(ctag, target) %>% count(ctag)

fo = file.path(dirw, '10.RData')
save(t_grn, file = fo)
#}}}

#{{{ kn1 and o2 TF targets
dirw = file.path(dird, '07.known.tf')
ti = tibble()
ctags = c("o2", 'kn1')
for (ctag in ctags) {
    fi = sprintf("%s/%s.tsv", dirw, ctag)
    ti0 = read_tsv(fi) %>% mutate(ctag = ctag) %>%
        distinct(ctag, regulator, target) %>%
        select(ctag, everything())
    ti = rbind(ti, ti0)
}

tr = ti %>% filter(regulator %in% tm$ogid, target %in% tm$ogid) %>%
    inner_join(tm[,c('ogid','gid')], by = c('regulator' = 'ogid')) %>%
    inner_join(tm[,c('ogid','gid')], by = c('target' = 'ogid')) %>%
    mutate(regulator = gid.x, target = gid.y) %>%
    select(-gid.x, -gid.y) %>%
    inner_join(tf, by = c("regulator"="gid"))
tr %>% count(ctag)
tr %>% distinct(ctag, regulator) %>% count(ctag)

t_tf = tr
fo = file.path(dirw, '10.RData')
save(t_tf, file = fo)
#}}}

#{{{ Top45 TFs and targets by Y1H
dirw = file.path(dird, '08.y1h.45')
fi = file.path(dirw, "y1h.targets.tsv")
ti = read_tsv(fi)
colnames(ti) = c("reg.v3", "reg.v4", "tgt.v3", "tgt.v4")
ti = ti %>% fill(reg.v3, reg.v4, .direction = 'down')
ti %>% distinct(reg.v3)
ti %>% distinct(tgt.v3)

tch = ti %>% distinct(reg.v3, reg.v4) %>% 
    left_join(tmr, by = c('reg.v3' = 'ogid')) %>%
    print(n = 45)
tch %>% filter(is.na(gid) | reg.v4 != gid)

tr = ti %>% filter(reg.v4 != 'none', !is.na(tgt.v4)) %>%
    transmute(reg = reg.v4, tgt = tgt.v4)
fr = file.path(dirw, '10.tsv')
write_tsv(tr, fr)
#}}}


