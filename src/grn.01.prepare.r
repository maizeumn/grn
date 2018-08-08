source("grn.fun.r")

#{{{ map known TF targets
dirg = '~/data/genome/B73'
fg = file.path(dirg, "51.gtb")
tg = read_tsv(fg, col_types = "ccciiccccccccccccc") %>%
    transmute(tid = id, gid = par, chrom = chr)
fm = file.path(dirg, "gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read_tsv(fm, col_names = F) %>%
    transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5)
dirw = file.path(dird, '07_known_tf')

#{{{ KN1
tag = 'KN1'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi)
ti = ti %>%
    transmute(ogid = `Gene ID`,
              binding = `High confidence binding loci within 10 kb`,
              fdr_ear = `FDR ear`,
              fdr_tassel = `FDR tassel`,
              fdr_SAM = `FDR sam`,
              fdr_leaf = `FDR leaf homo/wt`,
              fdr_leaf_het = `FDR leaf het/wt`
    ) %>%
    filter(binding == 'yes')

tia = ti %>% filter(fdr_ear < .01) %>% mutate(tag = 'KN1_ear')
tib = ti %>% filter(fdr_tassel < .01) %>% mutate(tag = 'KN1_tassel')
tic = ti %>% filter(fdr_leaf < .01) %>% mutate(tag = 'KN1_leaf')
ti2 = rbind(tia, tib, tic)
ti2 %>% count(tag)

ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1')
ti4 %>% count(tag)
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d033859'
tags = unique(ti4$tag)
for (tag in tags) {
    to = ti4 %>% filter(tag == !!tag) %>%
        distinct(gid) %>%
        transmute(reg.gid = tfid, tgt.gid = gid)
    fo = sprintf("%s/05.%s.tsv", dirw, tag)
    write_tsv(to, fo)
}
#}}}

#{{{ FEA4
tag = 'FEA4'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi, col_names = F)
ti2 = ti %>%
    transmute(ogid = X1, note = X2, DE = X3) %>%
    filter(DE == 'yes')
ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1') 
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d037317'
to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}

#{{{ O2 
# liftOver 01.v3.coord.bed $genome/B73/chain/AGP_v3_to_v4.bed.gz 02.v4.coord.bed 03.unmapped.bed
# slopBed -i 02.v4.coord.bed -g $genome/B73/15.sizes -b 10000 > 05.flank10k.bed
# intersectBed -a $genome/B73/v37/gene.bed -b 05.flank10k.bed -u > 06.ovlp.gid.bed
#fi = file.path(dirw, "o2/06.ovlp.gid.bed")
#ti = read_tsv(fi, col_names = F, col_types = c('ciic'))
tag = 'O2'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi, col_names = F) %>%
    transmute(ogid = X1)
ti3 = ti %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1') 
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d018971'
to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}

#{{{ RA1
tag = 'RA1'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi)
ti2 = ti %>% filter(!is.na(DEb)) %>% transmute(ogid = maize.gene.id)
nrow(ti2)
ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1') 
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d037317'
to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}

ctags = c('KN1_ear', 'KN1_tassel', 'KN1_leaf',
          'FEA4', 'O2', 'RA1')
to = tibble()
for (ctag in ctags) {
    fi = sprintf("%s/05.%s.tsv", dirw, ctag)
    ti = read_tsv(fi) %>% 
        mutate(ctag = ctag) %>%
        select(ctag, everything())
    to = rbind(to, ti)
}
t_gs = to
t_gs %>% count(ctag)
fo = file.path(dirw, '10.rda')
save(t_gs, file = fo)
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


