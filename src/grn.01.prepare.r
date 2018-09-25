source("grn.fun.r")
dirg = '~/data/genome/B73'
fg = file.path(dirg, "51.gtb")
tg = read_tsv(fg, col_types = "ccciiccccccccccccc") %>%
    transmute(tid = id, gid = par, chrom = chr)
fm = file.path(dirg, "gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read_tsv(fm, col_names = F) %>%
    transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5) %>%
    select(ogid, gid, type)

#{{{ map known TF targets
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
tid = ti %>% filter(fdr_ear<.01, fdr_tassel<.01, fdr_leaf<.01) %>% mutate(tag = 'KN1_any')
tie = ti %>% filter(fdr_ear<.01 | fdr_tassel<.01 | fdr_leaf<.01) %>% mutate(tag = 'KN1_all')
ti2 = rbind(tia, tib, tic, tid, tie)
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

to = ti4 %>% 
    distinct(gid) %>%
    transmute(reg.gid = tfid, tgt.gid = gid)
tag = 'KN1'
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
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
ti = read_tsv(fi, col_names = F)
ti2 = ti %>% transmute(ogid = X2, de1 = X5 ,de2 = X7) %>% 
    filter(de1 == 'yes' | de2 == 'yes')
nrow(ti2)
ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1') 
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d020430'
to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}

#{{{ HDA101
tag = 'HDA101'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi, col_names = F)
ti2 = ti %>% transmute(ogid = X1)
nrow(ti2)
ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 %>% filter(type == '1-to-1') 
sum(ti4$gid %in% tg$gid)

tfid = 'Zm00001d053595'
to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}
#}}}

#{{{ Go/CornCyc/#PPIM
dirg = '~/data/genome/B73/61_functional'
ctag = "GO"
fi = file.path(dirg, "02.go.gs.tsv")
ti = read_tsv(fi)
tig = ti %>% transmute(ctag = ctag, grp = goid, gid = gid, 
                       note = sprintf("%s_%s", evidence, gotype))

ctag = "CornCyc"
fi = file.path(dirg, "07.corncyc.tsv")
ti = read_tsv(fi)
tic = ti %>% transmute(ctag = ctag, grp = pid, gid = gid, note = pname)

t_grp = rbind(tig, tic)
grps = t_grp %>% count(grp) %>% filter(n > 4) %>% pull(grp)
t_grp_f = t_grp %>% filter(grp %in% grps)
t_grp
t_grp_f
t_grp_f %>% count(grp) %>% count(n) %>% print(n=20)

## create all possible regulation pairs
#to = tig %>% group_by(goid) %>%
#    do(data.frame(t(combn(.$gid, 2)))) %>%
#    ungroup() %>%
#    as_tibble() %>%
#    transmute(reg.gid = X1, tgt.gid = X2)
#to2 = to %>% rename(rid = reg.gid, tid = tgt.gid) %>%
#    transmute(reg.gid = tid, tgt.gid = rid)

#{{{ #PPIM
ctag = "PPIM"
dirg = '~/data/genome/B73/61_functional'
fi = file.path(dirg, "08.ppim.tsv")
ti = read_tsv(fi) %>% filter(type1 == '1-to-1', type2 == '1-to-1') 

to1 = ti %>% transmute(reg.gid = gid1, tgt.gid = gid2)
to2 = ti %>% transmute(reg.gid = gid2, tgt.gid = gid1)
to = to1 %>% bind_rows(to2)
fo = sprintf("%s/06_known_ggi/%s.tsv", dird, ctag)
#write_tsv(to, fo)
#}}}
#}}}

#{{{ build GRN gold-standard dataset
ctags = c('KN1_ear', 'KN1_tassel', 'KN1_leaf', 'KN1_any', 'KN1_all',
          'FEA4', 'O2', 'RA1', 'HDA101')
to = tibble()
for (ctag in ctags) {
    fi = sprintf("%s/07_known_tf/05.%s.tsv", dird, ctag)
    ti = read_tsv(fi) %>% 
        mutate(ctag = ctag) %>%
        select(ctag, everything())
    to = rbind(to, ti)
}
t_gs = to

#{{{ TF IDs
ff = '~/data/genome/Zmays_v4/61_functional/06.tf.tsv'
tf = read_tsv(ff)
tf
tf_ids = tf$gid
length(tf_ids)

tf_ids_n = t_gs %>% distinct(reg.gid) %>% filter(!reg.gid %in% tf_ids) %>% 
    pull(reg.gid)
tf_ids_n
tf_ids = c(tf_ids, tf_ids_n)
length(tf_ids)
#}}}

fo = file.path(dird, '09.gs.rda')
save(t_gs, t_grp, t_grp_f, tf_ids, file = fo)
#}}}

#{{{ Walley2016 and Huang2018 GRNs
studies = c(rep("huang", 4), rep("walley", 3))
tags = c('leaf', 'root', 'sam', 'seed',
         'rna', 'protein', 'all')
for (i in 1:length(studies)) {
    study = studies[i]; tag = tags[i]
    fi = sprintf("%s/05_previous_grns/%s_%s.txt", dird, study, tag)
    if(study == 'huang')
        ti = read_tsv(fi)[,1:3]
    else
        ti = read_tsv(fi, col_names = F)[,1:3]
    colnames(ti) = c('rid', 'tid', 'score')
    #
    tn = ti %>% 
        inner_join(tm, by = c('rid' = 'ogid')) %>%
        rename(reg.gid = gid, rtype = type) %>%
        inner_join(tm, by = c('tid' = 'ogid')) %>%
        rename(tgt.gid = gid, ttype = type) %>%
        filter(rtype == '1-to-1', ttype == '1-to-1') %>%
        select(reg.gid, tgt.gid, score)
    rids = unique(tn$reg.gid)
    tids = unique(tn$tgt.gid)
    #   
    fo = sprintf("%s/12_output/np%d.rda", dird, i)
    save(tn, rids, tids, file = fo)
    cat(i, "\n")
}
#}}}

#{{{ #Top45 TFs and targets by Y1H
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


