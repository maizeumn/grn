source("functions.R")
genome = 'B73'
x = load(file.path(dirg, genome, '55.rda'))
tm = v3_to_v4()

#{{{ map known TF targets
dirw = file.path(dird, '07_known_tf')
#{{{ KN1
tag = 'KN1'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi) %>%
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

#{{{ bZIP22
tag = 'bZIP22'
fi = sprintf("%s/raw/%s.tsv", dirw, tag)
ti = read_tsv(fi, col_names = F)
tfid = ti$X1[1]
ti2 = ti %>% transmute(ogid = X2)
nrow(ti2)
ti3 = ti2 %>% inner_join(tm, by = 'ogid')
ti3 %>% count(type)
ti4 = ti3 #%>% filter(type == '1-to-1')
sum(ti4$gid %in% t_gs$gid)

to = tibble(reg.gid = tfid, tgt.gid = unique(ti4$gid))
fo = sprintf("%s/05.%s.tsv", dirw, tag)
write_tsv(to, fo)
#}}}
#}}}

#{{{ Top45 TFs and targets by Y1H
dirw = file.path(dird, '08_y1h')
fi = file.path(dirw, "y1h.targets.tsv")
ti = read_tsv(fi)
colnames(ti) = c("reg.v3", "reg.v4", "tgt.v3", "tgt.v4")
ti = ti %>% fill(reg.v3, reg.v4, .direction = 'down')
ti %>% distinct(reg.v3)
ti %>% distinct(tgt.v3)

tch = ti %>% distinct(reg.v3, reg.v4) %>%
    left_join(tm, by = c('reg.v3' = 'ogid')) %>%
    print(n = 45)
tch %>% filter(is.na(gid) | reg.v4 != gid)

tn = ti %>% filter(reg.v4 != 'none', !is.na(tgt.v4)) %>%
    transmute(reg = reg.v4, tgt = tgt.v4)
fn = file.path(dirw, '10.tsv')
#write_tsv(tn, fn)


fi = file.path(dirw, 'phenolic.xlsx')
ti = read_xlsx(fi,
    col_names = c('gid_v3','gname','yeast_prom','prom_mplant',
                 'gid','chrom','start','end','orig','note',
                 'ref','x1','coexp','syntelog','sub1','sub2_ref')) %>%
    filter(row_number() > 1) %>%
    filter(!is.na(gid) & gid != 'Na') %>%
    select(gid,gname,everything()) %>%
    filter(gid_v3 != 'GRMZM2G049424')
ti %>% count(gid) %>% count(n)

fo = file.path(dirw, '01.rds')
res = list(reg.gids = unique(tn$reg), tgt.gids = ti$gid, tn = tn)
saveRDS(res, file=fo)
#}}}

#{{{ Walley2016 and Huang2018 GRNs
lift_previous_grn <- function(nid, study, tag, tm, dird = '~/projects/maize.grn/data') {
    #{{{
    fi = sprintf("%s/05_previous_grns/%s_%s.txt", dird, study, tag)
    if(study == 'huang')
        ti = read_tsv(fi)[,1:3]
    else
        ti = read_tsv(fi, col_names = F)[,1:3]
    colnames(ti) = c('rid', 'tid', 'score')
    tn = ti %>%
        inner_join(tm, by = c('rid' = 'ogid')) %>%
        rename(reg.gid = gid, rtype = type) %>%
        inner_join(tm, by = c('tid' = 'ogid')) %>%
        rename(tgt.gid = gid, ttype = type) %>%
        filter(rtype == '1-to-1', ttype == '1-to-1') %>%
        select(reg.gid, tgt.gid, score)
    rids = unique(tn$reg.gid)
    tids = unique(tn$tgt.gid)
    reg.mat = tn %>% spread(tgt.gid, score) %>%
        as.data.frame() %>% column_to_rownames(var = 'reg.gid')
    cat(sprintf("%s %s %s: %d edges, %d TFs, %d targets\n", nid, study, tag,
                nrow(tn), length(rids), length(tids)))
    fo = sprintf("%s/12_output/%s.rda", dird, nid)
    save(reg.mat, rids, tids, tn, file = fo)
    TRUE
    #}}}
}

tp = th %>% filter(nid %in% c("np16_1", sprintf("np18_%d", 1:4))) %>%
    transmute(nid=nid, study=str_replace(study,"\\d+$",''), tag=note)

pmap_lgl(tp, lift_previous_grn, tm)
#}}}

#{{{ functional annotation: GO CornCyc Y1H
dirg = '~/data/genome/B73/61_functional'
fi = file.path(dirg, "02.go.gs.tsv")
ti = read_tsv(fi)
go_hc = ti %>% mutate(note = str_c(evidence, gotype)) %>%
    mutate(ctag = 'GO_HC') %>% select(ctag, grp=goid, gid, note)

fi = file.path(dirg, "01.go.tsv")
ti = read_tsv(fi)
go = ti %>%
    filter(!ctag %in% c('aggregate','fanngo'), gotype=='P') %>%
    mutate(ctag = str_c("GO",ctag,sep="_")) %>%
    select(ctag, grp=goid, gid, note=goname)

ctag = "CornCyc"
fi = file.path(dirg, "07.corncyc.tsv")
ti = read_tsv(fi)
cc = ti %>% transmute(ctag = !!ctag, grp = pid, gid = gid, note = pname)

ctag = 'Y1H'
fi = '~/projects/grn/data/08_y1h/01.rds'
y1h = readRDS(fi)
y1h = tibble(ctag = !!ctag, grp = 'Y1H', gid = y1h$tgt.gids, note = '')

ctag = "PPIM"
fi = file.path(dirg, "08.ppim.tsv")
ti = read_tsv(fi)
ppi = ti %>% filter(type1 != '1-to-0', type2 != '1-to-0') %>%
    select(gid1, gid2)
ppic = ppi %>% mutate(grp = sprintf('ppi%d', 1:nrow(ppi))) %>%
    gather(tag, gid, -grp) %>%
    transmute(ctag='PPIM', grp=grp, gid=gid, note='')

fun_ann = rbind(go_hc, go, cc, y1h)
fun_ann %>% distinct(ctag,grp) %>% count(ctag)
#}}}

#{{{ known TF/target pairs
ctags = c('KN1_ear', 'KN1_tassel', 'KN1_leaf', 'KN1_any', 'KN1_all',
          'FEA4', 'O2', 'RA1', 'HDA101', 'bZIP22')
tf = tibble()
for (ctag in ctags) {
    fi = sprintf("%s/07_known_tf/05.%s.tsv", dird, ctag)
    ti = read_tsv(fi) %>%
        mutate(ctag = ctag) %>%
        select(ctag, everything())
    tf = rbind(tf, ti)
}
tf = tf %>%
    filter(! ctag %in% c("KN1_any","KN1_ear","KN1_tassel","KN1_leaf")) %>%
    mutate(ctag = ifelse(ctag=='KN1_all', 'KN1', ctag)) %>%
    mutate(binding = 1)
tf %>% count(ctag)
tfs = tf %>% distinct(ctag, reg.gid)
#}}}

#{{{ TF IDs
ff = '~/data/genome/Zmays_v4/61_functional/06.tf.tsv'
ti = read_tsv(ff)
tf_ids = ti$gid
length(tf_ids)
length(unique(tf_ids))
#
tf_ids_n = tf %>% distinct(reg.gid) %>% filter(!reg.gid %in% tf_ids) %>%
    pull(reg.gid)
tf_ids_n
tf_ids = c(tf_ids, tf_ids_n)
length(tf_ids)
#}}}

# build GRN gold-standard dataset
res = list(tf=tf, tfs=tfs, tf_ids=tf_ids, fun_ann=fun_ann, ppi=ppi)
fo = file.path(dird, '09.gs.rds')
saveRDS(res, file=fo)

