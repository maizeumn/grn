source("functions.R")
genome = 'B73'
genome_cfg = read_genome_conf(genome)
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

ctags = c('li2013','liu2017','wang2018')
hs = tibble(ctag=ctags) %>%
    mutate(fi=file.path('~/projects/genomes/data',ctag,'10.rds')) %>%
    mutate(data = map(fi, readRDS)) %>%
    mutate(hs = map(data, 'hs.tgt')) %>%
    select(ctag, hs) %>% unnest() %>%
    select(ctag,grp=qid,gid) %>% mutate(note=NA)

fun_ann = rbind(go_hc, go, cc, hs, y1h)
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
    filter(ctag != 'HDA101')
ctags_tf5 = tf %>% count(ctag) %>% arrange(n) %>% pull(ctag)
tf2 = tf %>% mutate(ctag = 'known_TFs') %>% distinct(ctag, reg.gid, tgt.gid)
tf = tf %>% bind_rows(tf2)
tf %>% count(ctag)
#}}}

#{{{ FunTFBS regulations
fi = file.path(dird, '03_tfbs', 'regulation_merged.txt')
ti = read_tsv(fi, col_names=c("o.reg.gid",'relation','o.tgt.gid','org','evi'))
ti %>% count(relation, org, evi)
ti %>% count(o.reg.gid)

to = ti %>% select(o.reg.gid, o.tgt.gid, evi) %>%
    inner_join(tm, by=c("o.reg.gid"="ogid")) %>%
    rename(reg.gid=gid, reg.type=type) %>%
    inner_join(tm, by=c("o.tgt.gid"="ogid")) %>%
    rename(tgt.gid=gid, tgt.type=type) %>%
    filter(reg.gid %in% genome_cfg$size.gene$gid, tgt.gid %in% genome_cfg$size.gene$gid) %>%
    distinct(reg.gid, tgt.gid, reg.type, tgt.type, evi)
to %>% count(reg.type,tgt.type)
to %>% count(reg.gid) %>% arrange(n)

tfbs = to %>% distinct(reg.gid, tgt.gid, evi) %>%
    transmute(ctag=evi, reg.gid=reg.gid, tgt.gid=tgt.gid)
ctags_tfbs = tfbs %>% count(ctag) %>% arrange(n) %>% pull(ctag)
tfbs %>% count(ctag)
tnk = tf %>% bind_rows(tfbs)
tnk %>% count(ctag)
#}}}

#{{{ FunTFBS (no use)
fi = file.path(dird, '03_tfbs', 'TFBS_from_FunTFBS_inProm.gff')
ti = read_tsv(fi, col_names=c('chr','src','type','beg','end','score','srd','phase','note'), col_types='cccdddccc')
ti %>% count(src,type)

to = ti %>% select(chr,beg,end,score,note) %>%
    separate(note, c('gid','tid','corr','pval','seq'), sep=';', extra='drop') %>%
    separate(gid, c('tag','gid'), sep='=', extra='drop') %>%
    separate(tid, c('tag','tid'), sep='=', extra='drop') %>%
    separate(corr, c('tag','corr'), sep='=', extra='drop') %>%
    separate(pval, c('tag','pval'), sep='=', extra='drop') %>%
    mutate(corr=as.numeric(corr), pval=as.numeric(pval)) %>%
    select(chr,beg,end,score,tid,gid,corr,pval)
#}}}

#{{{ TF IDs
ff = '~/data/genome/Zmays_v4/61_functional/06.tf.tsv'
ti = read_tsv(ff)
all_tf = ti
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

# optional: run grn.91.tf45.R

# build GRN gold-standard dataset
res = list(tnk=tnk, ctags_tf5=ctags_tf5, ctags_tfbs=ctags_tfbs,
           all_tf=all_tf, tf_ids=tf_ids, fun_ann=fun_ann, ppi=ppi)
fo = file.path(dird, '09.gs.rds')
saveRDS(res, file=fo)
ft = file.path(dird, '09.tf.txt')
write(tf_ids, file=ft)

