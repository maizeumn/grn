source("functions.R")
dirw = file.path(dird, '02_tfbs_pwm')
diri = file.path(dirw, 'raw')
tm = v3_to_v4()

#{{{ compile motifs
#{{{ process cis-BP motifs
diri1 = file.path(diri, 'cis_bp/TF_Information.txt')
ti = read_tsv(diri1) %>% select(Motif_ID,DBID,TF_Name,
    TF_Status,DBDs,  DBID_1,Motif_Type, MSource_Type,SR_Model) %>%
    filter(Motif_ID != '.') %>%
    mutate(fm = sprintf("%s/cis_bp/pwms_all_motifs/%s.txt", diri, Motif_ID)) %>%
    mutate(size=map_dbl(fm, file.size)) %>%
    select(-fm) %>% arrange(size) %>% filter(size > 12)
ti %>% print(width=Inf)
ti2 = ti %>% select(-size) %>% group_by(DBID) %>% slice(1) %>% ungroup()

ti3 = ti2 %>% inner_join(tm, by=c("DBID"='ogid')) %>% print(width=Inf)
ti3 %>% distinct(DBID, type) %>% count(type)
ti3 %>% distinct(gid, type) %>% count(type)
ti3 %>% count(type)
ti4 = ti3 %>% filter(type=='1-to-1') %>% select(-type)

motifs = ti %>% distinct(Motif_ID) %>% pull(Motif_ID)
fo = file.path(dirw, '02.cisbp.motif_ids.txt')
write(motifs, file=fo)

to1 = ti4 %>% mutate(ctag = 'cisbp') %>%
    select(ctag, gid, motif=Motif_ID, ogid=TF_Name, tf_status=TF_Status,
           method=Motif_Type, src=MSource_Type, src_id=DBID_1)
#}}}

#{{{ plantRegMap
fi = file.path(diri, 'Zma_TF_binding_motifs_information.txt')
ti = read_tsv(fi) %>% select(-Species, -Datasource) %>%
    inner_join(tm, by=c("Gene_id"='ogid')) %>% print(width=Inf)
ti %>% distinct(DBID, type) %>% count(type)
ti %>% distinct(gid, type) %>% count(type)
ti %>% count(type)
ti2 = ti %>% filter(type=='1-to-1') %>% select(-type)

motifs = ti %>% distinct(Matrix_id) %>% pull(Matrix_id)
fo = file.path(dirw, '12.plantregmap.motif_ids.txt')
write(motifs, file=fo)

to2 = ti2 %>% mutate(ctag = 'plantregmap') %>%
    select(ctag, gid, motif=Matrix_id, ogid=Gene_id, method=Method, src_id=Datasource_ID) %>%
    mutate(src_id = str_replace(src_id, 'transfer from ', '')) %>%
    mutate(src_id = str_replace(src_id, '\\(Arabidopsis thaliana\\)', ''))
#}}}

to = to1 %>% bind_rows(to2)
to %>% count(ctag, method, src)

fo = file.path(dirw, "50.tf.motif.tsv")
write_tsv(to, fo)
#}}}

#{{{ process fimo results and filter - need large memory (qsub2l)
fi = file.path(dirw, "50.tf.motif.tsv")
ti = read_tsv(fi) %>% select(ctag, tf=gid, mid=motif)
fi1 = file.path(dirw, '21_fimo_cisbp', 'fimo.tsv')
ti1 = read_tsv(fi1)
fi2 = file.path(dirw, '21_fimo_plantregmap', 'fimo.tsv')
ti2 = read_tsv(fi2)

ti12 = ti1 %>%
    select(mid=1,maid=2,gid=3,start=4,end=5,srd=6,score=7,pval=8,qval=9) %>%
    filter(!str_starts(mid, "#")) %>% select(-maid)
ti22 = ti2 %>%
    select(mid=1,maid=2,gid=3,start=4,end=5,srd=6,score=7,pval=8,qval=9) %>%
    filter(!str_starts(mid, "#")) %>%
    rename(mid = maid, maid = mid)
tx = ti22 %>% distinct(mid, maid) %>% group_by(mid) %>% slice(1)
ti22 = ti22 %>% inner_join(tx, by=c('mid','maid')) %>% select(-maid)

ti13 = ti12 %>% inner_join(ti, by='mid')
ti23 = ti22 %>% inner_join(ti, by='mid')

tb = ti13 %>% bind_rows(ti23) %>% select(ctag, tf, mid, everything())
tb %>% count(ctag, tf) %>% summary()

tbf = tb %>% select(-mid) %>% filter(pval < 1e-6)
tbf %>% count(ctag, tf) %>% summary()

# convert to genomic coordinate
tbs = tbf %>% distinct(gid, start, end, srd) %>% arrange(gid, start) %>% rename(chrom=gid)

fc = '~/projects/genome/data/Zmays_B73/50_annotation/16.promoter.2kb.chain'
tbs2 = liftover(tbs, fc)
tbs3 = tbs2 %>% rename(gid=chrom,gstart=start,gend=end,chrom=chrom2,start=start2,end=end2)

to = tbf %>% rename(gstart=start, gend=end) %>%
    inner_join(tbs3, by=c('gid','gstart','gend','srd')) %>%
    filter(!is.na(chrom)) %>% mutate(dist_tss = gend - 2000) %>%
    select(ctag, tf, chrom, start, end, srd, score, gid, dist_tss) %>%
    arrange(ctag, tf, chrom, start)
to %>% count(ctag, tf) %>% summary()
summary(to$dist_tss)

fo = file.path(dirw, '55.tfbs.fimo.tsv')
write_tsv(to, fo)
#}}}







