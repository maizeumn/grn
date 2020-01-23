source("functions.R")
dirw = file.path(dird, '04_tfbs')
diro = file.path(dird, '04_tfbs')

#{{{ read in
fi1 = file.path(dirw, '../02_tfbs_pwm/55.tfbs.fimo.tsv')
ti1 = read_tsv(fi1)
fi2 = file.path(dirw, '../03_tfbs_chipseq/Ricci2019.bed')
ti2 = read_tsv(fi2, col_names=c('chrom','start','end','srd','tf','tgt.gid','dist_tss')) %>%
    mutate(ctag='Ricci2019', start=start+1, tf=str_replace(tf, '-[0-9]+$', ''))

tb = ti1 %>% select(ctag,chrom,start,end,srd,tf, tgt.gid=gid,dist_tss) %>%
    bind_rows(ti2)
tb %>% distinct(ctag, tf) %>% count(ctag)
tb %>% count(ctag, tf) %>% summary()
#}}}

#{{{ assess UMR and ACR
acrE = read_regions('acrE')
acrL = read_regions('acrL')
umr = read_regions('umr')
tbs = tb %>% distinct(chrom, start, end, srd) %>% arrange(chrom, start)

o_acrE = intersect_once(tbs, acrE, bp_pct=.5)
o_acrL = intersect_once(tbs, acrL, bp_pct=.5)
o_umr = intersect_once(tbs, umr, bp_pct=.5)
tbs2 = tbs %>% mutate(acrE=o_acrE, acrL=o_acrL, umr=o_umr)
tbs2 %>% count(acrE, acrL, umr)

to = tb %>% inner_join(tbs2, by=c('chrom','start','end','srd'))

fo = file.path(dirw, '10.tfbs.acr.umr.tsv')
write_tsv(to, fo)
#fo = file.path(dirw, '10.tfbs.acr.umr.rds')
#saveRDS(to, fo)
#}}}

#{{{ generate TFBS and TF-tgt regulation
fi = file.path(dirw, '10.tfbs.acr.umr.tsv')
ti = read_tsv(fi) %>% filter(dist_tss>=-2000, dist_tss<=200)

t1 = ti %>% distinct(ctag, tf, tgt.gid)
t1 %>% count(ctag, tf) %>% summary()
t1 %>% count(ctag, tf) %>% arrange(-n)

t2 = ti %>% filter(umr) %>% distinct(ctag, tf, tgt.gid) %>%
    mutate(ctag = str_c(ctag, 'umr', sep="_"))
t2 %>% count(ctag, tf) %>% summary()

t3a = ti %>% filter(acrE) %>% distinct(ctag, tf, tgt.gid) %>%
    mutate(ctag = str_c(ctag, 'acrE', sep="_"))
t3a %>% count(ctag, tf) %>% summary()

t3b = ti %>% filter(acrL) %>% distinct(ctag, tf, tgt.gid) %>%
    mutate(ctag = str_c(ctag, 'acrL', sep="_"))
t3b %>% count(ctag, tf) %>% summary()

#{{{ [side project] Yi-Hsuan's P1 targets
d04 = file.path(dird, '04_tfbs')
fi = file.path(d04, "P1/00.xlsx")
ti = read_xlsx(fi) %>%
    select(chip_v2=1, chip_v4=2, dap_M=3, dap_deM=4, rna_mutant=5, rna_proto=6,
        chip_rna=7, dap_rna_proto=8, true_tgts=9) %>%
    gather(ctag, gid) %>%
    mutate(ctag = str_c("P1", ctag, sep="|")) %>%
    filter(!is.na(gid)) %>%
    mutate(tf = 'Zm00001d028842') %>% select(ctag, tf, tgt.gid=gid)
t9 = ti
fo = file.path(dirw, "P1/10.tsv")
#write_tsv(t9, fo)
#}}}

to = rbind(t1, t2, t3a, t3b, t9)
tof = to %>% count(ctag, tf) %>% filter(ctag=='P1|true_tgts' | n >= 30) %>% select(ctag, tf)
to = to %>% inner_join(tof, by=c('ctag','tf'))
to %>% count(ctag, tf) %>% arrange(n)
to %>% count(ctag, tf) %>% filter(str_starts(ctag, 'Ricci')) %>% mutate(x=map_int(tf, nchar)) %>% print(n=50)
to %>% count(ctag, tf) %>% filter(str_starts(ctag, 'P1'))

fo = file.path(diro, '15.regulations.tsv')
write_tsv(to, fo)
#}}}

