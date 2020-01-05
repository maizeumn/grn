source("functions.R")
dirw = file.path(dird, '91_recovery')

#{{{ read data
fi = file.path(dirw, '00.exp.mat.csv')
em = read_csv(fi) %>% rename(gid = 1) %>%
    gather(cond, cpm, -gid)

ft = file.path(dirw, '00.tf.tsv')
str_c_rm_na <- function(x, sep=',', collapse=',') {
    #{{{
    x2 = x[!is.na(x)]
    ifelse(length(x2)==0, '', str_c(x2, sep=sep, collapse=collapse))
    #}}}
}
tf = read_tsv(ft, col_names=F, col_types='ccccccccc') %>%
    select(fam=1, gid=2, gname=3, note=4) %>%
    mutate(gid = str_to_upper(gid)) %>%
    group_by(fam, gid) %>%
    summarise(gname = str_c_rm_na(gname), note=str_c_rm_na(note)) %>%
    ungroup()

fd = file.path(dirw, '00.decay.csv')
decay = read_csv(fd) %>% select(gid=1, tid=2, note=3, type=4, gname=5)
fr = file.path(dirw, '00.rbp.csv')
rbp = read_csv(fr) %>% select(gid=1, gname=2, type=3) %>%
    mutate(gid=str_to_upper(str_replace(gid, "\xca", '')))

nrow(tf)
sum(tf$gid %in% em$gid)
nrow(decay)
sum(decay$gid %in% em$gid)
nrow(rbp)
sum(rbp$gid %in% em$gid)

tg = em %>% distinct(gid) %>% mutate(type = 'gene') %>%
    mutate(type = ifelse(gid %in% tf$gid, str_c(type,'TF',sep=','), type)) %>%
    mutate(type = ifelse(gid %in% rbp$gid, str_c(type,'RBP',sep=','), type)) %>%
    mutate(type = ifelse(gid %in% decay$gid, str_c(type,'decay',sep=','), type))
tg %>% count(type)

res = list(em=em, tf=tf, rbp=rbp, decay=decay, tg=tg)
fo = file.path(dirw, '05.raw.rds')
saveRDS(res, file=fo)
#}}}

#{{{ filtering
fi = file.path(dirw, '05.raw.rds')
res = readRDS(fi)

min_cpm = 1
num_sam_on = 0
pct_sam_on = .05
min_var_p = 0
ems = res$em %>% rename(val=cpm) %>% group_by(gid) %>%
    summarise(nsam_on = sum(val >= min_cpm),
              psam_on = nsam_on/n(),
              val_sd = sd(val)) %>%
    ungroup()
gids = ems %>%
    filter(nsam_on >= num_sam_on,
           psam_on >= pct_sam_on,
           val_sd >= 0) %>% pull(gid)

ems2 = ems %>% filter(gid %in% gids)
min_sd = quantile(ems2$val_sd, min_var_p)
gids = ems2 %>% filter(val_sd >= as.numeric(min_sd)) %>% pull(gid)

tg = res$tg %>% filter(gid %in% gids)
tg %>% count(type)
em = res$em %>% filter(gid %in% gids) %>%
    mutate(cpm = asinh(cpm)) %>%
    spread(cond, cpm)

fo1 = file.path(dirw, '10.exp.mat.tsv')
write_tsv(em, fo1)
fo2 = file.path(dirw, '10.reg.txt')
regs = tg %>% filter(type != 'gene') %>% pull(gid)
write(regs, fo2)
#}}}

#{{{ get GRN output
fi = file.path(dirw, '05.raw.rds')
res = readRDS(fi)
fi = file.path(dirw, '13.grn.rds')
x = readRDS(fi)

tn = x$tn %>% select(reg.gid, tgt.gid, score, pcc, pval.p, spc, pval.s) %>%
    inner_join(res$tg, by=c('reg.gid'='gid')) %>% rename(reg.type=type) %>%
    inner_join(res$tg, by=c('tgt.gid'='gid')) %>% rename(tgt.type=type) %>%
    slice(1:1e6) %>% mutate(i = 1:1e6) %>%
    select(i, reg.gid, tgt.gid, score, reg.type, tgt.type, everything())
tn %>% slice(1:1e4) %>% count(reg.type)

fo = file.path(dirw, '15.grn.tsv')
write_tsv(tn, fo)
#}}}



