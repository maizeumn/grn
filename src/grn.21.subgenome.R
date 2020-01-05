source("functions.R")
diri = '~/projects/rnaseq'
dirw = file.path(dird, '21_subgenome')
#
ft = '~/projects/grn/data/09.tf.txt'
tf_ids = read_tsv(ft, col_names = F) %>% pull(X1)
t_tf0 = tibble(gid=tf_ids, tf = 'TF')
#
tsyn = read_syn(gcfg)

#{{{ TF - syntenic composition
tp0 = tsyn %>% left_join(t_tf0, by = 'gid') %>%
    replace_na(list(tf='non-TF'))
tps = tp0 %>% count(tf) %>% mutate(lab = sprintf("%s (N=%s)", tf, number(n)))
tp = tp0 %>%
    group_by(tf, ftype) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(lab = percent(freq, accuracy=2)) %>%
    arrange(tf, desc(ftype)) %>%
    group_by(tf) %>% mutate(y = cumsum(freq)) %>%
    mutate(y = y - freq/2) %>% ungroup()

cols5 = c('grey70', brewer.pal(6, 'Paired')[3:6])
p1 = ggplot(tp) +
    geom_bar(aes(x=tf, y=n, fill=ftype), stat='identity', position='fill', width=.8) +
    geom_text(aes(x=tf, y=y, label=lab), color='black') +
    scale_x_discrete(expand=c(0,0), breaks=tps$tf, labels=tps$lab) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(name='Syntenic status', values=cols5) +
    otheme(legend.pos = 'top', legend.dir='v', legend.title = F,
        xtick=T, ytick=F,
        xtext=T, ytext=F)
fo = file.path(diro, '01.tf.syn.pdf')
ggsave(p1, file=fo, width=4, height=6)
#}}}

#{{{ check correlation of TF paris
yid = 'rnc01'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

tpm = tm_m %>% select(gid, SampleID, CPM) %>%
    arrange(gid, SampleID) %>% group_by(gid) %>%
    summarise(cpm = list(CPM)) %>% ungroup()

tp = tsyn %>% filter(!is.na(maize1), !is.na(maize2), maize1 %in% tf_ids, maize2 %in% tf_ids) %>%
    inner_join(tpm, by=c('maize1'='gid')) %>% rename(cpm1=cpm) %>%
    inner_join(tpm, by=c('maize2'='gid')) %>% rename(cpm2=cpm) %>%
    mutate(res = map2(cpm1, cpm2, cor.test)) %>%
    mutate(pcc = map_dbl(res, 'estimate'), pval = map_dbl(res, 'p.value')) %>%
    select(maize1, maize2, pcc, pval)
summary(-log10(tp$pval))

tp = tsyn %>% filter(!is.na(maize1), !is.na(maize2), !maize1 %in% tf_ids, !maize2 %in% tf_ids) %>%
    inner_join(tpm, by=c('maize1'='gid')) %>% rename(cpm1=cpm) %>%
    inner_join(tpm, by=c('maize2'='gid')) %>% rename(cpm2=cpm) %>%
    mutate(res = map2(cpm1, cpm2, cor.test)) %>%
    mutate(pcc = map_dbl(res, 'estimate'), pval = map_dbl(res, 'p.value')) %>%
    select(maize1, maize2, pcc, pval)
summary(-log10(tp$pval))
#}}}

fv = sprintf("%s/rf.100k.rds", dirr)
ev = readRDS(fv)
tx = ev %>% select(nid, tn) %>% unnest()
tx2 = tx %>% count(nid, reg.gid)

tsyn = read_syn(gcfg, opt=2)

tt = tsyn %>% filter(maize1 %in% tf_ids | maize2 %in% tf_ids) %>%
    select(-sorghum) %>% distinct(maize1, maize2, ftype)

tx3 = tx2 %>% filter(nid == 'nc01')
tp = tt %>% left_join(tx3, by=c('maize1'='reg.gid')) %>% rename(n1=n) %>%
    left_join(tx3, by=c('maize2'='reg.gid')) %>% rename(n2=n) %>%
    replace_na(list(n1=0,n2=0))
tp %>% filter(ftype=='both-retained') %>% filter(n1>0, n2>0) %>%
    count(n1 > n2)

len_ovlp <- function(x1, x2) sum(x1 %in% x2)
tx3 = tx %>% filter(nid == 'nc01') %>% select(reg.gid, tgt.gid) %>%
    group_by(reg.gid) %>% summarise(tgts = list(tgt.gid)) %>% ungroup()
tp = tt %>% filter(ftype=='both-retained') %>%
    left_join(tx3, by=c('maize1'='reg.gid')) %>% rename(tgts1=tgts) %>%
    left_join(tx3, by=c('maize2'='reg.gid')) %>% rename(tgts2=tgts) %>%
    replace_na(list(n1=0,n2=0)) %>%
    mutate(n1 = map_int(tgts1, length)) %>%
    mutate(n2 = map_int(tgts2, length)) %>%
    mutate(n12 = map2_int(tgts1, tgts2, len_ovlp))

tp %>% filter(n1>50 & n2 >50) %>% select(maize1,maize2,n1,n2,n12) %>%
    print(n=400)


