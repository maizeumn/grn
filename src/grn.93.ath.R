source("functions.R")
dirw = file.path(dird, '93_ath')

#{{{ read data
fi = file.path(dirw, '00.exp.mat.tsv')
ti = read_tsv(fi)
em = ti %>% select(-Entrez) %>% select(gid=TAIRid, everything()) %>%
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

tg = em %>% distinct(gid) %>% mutate(type = 'gene') %>%
    mutate(type = ifelse(gid %in% tf$gid, str_c(type,'TF',sep=','), type))
tg %>% count(type)

res = list(em=em, tf=tf, tg=tg)
fo = file.path(dirw, '05.raw.rds')
saveRDS(res, file=fo)
#}}}

#{{{ filtering
fi = file.path(dirw, '05.raw.rds')
res = readRDS(fi)

#{{{ do not filter
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
#}}}
gids = unique(res$tg$gid)

tg = res$tg %>% filter(gid %in% gids)
tg %>% count(type)
em = res$em %>% filter(gid %in% gids) %>%
#    mutate(cpm = asinh(cpm)) %>%
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

tn = x$tn %>% mutate(direction=ifelse(pcc < 0, '-', '+')) %>%
    select(reg.gid, tgt.gid, score, direction)

fo = file.path(dirw, '15.grn.tsv')
write_tsv(tn, fo)
#}}}

#{{{ At PDI
#{{{ read
fi = file.path(dirw, '13.grn.rds')
x = readRDS(fi)
net_size = 10e6
tn = x$tn %>% select(reg.gid, tgt.gid, score) %>%
    filter(reg.gid != tgt.gid) %>%
    mutate(score = as.numeric(score)) %>%
    mutate(score = ifelse(is.na(score), 0, score)) %>%
    arrange(desc(score)) %>%
    filter(row_number() <= net_size)
#
fi = file.path(dirw, '00.pdi.tsv')
ti = read_tsv(fi, col_names=c('reg.gid','opt','tgt.gid')) %>% mutate(weight=T)
td = ti %>% filter(opt == 'DAPseq')
tc = ti %>% filter(opt == 'ChIP')
#}}}

#{{{ TF-centric
eval_tf <- function(tx, tids) { # score, weight
#{{{
    to = tx %>% arrange(-score) %>% full_join(tibble(tgt.gid=tids), by='tgt.gid') %>%
        replace_na(list(score = 0, weight = F))
    if(max(to$score) == 0) to$score[1] = 0.1
    resR = roc.curve(scores.class0=to$score, weights.class0=to$weight)
    resP = pr.curve(scores.class0=to$score, weights.class0=to$weight)
    scores1 = to %>% filter(weight) %>% pull(score)
    scores2 = to %>% filter(!weight) %>% pull(score)
    pval = NA
    if(length(scores1) > 0 & length(scores2) > 0)
        pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
    auroc = resR$auc
    auprc = resP$auc.integral
    tibble(auroc=auroc, auprc=auprc, pval=pval)
#}}}
}
fi = file.path(dirw, '05.raw.rds')
res = readRDS(fi)

#{{{ eval chip
ti = tc
rids = unique(ti$reg.gid)
tids = unique(ti$tgt.gid)
rids = rids[rids %in% tn$reg.gid]
length(rids)
#
tn1 = tn %>% filter(reg.gid %in% rids)
ti1 = ti %>% filter(reg.gid %in% rids)
to = tn1 %>%
    full_join(ti1, by = c('reg.gid','tgt.gid')) %>%
    mutate(grn=!is.na(score), pdi=!is.na(opt)) %>%
    replace_na(list(score=0, weight=F))
to %>% count(grn, pdi)
#
tp = to %>% group_by(reg.gid) %>% nest() %>%
    mutate(n_tgt = map_int(data, x <- function(tt) sum(tt$weight))) %>%
    mutate(res = map(data, eval_tf, tids=tids)) %>%
    select(reg.gid, n_tgt, res) %>% unnest()
tp %>% print(n=50)
toc = tp %>% mutate(opt = 'chip-seq')
#}}}
#{{{ eval dap
ti = td
rids = unique(ti$reg.gid)
tids = unique(ti$tgt.gid)
rids = rids[rids %in% tn$reg.gid]
length(rids)
#
tn1 = tn %>% filter(reg.gid %in% rids)
ti1 = ti %>% filter(reg.gid %in% rids)
to = tn1 %>%
    full_join(ti1, by = c('reg.gid','tgt.gid')) %>%
    mutate(grn=!is.na(score), pdi=!is.na(opt)) %>%
    replace_na(list(score=0, weight=F))
to %>% count(grn, pdi)
#
tp = to %>% group_by(reg.gid) %>% nest() %>%
    mutate(n_tgt = map_int(data, x <- function(tt) sum(tt$weight))) %>%
    mutate(res = map(data, eval_tf, tids=tids)) %>%
    select(reg.gid, n_tgt, res) %>% unnest()
tp %>% print(n=50)
tod = tp %>% mutate(opt = 'dap-seq')
#}}}
tp = rbind(toc, tod) %>% select(opt, everything()) %>% rename(rid=reg.gid) %>%
    inner_join(res$tf, by = c('rid'='gid'))

oids = tp %>% filter(pval<1e-5) %>% count(rid) %>% filter(n>1) %>% pull(rid)
tps = tp %>% filter(rid %in% oids)
p = ggplot(tp) +
    geom_point(aes(auroc, auprc, color=-log10(pval)), size=2) +
    geom_text_repel(data=tps, aes(auroc, auprc, label=rid), nudge_x=.01,nudge_y=.05, size=2) +
    scale_x_continuous(name='AUROC', expand= c(.03,.03)) +
    scale_y_continuous(name='AUPRC', expand= c(.03,.03)) +
    scale_color_viridis(name='-log10(P)', option='viridis', direction=-1) +
    facet_wrap(~opt, nrow=1, scale='free') +
    otheme(xtitle=T,ytitle=T, xtext=T,ytext=T, xtick=T, ytick=T,
        legend.pos = 'top.center.out', legend.dir='h', legend.title=T)
fo = file.path(dirw, '21.tf.pdf')
ggsave(fo, p, width=10, height=6)

tt = tp %>% filter(pval < 1e-10) %>%
    arrange(pval) %>%
    print(n=50)
tt %>% filter(rid %in% oids) %>% print(width=Inf)

tp %>% distinct(fam, rid) %>% count(fam) %>% mutate(prop=n/sum(n))%>% arrange(-n) %>% print(n=20)
tt %>% distinct(fam, rid) %>% count(fam) %>% mutate(prop=n/sum(n))%>% arrange(-n)

#{{{ # auc [ obsolete ]
plot_auc <- function(resR, resP, fo) {
#{{{
    tpr = resR$curve %>% as_tibble() %>%
        transmute(TPR = V1, FPR = V2, score = V3)
    tpp = resP$curve %>% as_tibble() %>%
        transmute(recall = V1, precision = V2, score = V3)
    p1 = ggplot(tpr) +
        geom_line(mapping = aes(x = TPR, y = FPR, color = 'auc')) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
        scale_x_continuous(name = 'FPR: FP/(FP+TN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
        scale_y_continuous(name = 'TPR: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
        scale_color_manual(values = pal_aaas()(2)) +
        otheme(strip.size = 9, legend.pos = 'none', margin = c(.5,.5,.1,.5),
               xtitle=T, ytitle=T, xtext=T, ytext=T,
               xgrid=T, ygrid=T, xtick=T, ytick=T)
    p2 = ggplot(tpp) +
        geom_line(mapping = aes(x = recall, y = precision, color = 'prc')) +
        geom_abline(slope = -1, intercept = 1, linetype = 'dotted') +
        scale_x_continuous(name = 'Recall: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
        scale_y_continuous(name = 'Precision: TP/(TP+FP)', breaks=c(.25,.5,.75), limits=c(0,1), expand=c(0,0)) +
        scale_color_manual(values = pal_aaas()(2)[2]) +
        otheme(strip.size = 9, legend.pos = 'none', margin = c(.1,.5,.5,.5),
               xtitle=T, ytitle=T, xtext=T, ytext=T,
               xgrid=T, ygrid=T, xtick=T, ytick=T)
    #
    ggarrange(p1, p2, nrow = 1, ncol = 2, widths = c(1,1))  %>%
        ggexport(filename = fo, width = 10, height = 5)
#}}}
}

resR = roc.curve(scores.class0=to$score, weights.class0=to$weight, curve=T)
resP = pr.curve(scores.class0=to$score, weights.class0=to$weight, curve=T)
scores1 = to %>% filter(weight) %>% pull(score)
scores2 = to %>% filter(!weight) %>% pull(score)
pval = NA
if(length(scores1) > 0 & length(scores2) > 0)
    pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
auroc = resR$auc
auprc = resP$auc.integral

fo = file.path(dirw, "25.auc.pdf")
plot_auc(resR, resP, fo)
#}}}
#}}}

#{{{ gene-centric
eval_tgt <- function(tx) { # score, weight
#{{{
    to = tx
    if(max(to$score) == 0) to$score[1] = 0.1
    scores1 = to %>% filter(weight) %>% pull(score)
    scores2 = to %>% filter(!weight) %>% pull(score)
    pval = NA
    if(length(scores1) > 0 & length(scores2) > 0)
        pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
    pval
#}}}
}
tb = rbind(tc,td) %>% distinct(reg.gid, tgt.gid) %>%
    count(tgt.gid) %>% rename(tid = tgt.gid, nbind = n)
ti = rbind(tc,td) %>% distinct(reg.gid,tgt.gid,weight)
rids = unique(ti$reg.gid)
tids = unique(ti$tgt.gid)
length(rids)

#{{{ target bind
tps = tibble(stat = c('mean','median'),
             val = c(mean(tb$nbind), median(tb$nbind))) %>%
    mutate(lab = str_c(stat, number(val), sep='='))
p = ggplot(tb, aes(x=nbind)) +
    #geom_area(stat='bin', bins=50, color='black', fill=pal_npg()(1), alpha=.2) +
    geom_histogram(color='gray', fill='white', bins=50) +
    #geom_density(aes(y=..count..), alpha=.2, fill=pal_npg()(1)) +
    geom_vline(data=tps, aes(xintercept=val), linetype='dashed') +
    geom_text_repel(data=tps, aes(val, y=100, label=lab), nudge_y=20, size=3) +
    scale_x_continuous(name='# TFs binding to the target (PDI)', expand=expand_scale(mult=c(.03,.03))) +
    scale_y_continuous(name='# target genes', expand=expand_scale(mult=c(.03,.03))) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T)
fo = file.path(dirw, '25.target.nbind.pdf')
ggsave(fo, p, width=8, height=5)
#}}}

tp = crossing(reg.gid = rids, tgt.gid=tids) %>%
    left_join(ti, by=c('reg.gid','tgt.gid')) %>% replace_na(list(weight=F)) %>%
    left_join(tn, by=c('reg.gid','tgt.gid')) %>%
    replace_na(list(score=0))

tp2 = tp %>% mutate(pdi=ifelse(weight, 'py', 'pn')) %>%
    group_by(tgt.gid, pdi) %>%
    summarise(scores = list(score)) %>%
    ungroup() %>% spread(pdi, scores) %>%
    mutate(res = map2(py, pn, wilcox.test, alternative='greater')) %>%
    mutate(pval = map_dbl(res, 'p.value')) %>%
    mutate(padj = p.adjust(pval, 'BH')) %>%
    mutate(logp = -log10(padj))
summary(tp2$logp)

tp3 = tp2 %>% filter(padj < .01)

tpn = tn %>% rename(tid = tgt.gid) %>%
    filter(score > .01) %>%
    count(tid) %>% rename(n_tf = n) %>%
    full_join(tb, by='tid') %>%
    replace_na(list(nbind=0, n_tf=0))

p = ggplot(tpn) +
    geom_point(aes(nbind, n_tf), size=1) +
    scale_x_continuous(name='# PDIs', expand= c(.03,.03)) +
    scale_y_continuous(name='centrality', expand= c(.03,.03)) +
    scale_color_viridis(name='# PDIs', option='viridis', direction=-1) +
    otheme(xtitle=T,ytitle=T, xtext=T,ytext=T, xtick=T, ytick=T,
        legend.pos = 'top.center.out', legend.dir='h', legend.title=T)
fo = file.path(dirw, '25.target.pdf')
ggsave(fo, p, width=7, height=7)

#}}}
#}}}



