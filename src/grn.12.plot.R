source("functions.R")
diri = '~/projects/rnaseq'
dirw = file.path(dird, '14_eval_sum')
gopt = 'rf'

#{{{ general stats
fv = sprintf("%s/%s.100k.rds", dirr, gopt)
ev = readRDS(fv)
eopt = 'tf'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_tf = readRDS(fv)

#{{{ mse stats
tp = ev %>%
    select(nid, mse) %>% unnest() %>%
    inner_join(th, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts0)))
tps = tp %>% group_by(nid) %>%
    summarise(n=n(), n_hq = sum(mse >= .8),
              q25=quantile(mse,.25),
              q50=quantile(mse,.5),
              q75=quantile(mse,.75)) %>% ungroup() %>%
    mutate(lgd2 = sprintf("N=%d", n)) %>%
    inner_join(th, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts0)))

p1 = ggplot(tp) +
    geom_violin(aes(x=lgd, y=mse)) +
    geom_errorbar(data=tps, aes(x=lgd, ymin=q25, ymax=q75), width=.2) +
    geom_point(data=tps, aes(x=lgd, y=q50)) +
    coord_flip() +
    scale_x_discrete(name = '# TFs in network') +
    scale_y_continuous(name = '# Targets in network') +
    scale_color_aaas(name = 'Network type') +
    scale_size(name = 'Median degree per TF') +
    otheme(legend.pos = 'top.left', legend.dir = 'v', legend.title = T,
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T) +
    theme(axis.text.y = element_text(color=rev(nid_cols0)))
fo = file.path(dirw, '03.mse.xgb.pdf')
ggsave(p1, file=fo, width=6, height=8)
#}}}

#{{{ r2 eval stats
tp0 = readRDS(fi_em)
n_gene= 1000
t_oob = ev_tf %>% select(nid, oob) %>% unnest() %>%
    arrange(nid, desc(oob)) %>% group_by(nid) %>%
    filter(row_number() <= n_gene) %>%
    #filter(!gid %in% tsyn$gid) %>%
    ungroup() %>% select(-oob)

th0 = th %>% select(nid,lgd,col)
tp = tp0 %>%
    inner_join(t_oob, by = c('nid','gid')) %>%
    rename(aid=nid, bid=nid_b) %>% group_by(aid, bid) %>%
    summarise(dist = 1 - max(0,quantile(score, .75))) %>%
    ungroup() %>%
    inner_join(th0, by = c('aid'='nid')) %>%
    rename(algd=lgd, acol=col) %>%
    inner_join(th0, by = c('bid'='nid')) %>%
    rename(blgd=lgd, bcol=col)

p = heatmap_hc(tp, leg='Q50 R2 score', top=.4, bottom=5, ratio=4)
fo = file.path(dirw, '03.stat.eval.pdf')
p %>% ggexport(filename = fo, width = 10, height = 8)
#}}}

#{{{ topology stats
tp1 = tp0 %>%
    mutate(deg.reg.q50 = map_dbl(deg.reg, .f <- function(x) median(x$n))) %>%
    select(nid, n.reg, n.tgt, deg.reg.q50) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts)))
#
ymax = 150
tp2 = tp0 %>%
    mutate(res = map(deg.reg, .f <- function(x) desc_stat(x$n))) %>%
    select(nid, res) %>%
    unnest() %>%
    spread(statK, statV) %>%
    mutate(q95 = ifelse(q95 > ymax, ymax, q95)) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts)))

p1 = ggplot(tp1) +
    geom_point(aes(x=n.reg, y=n.tgt, size=deg.reg.q50, color=net_type)) +
    geom_text_repel(aes(x=n.reg, y=n.tgt, label=lgd, color=net_type), size=2.5) +
    scale_x_continuous(name = '# TFs in network') +
    scale_y_continuous(name = '# Targets in network') +
    scale_color_aaas(name = 'Network type') +
    scale_size(name = 'Median degree per TF') +
    otheme(legend.pos = 'top.left', legend.dir = 'v', legend.title = T,
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T)
p2 = ggplot(tp2) +
    geom_errorbar(aes(x=lgd, y=mean, ymin=q5, ymax=q95, color=net_type), linetype='solid', size=.2, width=.4) +
    geom_crossbar(aes(x=lgd, ymin=q25, y=q50, ymax=q75, color=net_type), fill='white', alpha=1, width=.6) +
    scale_x_discrete(expand = c(.02,0)) +
    scale_y_continuous(name = '# Targets per TF', limits=c(0,ymax), expand=expand_scale(mult=c(.01,.03))) +
    coord_flip() +
    scale_color_aaas() +
    otheme(legend.pos = 'top.center.out', legend.dir='h',
        xtick=T, ytick=T,
        xtitle=T, ytitle=F, xtext=T, ytext=T) +
    theme(axis.text.y = element_text(color=rev(nid_cols)))
fo = file.path(dirw, '03.stat.pdf')
ggpubr::ggarrange(p1, p2,
    nrow=1, ncol=2, widths = c(2.5,2), heights = c(1),
    labels=LETTERS[1:2]) %>%
    ggpubr::ggexport(filename = fo, width = 10, height = 8)
#}}}

#{{{ network clustering
tu = ev %>% select(nid, tn) %>% unnest() %>%
    select(nid, reg.gid, tgt.gid, score) %>%
    filter(nid %in% th1$nid) %>%
    group_by(nid) %>% slice(1:1e5) %>% ungroup() %>%
    mutate(reg.tgt = str_c(reg.gid, tgt.gid, sep = '-')) %>%
    select(-reg.gid, -tgt.gid)
reg.tgts = tu %>% count(reg.tgt) %>% filter(n>=2) %>% pull(reg.tgt)
length(reg.tgts)
tu = tu %>% filter(reg.tgt %in% reg.tgts)
#
tuw = tu %>% spread(reg.tgt, score)
#tuw1 = tu %>% spread(reg.tgt, score) %>% replace(is.na(.), 0)
mat = as.matrix(tuw[,-1])
rownames(mat) = tuw$nid
mat[is.na(mat)] = 0
dist_mat = dist(mat, method = 'binary')

th0 = th1 %>% select(nid, lgd, col)
#dis = daisy(tuw[,-1], metric = 'gower')
tp = as.matrix(dist_mat)
colnames(tp) = tuw$nid
tp = as_tibble(tp) %>% mutate(aid = tuw$nid) %>%
    gather(bid, dis, -aid) %>%
    mutate(dis = ifelse(aid==bid, NA, dis)) %>%
    mutate(dis = pmax(.5, dis)) %>%
    rename(dist = dis) %>%
    inner_join(th0, by=c('aid'='nid')) %>% rename(algd=lgd, acol=col) %>%
    inner_join(th0, by=c('bid'='nid')) %>% rename(blgd=lgd, bcol=col)
#
p = heatmap_hc(tp, top=.4, bottom=5, ratio=4)
fo = file.path(dirw, '03.hc.pdf')
p %>% ggexport(filename=fo, width=10, height=8)
#}}}
#}}}

#{{{ TF stats
t_tf0 = gs$all_tf %>% mutate(tf = 'TF')
t_tf = t_tf0 %>% inner_join(tsyn, by = 'gid')

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
fo = file.path(dirw, 'tf.pdf')
ggsave(p1, file=fo, width=4, height=6)
#}}}

ev_tf %>% filter(nid=='nc03') %>% select(nid, tn) %>% unnest() %>%
    select(gid = reg.gid, score) %>% inner_join(t_tf, by='gid') %>%
    group_by(ftype) %>%
    summarise(n = n(), s25=quantile(score,.25), s50=median(score),
              s75 = quantile(score,.75)) %>% ungroup()


t_deg = ev_tf %>% filter(nid=='n17a') %>% select(nid, tn) %>% unnest() %>%
    select(gid = reg.gid, score) %>% count(gid) %>% rename(deg = n)
t_tf %>%
    left_join(t_deg, by = 'gid') %>% replace_na(list(deg = 0)) %>%
    mutate(degbin = cut(deg, breaks=c(0,1,3,10,100,Inf), include.lowest=T, right=F)) %>%
    group_by(ftype, degbin) %>%
    summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
    ungroup() %>% select(-n) %>% spread(degbin, freq)

#}}}

#{{{ ## Y1H eval
dirw = file.path(dird, '08_y1h')
fi = file.path(dirw, '01.rds')
y1h = readRDS(fi)
tn = y1h$t_y1h %>% transmute(reg.gid=reg, tgt.gid=tgt, Y1H='Yes')

tp = ev_tf %>% select(nid, ystat) %>% unnest() %>%
    filter(net_size == 50000) %>% select(-net_size) %>%
    unnest() %>%
    left_join(tn, by = c('reg.gid','tgt.gid'))
tp %>% count(reg.gid)
tp %>% count(tgt.gid)
tp %>% count(reg.gid,tgt.gid)
tp %>% filter(!is.na(Y1H)) %>% count(reg.gid,tgt.gid) %>% count(n)
tp %>% count(reg.gid,tgt.gid) %>% count(n)
tp %>% distinct(reg.gid,tgt.gid,Y1H) %>% count(Y1H)

to = tp %>%
    inner_join(th, by='nid') %>%
    select(reg.gid, tgt.gid, Y1H, txt) %>%
    mutate(txt = factor(txt, levels=th$txt)) %>%
    mutate(s = T) %>%
    spread(txt, s)
fo = file.path(dirw, '11.support.edges.tsv')
write_tsv(to, fo, na='')
#}}}

#{{{ known TF / TFBS AUROC
eopt = 'tf'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_tf = readRDS(fv)

#{{{ obsolete - selected roc/prc plot
nids_tf = c("n16b","n16c","n99b_1","nc03")
cols.dev = c(pal_npg()(4)[2:4], brewer.pal(6,"Paired")[6])
cols.dev = c(pal_npg()(8))

tp1 = tp0 %>% select(nid, ctag, roc) %>%
    filter(nid %in% nids_tf) %>% unnest() %>%
    mutate(ctag = sprintf("%s AUROC", ctag)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts)))
tp2 = tp0 %>% select(nid, ctag, prc) %>%
    filter(nid %in% nids_tf) %>% unnest() %>%
    mutate(ctag = sprintf("%s AUPR", ctag)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels = rev(nid_txts)))
p1 = ggplot(tp1) +
    geom_line(mapping = aes(x = TPR, y = FPR, color = lgd)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'FPR: FP/(FP+TN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'TPR: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    otheme(strip.size = 9, legend.pos = 'none', margin = c(.5,.5,.1,.5),
           xtitle=T, ytitle=T, xtext=T, ytext=T,
           xgrid=T, ygrid=T, xtick=T, ytick=T)
p2 = ggplot(tp2) +
    geom_line(mapping = aes(x = recall, y = precision, color = lgd)) +
    #geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'Recall: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'Precision: TP/(TP+FP)', breaks=c(.25,.5,.75), limits=c(0,1), expand=c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    otheme(strip.size = 9, margin = c(.1,.5,.5,.5),
           xtitle=T, ytitle=T, xtext=T, ytext=T,
           xgrid=T, ygrid=T, xtick=T, ytick=T) +
    theme(legend.position = c(1,.25), legend.justification = c(1,0)) +
    guides(direction = 'vertical', color = guide_legend(ncol = 1, byrow = F))
#
fo = file.path(dirw, "05.roc_prc.dev.pdf")
ggpubr::ggarrange(p1, p2, nrow = 2, ncol = 1, heights = c(1,1))  %>%
    ggpubr::ggexport(filename = fo, width = 10, height = 5)
#}}}

#{{{ use TF mutant RNA-Seq to validate
fd = file.path(dird, '07_mutants', 'degs.rds')
ds = readRDS(fd) %>% filter(gene_alias != 'P1')
dss = ds %>% unnest() %>% group_by(gene_alias, Tissue) %>%
    summarise(n_tot=n(), n_de=sum(padj<.01), prop_de=n_de/n_tot) %>%
    ungroup() %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags = dss$ctag

cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
tv = ev_tf %>% select(nid, ko) %>% unnest() %>% filter(!is.na(pval)) %>%
    inner_join(dss, by=c('gene_alias','Tissue'))
tp = tv %>%
    mutate(lab = ifelse(pval<.05, number(-log10(pval),accuracy=2), '')) %>%
    mutate(pval=-log10(pval)) %>%
    select(nid,ctag,pval,lab) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th1, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(th1$lgd)))
pval.max = max(tp$pval)
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=pval)) +
    geom_tile() +
    geom_text(aes(label=lab, col=abs(pval-pval.max/2)<pval.max*3/8), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="-log10(P-value)", begin = .4) +
    scale_color_manual(values = c('white','black')) +
    scale_fill_gradientn(name = '-log10(P-value)', colors = cols100) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(th1$col)))
fo = sprintf('%s/05.tf.pval.pdf', dirw)
ggsave(p1, file = fo, width = 8, height = 7)
#}}}

#{{{ aupr/auroc bar-plot
levs = c("AUROC", "AUPR")
tp = ev_tf %>% select(nid, ko) %>% unnest() %>%
    inner_join(dss, by=c('gene_alias','Tissue')) %>%
    filter(ctag %in% ctags) %>%
    select(nid,ctag,auroc,auprc) %>%
    mutate(auroc = as.numeric(auroc), auprc = as.numeric(auprc)) %>%
    rename(AUPR = auprc, AUROC = auroc) %>%
    gather(type, auc, -nid,  -ctag) %>%
    mutate(auc = as.numeric(auc)) %>% filter(!is.na(auc)) %>%
    group_by(type) %>% mutate(auc.norm = scale(auc)) %>% ungroup() %>%
    mutate(lab = str_remove(sprintf("%.03f", auc), '^0+')) %>%
    mutate(type = factor(type, levels = levs)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th1, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(th1$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auc.norm)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="Area Under Curve (AUC)", begin = .4) +
    scale_fill_gradientn(name = 'Area Under Curve (AUC)', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    facet_wrap(~type, nrow = 1) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(th1$col)))
fp = sprintf('%s/05.tf.auc.pdf', dirw)
ggsave(p1, file = fp, width = 11, height = 7)
#}}}


#{{{ use TFBS to validate
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
tv = ev_tf %>% select(nid, tfbs) %>% unnest() %>% filter(!is.na(pval))
tp = tv %>%
    mutate(lab = ifelse(pval<.05, number(-log10(pval),accuracy=2), '')) %>%
    mutate(pval=-log10(pval)) %>%
    select(nid,ctag,pval,lab) %>%
    #mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th1, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(th1$lgd)))
pval.max = max(tp$pval)
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=pval)) +
    geom_tile() +
    geom_text(aes(label=lab, col=abs(pval-pval.max/2)<pval.max*3/8), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="-log10(P-value)", begin = .4) +
    scale_color_manual(values = c('white','black')) +
    scale_fill_gradientn(name = '-log10(P-value)', colors = cols100) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(th1$col)))
fo = sprintf('%s/06.tfbs.pval.pdf', dirw)
ggsave(p1, file = fo, width = 5, height = 7)


levs = c("AUROC", "AUPR")
tp = ev_tf %>% select(nid, tfbs) %>% unnest() %>%
    select(nid,ctag,auroc,auprc) %>%
    mutate(auroc = as.numeric(auroc), auprc = as.numeric(auprc)) %>%
    rename(AUPR = auprc, AUROC = auroc) %>%
    gather(type, auc, -nid,  -ctag) %>%
    mutate(auc = as.numeric(auc)) %>% filter(!is.na(auc)) %>%
    group_by(type) %>% mutate(auc.norm = scale(auc)) %>% ungroup() %>%
    mutate(lab = str_remove(sprintf("%.03f", auc), '^0+')) %>%
    mutate(type = factor(type, levels = levs)) %>%
    inner_join(th1, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(th1$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auc.norm)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="Area Under Curve (AUC)", begin = .4) +
    scale_fill_gradientn(name = 'Area Under Curve (AUC)', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    facet_wrap(~type, nrow = 1) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(th1$col)))
fp = sprintf('%s/06.tfbs.auc.pdf', dirw)
ggsave(p1, file = fp, width = 6, height = 6)
#}}}
#}}}

#{{{ evaluate using GO/CornCyc
eopt = 'go'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_go = readRDS(fv)
net_sizes = c(1e4,5e4,1e5,5e5)
net_sizes = c(5e4,5e5)
net_size_map = c('1e4'=1e4, '5e4'=5e4, '1e5'=1e5)

#{{{ enrichment
tp1 = ev_go %>% select(nid, enrich) %>% unnest()
tp2 = ev_go %>% select(nid, enrich_grp) %>% unnest() %>%
    filter(n >= 10) %>%
    group_by(nid, net_size, ctag) %>%
    summarise(n_grp=length(grp), n_grp_sig=sum(pval<.05)) %>%
    ungroup() %>%
    mutate(sigtxt = str_c(n_grp_sig,n_grp,sep='/')) %>%
    select(nid,net_size,ctag,sigtxt)
ctags = c("GO_HC", "GO_arabidopsis","CornCyc")
fo = file.path(dirw, "08.go.pdf")
ctags = c("li2013","liu2017","wang2018")
fo = file.path(dirw, "08.go.eQTL.pdf")
tp = tp1 %>% left_join(tp2, by=c('nid','net_size','ctag')) %>%
    mutate(sig = ifelse(pval<.05, 1, 0)) %>%
    filter(ctag %in% ctags) %>%
    mutate(fc = log2(fc)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th1, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(th1$lgd))) %>%
    filter(net_size %in% net_sizes) %>%
    mutate(net_size=factor(net_size,levels=net_sizes))
tpi = tp %>% filter(sig == 0)
tps = tp %>% distinct(lgd, col) %>% arrange(lgd)
#
p1 = ggplot(tp, aes(lgd, fc)) +
    geom_point(aes(color=net_size, shape=net_size), size=2) +
    geom_point(data=tpi, aes(x=lgd,y=fc), shape=4, size=2) +
    geom_text_repel(aes(label=sigtxt), size=2) +
    #geom_hline(yintercept = 1, alpha= .5, linetype='dotted') +
    scale_x_discrete(expand = expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name = 'log2 Fold Enrichment', expand = c(.05,0)) +
    coord_flip() +
    facet_grid(.~ctag, scale='free') +
    scale_color_npg() +
    scale_shape_manual(values=c(2,1)) +
    scale_alpha(range = c(.4,1), breaks=c(0,1), labels=c("pval > 0.05", "pval < 0.05")) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           strip.size = 8,
           xtitle=T, xtext=T, ytext=T, ygrid=T, xtick=T, ytick=T) +
    theme(axis.text.y = element_text(color=rev(th1$col))) +
    theme(legend.box = "horizontal") +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.3)) +
    guides(color = guide_legend("network size:", nrow=1, order=1),
           shape = guide_legend("network size:", nrow=1, order=1),
           alpha = guide_legend("signifance:", nrow=1, order=2))
    #theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}

#{{{ GO heatmap
gotag = 'GO_uniprot.plants'
goname = gs$fun_ann %>% filter(ctag==gotag) %>% distinct(grp,note) %>%
    mutate(note = str_sub(note, 1, 65))
tx = ev_go %>% select(nid,enrich_grp) %>%
    unnest() %>%
    filter(net_size==5e4,ctag==gotag) %>%
    select(-net_size, -ctag) %>%
    filter(n>=50, pval<.01, fc>=2) %>%
    mutate(fc = log2(fc)) %>%
    inner_join(th1, by = 'nid')

tx %>% inner_join(th1, by = 'nid') %>%
    count(grp,net_type) %>% spread(net_type,n) %>%
    select(-ril,-liftover) %>% inner_join(goname, by='grp') %>%
    print(n=40)

#{{{ #targets per TF histogram
tps = tp %>% count(grp)
p = ggplot(tps, aes(x=nn)) +
    geom_histogram() +
    scale_x_continuous(name = '# Enriched GRNs', expand=c(0,0)) +
    scale_y_continuous(name = '# GO Terms', expand=expand_scale(mult=c(0,.05))) +
    scale_color_aaas() +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '14.2.dist.pdf')
ggsave(p, file=fo, width=6, height=4)
#}}}

require(cluster)
require(ape)
grps = tx %>% count(grp) %>% filter(n>=8) %>% pull(grp)
length(grps)
tp = tx %>% filter(grp %in% grps)

#{{{ GO term order
e = tp %>% select(nid,grp,fc) %>% spread(grp, fc) %>% select(-nid)
e[is.na(e)] = NA
dim(e)
#
cor_opt = "pearson"
hc_opt = "ward.D"
# edist <- as.dist(1-cor(e, method = cor_opt))
edist = daisy(t(e), metric = 'gower')
edist[is.na(edist)] = 0
ehc = hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
go_names = ehc$labels[ehc$order]
#}}}
#{{{ GRN order
e = tp %>% select(nid,grp,fc) %>% spread(nid, fc) %>% select(-grp)
e[is.na(e)] = NA
dim(e)
#
cor_opt = "pearson"
hc_opt = "ward.D"
#edist <- as.dist(1-cor(e, method = cor_opt))
edist = daisy(t(e), metric = 'gower')
edist[is.na(edist)] = 0
ehc <- hclust(edist, method = hc_opt)
ntree = as.phylo(ehc)
s_nids = ehc$labels[ehc$order]
s_ncfg = th1 %>% mutate(nid=factor(nid, levels=s_nids)) %>% arrange(nid)
s_lgds = s_ncfg %>% pull(lgd)
s_cols = s_ncfg %>% pull(col)
#}}}

#{{{
tp = tp %>% mutate(lgd = factor(lgd, levels = s_lgds)) %>%
    mutate(grp = factor(grp, levels = go_names))
tpg = tp %>% distinct(grp) %>% inner_join(goname,by='grp') %>%
    mutate(grp = factor(grp, levels = go_names)) %>%
    arrange(grp)
tpn = tp %>% distinct(nid,lgd,col) %>% rename(taxa=nid) %>% arrange(lgd)
p1 = ggplot(tp) +
    geom_tile(aes(x = lgd, y = grp, fill = fc)) +
    #geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    #geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    #geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_discrete(position='bottom', expand=c(0,0)) +
    scale_y_discrete(breaks=tpg$grp, labels=tpg$note, expand = c(0,0)) +
    scale_fill_viridis(name='log2FC', direction=-1) +
    otheme(legend.pos='right', legend.dir='v', margin=c(.2,.2,.2,.2),
           xtick=T,ytick=T,ytext = T, xgrid=T,ygrid=T) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x=element_text(color=s_cols,angle=45,size=7,hjust=1,vjust=1)) +
    theme(legend.title = element_text())
#
require(ggdendro)
margin = c(0,3.2,0,16.5)
p2 = ggdendrogram(ehc, labels=F, leaf_labels=F, theme_dendro=F) +
    theme_void() +
    theme(plot.margin = unit(margin, "lines"))
fo = file.path(dirw, '09.go.heat.pdf')
ggpubr::ggarrange(p2, p1,
    nrow=2, ncol=1, heights = c(1,11)) %>%
    ggpubr::ggexport(filename = fo, width = 10, height = 12)
#}}}

fo = file.path(dirw, "09.gotree.1.pdf")
pdf(fo, width=7, height=3)
dev.off()

p1 = ggtree(ntree, layout = 'rectangular') +
    #geom_tiplab(color = 'black') +
    scale_x_continuous(expand = expand_scale(.9,.03)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tpn + geom_tiplab(aes(label=lgd)) +
    scale_color_aaas()
fo = file.path(dirw, "09.gotree.1.pdf")
ggsave(p1, filename = fo, width=7, height=8)
#}}}
#}}}

#{{{ case study # obsolete
t_link = gs$tf %>% filter(reg.gid == gid) %>% select(reg.gid,tgt.gid,binding)

t_link_sup = t_net %>% unnest() %>%
    right_join(t_link, by = c("reg.gid",'tgt.gid'))

gs = read_gs()
gid = gs$tfs %>% filter(ctag=='O2') %>% pull(reg.gid)
tgts = gs$tf %>% filter(reg.gid==gid) %>% pull(tgt.gid)
gcfg$gene.desc %>% filter(id %in% tgts) %>% print(n=40)
tq = tibble(gid = c(gid, tgts), tf = c(T, rep(F,length(tgts)))) %>%
    inner_join(gcfg$gene.desc, by = c('gid'='id')) %>%
    mutate(gname = ifelse(is.na(note1),
                          ifelse(is.na(note2), gid, note2), note1)) %>%
    select(gid, tf, gname)

study='mec03'
study='me16a'
th = rnaseq_sample_meta(study)
tm = rnaseq_sample_cpm(study)
ths = th %>%
    distinct(Tissue,Genotype,Treatment) %>%
    arrange(Tissue,Genotype,Treatment) %>%
    #mutate(txt = sprintf("%s--%s",Tissue,Treatment))
    mutate(txt = Genotype)

tp = tm %>% filter(gid %in% tq$gid) %>%
    select(gid,SampleID,CPM) %>%
    inner_join(th, by='SampleID') %>%
    group_by(gid,Genotype,Tissue,Treatment) %>%
    summarise(CPM = asinh(mean(CPM))) %>% ungroup() %>%
    #mutate(txt = sprintf("%s--%s",Tissue,Treatment)) %>%
    mutate(txt = Genotype) %>%
    inner_join(tq, by='gid') %>%
    mutate(gid = factor(gid, levels=tq$gid)) %>%
    mutate(txt = factor(txt, levels=ths$txt))

p = ggplot(tp) +
    geom_tile(aes(x = gid, y = txt, fill = CPM)) +
    #geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    #geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    #geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_discrete(expand = c(0,0), breaks=tq$gid, labels=tq$gname) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    #scale_fill_viridis() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xtext=T, ytext = T) +
    theme(axis.text.x=element_text(angle=30,size=8,hjust=1,vjust=1)) +
    theme(axis.text.y=element_text(size=7)) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 7))
fo = sprintf("%s/15.cpm.2.pdf", dirw)
ggsave(p, file = fo, width = 8, height = 15)
#
#}}}

#{{{ cross-network consistency
tp9 = ev_tf %>%
    filter(!str_detect(nid, '^np')) %>%
    select(nid, tn) %>% unnest() %>%
    mutate(p.drc = ifelse(pcc < 0, -1, 1)) %>%
    distinct(nid,reg.gid,tgt.gid,p.drc) %>%
    group_by(reg.gid,tgt.gid) %>%
    summarise(nnet = n(), m.p.drc=mean(p.drc)) %>%
    ungroup()

tpx = tp9 %>% filter(nnet >= 3) %>% count(reg.gid) %>% rename(deg = n) %>%
    left_join(tsyn, by = c('reg.gid'='gid'))
tpx %>% group_by(deg, ftype) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% ungroup() %>%
    select(-n) %>% spread(ftype, freq) %>% print(n=50)


types = c("consis. +", "consis. -", "mix of +/-")
tp = tp9 %>%
    filter(nnet >= 2) %>%
    mutate(x = ifelse(nnet >= 10, 10, nnet)) %>%
    mutate(type=ifelse(m.p.drc >= .8, 'consis. +', ifelse(m.p.drc <= -.8, 'consis. -', 'mix of +/-'))) %>%
    mutate(xlab = ifelse(x==10, '10+', x)) %>%
    mutate(type = factor(type, levels=types)) %>%
    count(x, type, xlab)
tps = tp %>% group_by(x, xlab) %>%
    summarise(n = sum(n)) %>% ungroup() %>%
    mutate(lab=sprintf("N=%d",n)) %>% arrange(x)
p1 = ggplot(tp, aes(x=x, y=n)) +
    geom_bar(aes(fill=type), stat='identity', position='fill', alpha=.5, width=.75) +
    scale_x_reverse(name='# shared networks', breaks=tps$x, labels=tps$xlab,
                       sec.axis=sec_axis(trans=~.,breaks=tps$x,labels=tps$lab),
                       expand = c(0,0)) +
    scale_y_continuous(name='# predicted TF-target pairs', expand=expand_scale(mult=c(0,0))) +
    #scale_fill_futurama() +
    scale_fill_aaas() +
    coord_flip() +
    otheme(ytitle=T, xtext=T, ytext=T, ytick=T, ygrid=F,
           legend.pos = 'top.center.out', legend.dir = 'h')
fo = file.path(dirw, '10.cross.net.pdf')
ggpubr::ggarrange(p1,
    nrow=1, ncol=1, heights = c(2,2),
    align='v') %>%
    ggpubr::ggexport(filename = fo, width = 6, height = 4)
#}}}

#{{{ evaluate using briggs data
#{{{ read in
#fi = file.path(dirw, "../13_eval/n17a_br.rds")
#x=readRDS(fi)
ev_br = readRDS(fi_br) %>%
    mutate(p.drc = ifelse(is.na(p.drc), b.drc, p.drc)) %>%
    mutate(tgt.DE = ifelse(tgt.DE == 'non_DE', 'non_DE', 'DE'))
des = c("non_DE","DE1-2","DE2-4","DE4+","SPE")
#
br = read_briggs()
tiss = br$tissues
tiss_exc = c("seedlingroot_11DAS","tasselstem_0DAP","endosperm_14DAP",
    "internode_v12","husk_0DAP","seedlingmeristem_11DAS")
tiss_inc = c('seedlingleaf_11DAS','ear_v14','embryo_27DAP','kernel_14DAP','root_0DAP')
tiss1 = tiss[! tiss %in% tiss_exc]
nids_hc = c('nc03',
            'n13c','n14a','n15a','n16a','n18d','n19a',
            'n17a','n18a_1','n18a_2','n18a_3','n18a_4','n18a_5','n18a_6',
            'n99c',
            'nc04')
#}}}

#{{{ show SPE is better than DE in predicting target DE
tp0 = ev_br %>%
    filter(tissue %in% tiss_inc) %>%
    group_by(nid, tissue, reg.DE) %>%
    summarise(nl = n(), n.reg = length(unique(reg.gid)),
              nl.tgt.de = sum(tgt.DE != 'non_DE' & p.drc == b.drc),
              pl.tgt.de = nl.tgt.de / nl) %>%
    ungroup()
#
tp = tp0 %>%
    mutate(txt = sprintf("%d", nl)) %>%
    mutate(txt = sprintf("%d\n%.0f%%", nl, pl.tgt.de*100)) %>%
    mutate(lab = str_remove(sprintf("%.02f", pl.tgt.de), '^0+')) %>%
    rename(tag = reg.DE) %>%
    mutate(tag = factor(tag, levels = des)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts)))
leg = 'Prop. TF-target pairs showing predicted direction of change'
#
p1 = ggplot(tp, aes(tag, lgd)) +
    geom_tile(aes(fill = pl.tgt.de)) +
    geom_text(aes(label = txt), size=2, lineheight=.8) +
    scale_x_discrete(name = 'TF change', expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_fill_gradientn(name=leg, colors=pal_gradient(reverse=T)) +
    facet_wrap(~tissue, nrow=1) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           strip.size = 8, margin = c(1.3, .1, .1, .1),
           xtitle=T, xtext=T, xtick=T, ytext=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
    theme(axis.text.y = element_text(color = rev(nid_cols))) +
    guides(direction = 'horizontal')
p2 = ggplot(tp, aes(lgd, pl.tgt.de)) +
    geom_point(aes(color=tissue, shape=tissue), size = 1.5) +
    scale_x_discrete(name = '', expand = expand_scale(mult=c(.02,.02))) +
    scale_y_continuous(name = leg, breaks=c(.25,.5,.75), expand = expand_scale(mult=c(0,.02))) +
    scale_color_aaas() +
    scale_shape_manual(values = c(0,1,2,3,6)) +
    coord_flip() +
    facet_wrap(~tag, nrow=1) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           strip.size = 8, margin = c(.8, .1, .1, .1),
           xtitle=T, xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    theme(axis.text.y = element_text(color = rev(nid_cols))) +
    guides(direction = 'horizontal')
fp = sprintf("%s/11.de.spe.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 8)
fp = sprintf("%s/11.de.spe.2.pdf", dirw)
ggsave(p2, filename = fp, width = 10, height = 8)
#}}}

#{{{ SPE fold change - heatmap
tp0 = ev_br %>%
    filter(!tissue %in% tiss_exc) %>%
#    filter(!str_detect(nid, '^np')) %>%
    group_by(nid, tissue, reg.DE) %>%
    summarise(nl = n(), pl.tgt.de = sum(tgt.DE=='DE' & p.drc==b.drc)/nl) %>%
    ungroup() %>%
    inner_join(br$des, by = c("tissue"="Tissue")) %>%
    mutate(fc = pl.tgt.de / propDE)

tp = tp0 %>%
    filter(reg.DE == 'SPE') %>%
    #mutate(fc = ifelse(nl.s >= 10, fc, NA)) %>%
    mutate(txt = sprintf("%d\n%.0f%%", nl, pl.tgt.de*100)) %>%
    mutate(lab = str_remove(sprintf("%.02f", pl.tgt.de), '^0+')) %>%
    mutate(tissue = factor(tissue, levels = tiss1)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts)))
leg = 'Fold Enrichment in Prop. Target DE'
p1 = ggplot(tp, aes(x=tissue, y=lgd)) +
    geom_tile(aes(fill = fc)) +
    geom_text(aes(label = txt), size=2, lineheight=.8, color='black') +
    scale_x_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_fill_gradientn(name=leg, colors = pal_gradient(reverse=T)) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           margin = c(1.3, .1, .1, .1),
           xtitle=T, xtext=T, xtick=T, ytext=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
    theme(axis.text.y = element_text(color = rev(nid_cols))) +
    guides(direction = 'horizontal', fill=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12.br.heat.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 9)

tp0 = tp %>%
    mutate(lgd = factor(lgd, levels=rev(th2$lgd)))
p1 = ggplot(tp0, aes(x=lgd)) +
    geom_point(aes(y=fc, color=tag, shape=tag), size=1) +
    scale_x_discrete(expand = c(.01,0)) +
    scale_y_continuous(name = 'Enrichment (fold change) in Target DE',
                       breaks = c(1,2,4),
                       expand = expand_scale(mult=c(0,.1))) +
    scale_fill_d3() +
    scale_shape_manual(values=pal_shapes()) +
    scale_color_d3() +
    coord_flip() +
    facet_wrap(~tissue, nrow=2) +
    otheme(legend.pos = 'top.center.out',
           margin = c(.2,.2,.2,1),
           strip.size = 9,
           xtitle=T, xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(axis.text.y=element_text(color=rev(th$col))) +
    theme(strip.text.y = element_text(angle=0,hjust=.5)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', color=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12.br.1.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}

#{{{ #compare subsets of networks
nids1 = ncfg %>% filter(net_type == 'tissue') %>% pull(nid)
nids2a = ncfg %>% filter(str_detect(nid, '^n17a')) %>% pull(nid)
nids2b = ncfg %>% filter(str_detect(nid, '^n18a')) %>% pull(nid)
nids2c = ncfg %>% filter(str_detect(nid, '^n99a')) %>% pull(nid)
nids2d = ncfg %>% filter(str_detect(nid, '^n99c')) %>% pull(nid)
nids3 = ncfg %>% filter(net_type =='genotype', !str_detect(nid, '_[1-9]$')) %>% pull(nid)
nids4 = ncfg %>% filter(net_type == 'liftover') %>% pull(nid)
tz = tibble(pid = c('1.tissue',
                    '2a.lin2017','2b.kremling','2c.kaeppler','2d.biomap',
                    '3.geno', '4.liftover'),
    nids = list(nids1, nids2a, nids2b, nids2c, nids2d, nids3, nids4))

for (i in 1:7) {
pid = tz$pid[i]
nids = tz$nids[[i]]
tp1 = tp0 %>% filter(nid %in% nids) %>%
    mutate(tag = factor(reg.DE, levels=rev(des))) %>%
    mutate(tissue = factor(tissue, levels=tiss1)) %>%
    mutate(lab = sprintf("%.02f", fc)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts)))
tps = tp1 %>% distinct(lgd,col) %>% arrange(lgd)
p1 = ggplot(tp1) +
    #geom_jitter(aes(y=fc, color=lgd, shape=lgd), width=.4, size=1.5) +
    geom_tile(aes(tissue,lgd,fill=fc)) +
    geom_text(aes(tissue,lgd,label=lab), size=2, color='black') +
    scale_x_discrete(expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand = expand_scale(mult=c(0,0))) +
    #scale_fill_viridis(name='Fold Change', direction=1) +
    scale_fill_gradientn(name = 'Fold Enrichment in Prop. Target DE', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    facet_grid(tag~.) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           margin = c(1.3, .1, .1, .1), strip.size=8,
           xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
    theme(axis.text.y = element_text(color = tps$col)) +
    guides(direction = 'horizontal', color=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12_br/%s.pdf", dirw, pid)
ggsave(p1, filename = fp, width = 6, height = 8)
}
#}}}

#{{{ #root/endosperm/sam networks
tiss1 = c('endosperm_27DAP','kernel_14DAP')
tiss2 = c("root_0DAP","radicle_root")
tiss3 = c('seedlingmeristem_11DAS','seedlingleaf_11DAS','coleoptile_tip')
tiss3 = c('seedlingleaf_11DAS','coleoptile_tip')
tiss0 = c(tiss1, tiss2, tiss3)

tp = tp0 %>% filter(tissue %in% tiss0) %>%
    mutate(tag = factor(reg.DE, levels = des)) %>%
    mutate(tissue = factor(tissue, levels = tiss0)) %>%
    mutate(lab = sprintf("%.02f", fc)) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts)))
tps = tp %>% distinct(lgd,col) %>% arrange(lgd)
p1 = ggplot(tp) +
    geom_tile(aes(tag,lgd,fill=fc)) +
    geom_text(aes(tag,lgd,label=lab), size=2, color='black') +
    scale_x_discrete(expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand = expand_scale(mult=c(0,0))) +
    #scale_fill_viridis(name='Fold Change', direction=1) +
    scale_fill_gradientn(name = 'Fold Enrichment in Prop. Target DE', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    facet_grid(.~tissue) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           margin = c(1.3, .1, .1, .1), strip.size=8,
           xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
    theme(axis.text.y = element_text(color = tps$col)) +
    guides(direction = 'horizontal', color=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12.br.tis.pdf", dirw, pid)
ggsave(p1, filename = fp, width = 9, height = 7)
#}}}


#{{{ #how often is a link observed in multiple tissues? is it consistent?
de_map = des
names(de_map) = 1:length(des)
tp1 = ev_br %>%
    filter(!tissue %in% tiss_exc) %>%
    filter(nid %in% ncfg$nid) %>%
    filter(reg.DE != 'non_DE', tgt.DE != 'non_DE') %>%
    mutate(tag = as.numeric(factor(reg.DE, levels = des))) %>%
    group_by(nid, reg.gid, tgt.gid) %>%
    summarise(nt = n(),
              tag = max(tag),
              p.drc = mean(p.drc),
              m.drc = mean(b.drc)) %>%
    ungroup() %>%
    mutate(tag = de_map[tag]) %>%
    filter(nt >= 5)

tp = tp1 %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=nid_txts))
tps = tp %>% count(lgd, tag)
p = ggplot(tp) +
    geom_violin(aes(x=tag,y=m.drc,fill=tag), position='dodge') +
    geom_text(data=tps, aes(x=tag, y=1.05, label=n), vjust=0, size=2.5) +
    scale_y_continuous(name = 'TF-target regulatory direction', expand=expand_scale(mult=c(.05,.15))) +
    scale_fill_npg(name = 'TF DE class') +
    facet_wrap(~lgd, scale='free', ncol=5) +
    otheme(ytitle=T, ytext=T, xtick=F, ytick=T, ygrid=T,
           legend.pos='bottom.right', legend.dir='v', legend.title=T,
           strip.size=7)
fo = file.path(dirw, '16.drc.tissue.pdf')
ggsave(p, file=fo, width=8, height=10)
#}}}

#{{{ #consis. +/- TF-target pairs for SPE TFs
tags = c("consis. +", "mixed", "consis. -")
tp = tp1 %>% filter(tag == 'SPE') %>%
    mutate(tag = ifelse(m.drc >= .8, 'consis. +',
                        ifelse(m.drc <= -.8, 'consis. -', 'mixed'))) %>%
    mutate(tag = factor(tag, levels=tags))
tp %>% distinct(reg.gid, tgt.gid, tag) %>% count(tag)
tp %>% distinct(reg.gid, tgt.gid, tag) %>% count(0)
tp %>% count(nid, tag) %>% group_by(tag) %>% summarise(nl.min=min(n), nl.max=max(n), nl.avg=mean(n), nl.sum=sum(n))
tp %>% filter(nid=='n99a') %>% count(tag) %>% mutate(prop=n/sum(n))
#
tp = tp %>%
    count(nid, tag) %>% mutate(nl = n) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts))) %>%
    mutate(x = as.numeric(lgd))
tps = tp %>% group_by(nid) %>% summarise(lab=sprintf("N=%d",sum(nl))) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts))) %>%
    mutate(x = as.numeric(lgd)) %>%
    arrange(x)
p = ggplot(tp) +
    geom_bar(aes(x=x,y=nl,fill=tag), alpha=.5,stat='identity', position='fill') +
    #geom_text(data=tps, aes(x=lgd, y=1.01, label=nl), hjust=0, size=2.5) +
    scale_x_continuous(breaks=tps$x, labels=tps$lgd, expand=c(0,0),
                       sec.axis=sec_axis(trans=~.,breaks=tps$x,labels=tps$lab)) +
    scale_y_continuous(expand=expand_scale(mult=c(.0,.0))) +
    scale_fill_aaas() +
    coord_flip() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T,
           legend.pos = 'top.center.out', legend.dir = 'h') +
    theme(axis.text.y = element_text(color=rev(nid_cols)))
fo = file.path(dirw, '16.spe.pdf')
ggsave(p, file=fo, width=6, height=6)
#}}}

#{{{ if a TF has multiple targets, does it show consis. +/- effect?
de_map = des
names(de_map) = 1:length(des)
tp1 = ev_br %>%
    filter(!tissue %in% tiss_exc) %>%
    filter(reg.DE != 'non_DE', tgt.DE != 'non_DE') %>%
    filter(p.drc==b.drc) %>%
    mutate(tag = as.numeric(factor(reg.DE, levels = des))) %>%
    group_by(nid, reg.gid, tgt.gid) %>%
    summarise(n.tiss = n(),
              tag = max(tag),
              m.drc = mean(b.drc)) %>%
    ungroup() %>%
    mutate(m.drc = ifelse(m.drc < 0, -1, 1)) %>%  #trick for liftover nets only
    mutate(tag = de_map[tag]) %>%
    ungroup()
tp1 %>% count(m.drc)
tp1 %>% count(n.tiss)

tp2 = tp1 %>%
    filter(nid %in% nids_hc) %>%
    filter(tag %in% c('SPE')) %>%
    filter(n.tiss >= 3) %>%
    group_by(nid, reg.gid) %>%
    summarise(n.tgt = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    mutate(tf.type=ifelse(m.drc>=.8,'acti',ifelse(m.drc<=-.8,'repr','mix')))
tp2 %>% count(nid, tf.type)
tp2 %>% filter(n.tgt >= 3) %>% count(nid, tf.type)
#
tp3 = tp1 %>%
    filter(nid %in% nids_hc, tag %in% c("SPE"), n.tiss >= 3) %>%
    group_by(reg.gid, tgt.gid) %>%
    summarise(n.net = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    mutate(pair.type=ifelse(m.drc>=.8,'+',ifelse(m.drc<=-.8,'-','mix')))
tp3 %>% count(m.drc)
tp3 %>% filter(n.net >= 2) %>% count(m.drc)
tp3 %>% group_by(n.net) %>%
    summarise(n.link = n(), p.neg=sum(pair.type=='-')/n()) %>% ungroup()

tp4 = tp3 %>% filter(n.net >= 2) %>%
    group_by(reg.gid) %>%
    summarise(n.tgt = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    mutate(tf.type=ifelse(m.drc>=.8,'acti',ifelse(m.drc<=-.8,'repr','mix')))
tp4 %>% filter(n.tgt >= 3) %>% count(tf.type)
#}}}

val.br = list(tf.tgt=tp3, tf=tp4)
fo = file.path(dirw, '02.valid.br.rds')
saveRDS(val.br, file=fo)

#{{{ #multi-tgt plot
tp4 = tp3 %>%
    filter(n.tgt >= 3) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts))) %>%
    mutate(x = as.numeric(lgd))
tp2 %>% distinct(reg.gid) %>% count()
tp2 %>% count(nid) %>% group_by(0) %>% summarise(n.min=min(n), n.max=max(n))

tp = tp2
tps = tp %>% group_by(nid) %>% summarise(lab=sprintf("N=%d",n())) %>%
    inner_join(ncfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts))) %>%
    mutate(x = as.numeric(lgd)) %>%
    arrange(x)
tpz = tibble(xb=c(-Inf,-Inf), xe=c(Inf,Inf), yb=c(-Inf,.8), ye=c(-.8,Inf),
    type=c('repressor','activator'))
p = ggplot() +
#    geom_text(data=tps, aes(x=lgd, y=1.05, label=n.tf), hjust=0, size=2.5) +
    geom_rect(data=tpz, aes(xmin=xb,xmax=xe,ymin=yb,ymax=ye,fill=type),alpha=.5) +
    geom_jitter(data=tp, aes(x=x,y=m.drc), height=0,width=.4, size=1, shape=4) +
    scale_x_continuous(breaks=tps$x, labels=tps$lgd, expand=c(.02,0),
                       sec.axis=sec_axis(trans=~.,breaks=tps$x,labels=tps$lab)) +
    scale_y_continuous(expand=expand_scale(mult=c(.02,.02))) +
    scale_fill_tron() +
    coord_flip() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, xgrid=T, ygrid=T,
           legend.pos = 'top.center.out', legend.dir = 'h') +
    theme(axis.text.y = element_text(color=rev(nid_cols)))
fo = file.path(dirw, '17.drc.tgt.pdf')
ggsave(p, file=fo, width=6, height=8)
#}}}
#}}}

#{{{ eval using biomap data
#{{{ read
ncfg = th %>% filter(!str_detect(nid, "^n99[abc]")) %>%
   select(nid, net_type, sample_size, col, lgd)
nids = ncfg$nid
nid_txts = ncfg$lgd
nid_cols = ncfg$col
#
ev_bm = readRDS(fi_bm)
tiss1 = c("Seedling","Root","Internode","Leaf","Endosperm")
ev_bm0 = ev_bm

adjustp <- function(td)
    td %>%
        mutate(bm.pval.spc = p.adjust(bm.pval.spc, 'bonferroni')) %>%
        mutate(bm.pval.spe = p.adjust(bm.pval.spe, 'BH'))
ev_bm = ev_bm0 %>%
    group_by(nid, Tissue, simu) %>%
    nest() %>%
    #mutate(data = map(data, adjustp)) %>%
    unnest() %>%
    mutate(bm.spc=ifelse(is.na(bm.pval.spc),NA, ifelse(bm.pval.spc<.05,T,F)),
           bm.spe=ifelse(is.na(bm.pval.spe),NA, ifelse(bm.pval.spc<.05,T,F)))
#
nids_hc = c('nc03',
            'n13c','n14a','n15a','n16a','n18d','n19a',
            'n17a','n18a_1','n18a_2','n18a_3','n18a_4','n18a_5','n18a_6',
            'n18d', 'n19a',
            'nc04')
tz = tibble(
    tissue=tiss1,
    nids = list(
                c('nc03','n13c','n14a','n15a','n18a_5','n19a'), #seedling
                c('nc03','nc04','n15a','n17a','n18a_1','n18d','np18_3'), #root
                c('nc03','n15a','n17a','n18a_2','n18a_4'), #internode
                c('nc01','nc03','n14a','n15a','n17a','n18a_2','n19a'), #leaf
                c('nc02','nc03','n14a','n16a','n18a_3','np18_4') #endosperm
    )
)
#}}}

#{{{ GRN eval
tp0 = ev_bm %>% rename(tissue=Tissue) %>%
    group_by(nid, tissue, simu) %>%
    summarise(size.spc = sum(!is.na(bm.spc)),
              p.sig.spc = sum(bm.spc, na.rm=T)/size.spc,
              size.spe = sum(!is.na(bm.spe)),
              p.sig.spe = sum(bm.spe,na.rm=T)/size.spe) %>%
    ungroup() %>%
    mutate(size = size.spe, p.sig = p.sig.spe) %>%
    mutate(size = size.spc, p.sig = p.sig.spc) %>%
    group_by(nid, tissue) %>%
    summarise(size = size[which(!simu)],
              p = p.sig[which(!simu)],
              fc = p.sig[which(!simu)]/p.sig[which(simu)]) %>%
              #fc = p.sig[which(!simu)]) %>%
    ungroup()
fo = file.path(dirw, '21.bm.heat.spe.pdf')
fo = file.path(dirw, '21.bm.heat.spc.pdf')
#
tp = tp0 %>%
    mutate(txt = sprintf("%d", size)) %>%
    mutate(lab = str_remove(sprintf("%.02f", fc), '^0+')) %>%
    mutate(lab = str_remove(sprintf("%d/%.02f/%.02f", size, p, fc), '^0+')) %>%
    mutate(tissue = factor(tissue, levels = tiss1)) %>%
    inner_join(ncfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(nid_txts)))
#tp %>% print(n=20)
leg = 'Fold Enrich. in Prop. Target Sig. Corr.'
p1 = ggplot(tp, aes(x=tissue, y=lgd)) +
    geom_tile(aes(fill = fc)) +
    geom_text(aes(label = lab), size=2.5, color='black') +
    scale_x_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    #scale_fill_viridis(name='Fold Change', direction=1) +
    scale_fill_gradientn(name = leg, colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           margin = c(1.3, .1, .1, .1),
           xtitle=T, xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    #theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
    theme(axis.text.y = element_text(color = rev(nid_cols))) +
    guides(direction = 'horizontal', fill=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
ggsave(p1, file=fo, width=6, height=7)
#}}}

#{{{ how often is a link observed in multiple networks
tp1 = ev_bm %>%
    mutate(supp = bm.spc) %>%
    #inner_join(unnest(tz), by=c('nid'='nids','Tissue'='tissue')) %>%
    filter(nid %in% nids_hc) %>%
    filter(simu==F, !is.na(supp)) %>%
    group_by(nid,reg.gid,tgt.gid) %>%
    summarise(n.tis = n(),
              tag = ifelse(sum(supp)>=1, 'yes','no'),
              m.drc = mean(p.drc)) %>%
    ungroup()
tp1 %>% count(nid, m.drc)
#
tp2 = tp1 %>%
    filter(n.tis >= 2) %>%
    group_by(reg.gid, tgt.gid) %>%
    summarise(n.net = n(),
              tag = ifelse(sum(tag=='yes')>=1, 'yes', 'no'),
              m.drc = mean(m.drc)) %>%
    #mutate(m.drc = ifelse(m.drc < 0, -1, 1)) %>%
    ungroup() %>%
    mutate(pair.type=ifelse(m.drc>=.8,'+',ifelse(m.drc<=-.8,'-','mix')))
tp2 %>% group_by(tag) %>%
    summarise(p.sigleton = sum(n.net==1)/n()) %>%
    ungroup() %>% print(n=50)
# negative regulation is network specific
tp2 %>% group_by(tag, n.net) %>%
    summarise(p.neg = sum(m.drc<=.8)/n(),
              p.pos = sum(m.drc>=.8)/n()) %>%
    ungroup() %>% print(n=50)
#
tp3 = tp2 %>% filter(tag=='yes') %>% select(-tag)
tp3 %>% count(pair.type)
tp3 %>% filter(n.net>=3) %>% count(pair.type)
#
# TF level
tp4 = tp3 %>%  filter(n.net >= 3) %>%
    group_by(reg.gid) %>%
    summarise(n.tgt = n(), m.drc = mean(m.drc)) %>%
    ungroup() %>%
    mutate(tf.type=ifelse(m.drc>=.8,'acti',ifelse(m.drc<=-.8,'repr','mix')))
tp4 %>% count(tf.type)
tp4 %>% filter(n.tgt >= 3) %>% count(tf.type)

#{{{ GO enrich
go = get_go()

go_pos = go$go %>% filter(ctag=='aggregate') %>% filter(str_detect(goname, 'positive regulation of transcription'))
go_neg = go$go %>% filter(ctag=='aggregate') %>% filter(str_detect(goname, 'negative regulation of transcription'))

reg.gids.pos = tp2 %>% filter(m.drc==1) %>% distinct(reg.gid) %>% pull(reg.gid)
reg.gids.neg = tp2 %>% filter(m.drc==-1) %>% distinct(reg.gid) %>% pull(reg.gid)
go_enrich(reg.gids.pos, go_pos)

reg.gids.neg = tp3 %>% filter(tag.spc=='spc',nnet>=3,m.drc==-1) %>% distinct(reg.gid) %>% pull(reg.gid)
reg.gids.neg = tp4 %>% filter(n.tgt>=3,m.drc<=0) %>% distinct(reg.gid) %>% pull(reg.gid)
go_enrich(reg.gids.neg, go_neg)
#}}}
#}}}


val.bm.spc = list(tf.tgt = tp3, tf = tp4)
val.bm.spe = list(tf.tgt = tp3, tf = tp4)

fo = file.path(dirw, '02.valid.bm.spc.rds')
saveRDS(val.bm.spc, file=fo)
fo = file.path(dirw, '02.valid.bm.spe.rds')
saveRDS(val.bm.spe, file=fo)
#}}}

#{{{ # check & plot high-conf TF-target pairs
tf1 = val.br$tf %>% filter(n.tgt>=3) %>% pull(reg.gid)
tf2 = val.bm.spc$tf %>% filter(n.tgt>=3) %>% pull(reg.gid)
tf3 = val.bm.spe$tf %>% filter(n.tgt>=3) %>% pull(reg.gid)
length(tf1)
length(tf2)
length(tf3)
sum(tf1 %in% tf2)
sum(tf1 %in% tf3)
sum(tf2 %in% tf3)

val = val.bm.spe; pre = 'bm.spe'
val = val.bm.spc; pre = 'bm.spc'
val = val.br; pre = 'br'
tp = val$tf.tgt %>%
    count(n.net, pair.type)
tps = tp %>% mutate(n.pair = n) %>% group_by(n.net) %>%
    summarise(lab=sprintf("%d [N=%d]", n.net[1], sum(n.pair))) %>% ungroup()
p1 = ggplot(tp) +
    geom_bar(aes(x=n.net,y=n,fill=pair.type), stat='identity',position='fill', width=.8) +
    scale_x_continuous(name='# shared networks', breaks=tps$n.net, labels=tps$lab, expand=expand_scale(mult=c(.02,.02))) +
    scale_y_continuous(expand=expand_scale(mult=c(0,0))) +
    scale_fill_npg(name='regulatory direction') +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, xgrid=T, ygrid=T,
           legend.pos = 'top.center.out', legend.dir = 'h', legend.title=T) +
    theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1))
#
tp = val$tf %>% filter(n.tgt >=3)
xlab = sprintf("Mean regulatory direction [N=%d TFs]",nrow(tp))
tpz = tibble(xb=c(-Inf,-Inf), xe=c(Inf,Inf), yb=c(-Inf,.8), ye=c(-.8,Inf),
    type=c('repressor','activator'))
p2 = ggplot() +
#    geom_text(data=tps, aes(x=lgd, y=1.05, label=n.tf), hjust=0, size=2.5) +
    geom_rect(data=tpz, aes(xmin=xb,xmax=xe,ymin=yb,ymax=ye,fill=type),alpha=.5) +
    geom_point(data=tp, aes(x=reg.gid,y=m.drc,size=n.tgt), shape=1) +
    scale_x_discrete(expand=c(.02,0)) +
    scale_y_continuous(name=xlab,limits=c(-1,1),expand=expand_scale(mult=c(.02,.02))) +
    scale_size('# targets') +
    scale_fill_tron(name='direction') +
    coord_flip() +
    otheme(xtext=T, ytext=F, xtitle=T, xtick=T, ytick=F, xgrid=T, ygrid=F,
           legend.pos = 'top.center.out', legend.dir = 'h', legend.title=T,
           margin=c(3,.1,.1,.1)) +
    guides(size=guide_legend(direction='horizontal',title.position='left'))
#
fo = sprintf("%s/23.%s.pdf", dirw, pre)
ggpubr::ggarrange(p1, p2,
    nrow=1, ncol=2, widths = c(2,2), heights = c(1),
    labels=LETTERS[1:2]) %>%
    ggpubr::ggexport(filename = fo, width = 8, height = 5)
#}}}

t_size = gcfg$chrom %>% select(chrom,size=end)
tx = gcfg$chrom
gstart = flattern_gcoord(tx %>% select(chrom,pos=start), t_size)
gend = flattern_gcoord(tx %>% select(chrom,pos=end), t_size)
tx = tx %>% mutate(start=gstart, end=gend, pos=(start+end)/2)
#
nids_hc = c('nc03',
            'n13c','n14a','n15a','n16a','n18d','n19a',
            'n17a','n18a_1','n18a_2','n18a_3','n18a_4','n18a_5','n18a_6',
            'nc04', 'n13a', 'n18c')

#{{{ check for overlap w. trans- hotspots
fun_ann_note = gs$fun_ann %>% distinct(ctag,grp,note)
t_gl = gcfg$loc.gene %>% group_by(gid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end)) %>%
    ungroup() %>%
    mutate(pos=(start+end)/2) %>% select(gid,chrom,pos)

qtags = c('li2013','liu2017','wang2018')
hs = tibble(qtag=qtags) %>%
    mutate(fi=sprintf('~/projects/genome/data2/%s/10.rds', qtag)) %>%
    mutate(data=map(fi, readRDS)) %>%
    mutate(data=map(data, 'hs')) %>%
    select(qtag, data) %>% unnest() %>%
    select(qtag,qid,qchrom,qpos,n.tgt)

nid = 'n18a_6'
nid = 'n13a'

#{{{
ztag = 'GO_uniprot.plants'
ztag = 'CornCyc'
tfid = 'Zm00001d026147'
tgid = 'Zm00001d017077'
tgid = 'Zm00001d026141'
tz2 = ev_go %>% filter(nid %in% nids_hc) %>%
    select(nid, enrich_reg) %>% unnest() %>%
    filter(ctag == ztag, n >= 5, net_size==50000, pval < .05) %>%
    select(nid, reg.gid, n, fc, grp=max.grp, max.grp.size) %>%
    count(reg.gid, grp) %>% filter(n>=1) %>%
    #arrange(grp,desc(fc)) %>%
    filter(reg.gid == tfid) %>%
    inner_join(fun_ann_note, by='grp')
tno = ev_tf %>% filter(nid==!!nid) %>% pull(tn)
tno[[1]] %>% filter(reg.gid==tfid) %>% select(tgt.gid) %>%
    inner_join(gcfg$gene.desc, by=c('tgt.gid'='id')) %>%
    select(-note2) %>% arrange(tgt.gid) %>% print(n=30)

nid = 'n18a_5'
fem = sprintf("%s/../11_exp_mat/%s.tsv", dirw, nid)
tem1 = read_tsv(fem)
nid = 'nc03'
fem = sprintf("%s/../11_exp_mat/%s.tsv", dirw, nid)
tem2 = read_tsv(fem)
nid = 'n18a_6'
fem = sprintf("%s/../11_exp_mat/%s.tsv", dirw, nid)
tem3 = read_tsv(fem)
nid = 'n18c'
fem = sprintf("%s/../11_exp_mat/%s.tsv", dirw, nid)
tem4 = read_tsv(fem)
nid = 'n13a'
fem = sprintf("%s/../11_exp_mat/%s.tsv", dirw, nid)
tem5 = read_tsv(fem)

tem = tem1
cor.test(as.numeric(tem[tem$gid==tfid,]), as.numeric(tem[tem$gid==tgid,]))
tem = tem2
cor.test(as.numeric(tem[tem$gid==tfid,]), as.numeric(tem[tem$gid==tgid,]))
tem = tem3
cor.test(as.numeric(tem[tem$gid==tfid,]), as.numeric(tem[tem$gid==tgid,]))
tem = tem4
cor.test(as.numeric(tem[tem$gid==tfid,]), as.numeric(tem[tem$gid==tgid,]))
tem = tem5
cor.test(as.numeric(tem[tem$gid==tfid,]), as.numeric(tem[tem$gid==tgid,]))
#}}}

nids_hc = nid
tz = ev_go %>% filter(nid %in% nids_hc) %>%
    select(nid, enrich_reg) %>% unnest() %>%
    filter(ctag %in% qtags, n >= 10, net_size==50000, pval < .05) %>%
    select(ctag, nid, reg.gid, n, fc, grp=max.grp, max.grp.size) %>%
    #count(reg.gid, grp) %>% filter(n>=1) %>%
    arrange(ctag,grp,desc(fc)) %>%
    left_join(t_gl, by=c('reg.gid'='gid')) %>%
    rename(gchrom=chrom,gpos=pos) %>%
    left_join(hs, by=c('ctag'='qtag','grp'='qid')) %>%
    mutate(hit=ifelse(gchrom==qchrom & abs(gpos-qpos)<=5e7,'hit','non-hit'))
#tz %>% select(-grp,-max.grp.size) %>% print(n=50)
tz0 = tz

#{{{ plot
ti = tz
gpos = flattern_gcoord(ti %>% select(chrom=gchrom, pos=gpos), t_size)
qpos = flattern_gcoord(ti %>% select(chrom=qchrom, pos=qpos), t_size)
tp = ti %>% mutate(gpos=!!gpos, qpos=!!qpos)
tps = tibble(reg.gid='Zm00001d026147', gname='R1')
tps = tps %>% inner_join(tp, by = 'reg.gid')
#
p1 = ggplot(tp) +
    geom_point(aes(x=qpos, y=gpos, color=max.grp.size, shape=hit), size=1) +
    geom_text_repel(data=tps, aes(qpos,gpos,label=gname), nudge_x=-150000000, direction='y', segment.size=.2, size=3) +
    geom_vline(xintercept = tx$start, alpha=.1) +
    geom_vline(xintercept = tx$end, alpha=.1) +
    geom_hline(yintercept = tx$start, alpha=.1) +
    geom_hline(yintercept = tx$end, alpha=.1) +
    #geom_abline(intercept = 0, slope = 1, alpha=.1) +
    scale_x_continuous(name='trans-eQTL hotsplot position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_y_continuous(name='enriched TF position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_shape_manual(values=c(16,4), guide=F) +
    facet_wrap(~ctag,nrow=1) +
    scale_color_viridis(name='N_targets',option = "plasma") +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, legend.title=T,
         legend.pos='top.center.out', legend.dir='h')
fp = sprintf("%s/32.hs.pdf", dirw, qtag)
ggsave(p1, filename = fp, width = 12, height = 4.5)
#}}}

# save to file
to = tz %>% filter(hit == 'hit') %>%
    group_by(reg.gid) %>%
    summarise(fc = max(fc), n_qtag = n(),
              max.grp.size = max(max.grp.size),
              qtags = paste(ctag, collapse=',')) %>% ungroup() %>%
    filter(max.grp.size >= 5) %>%
    arrange(desc(n_qtag), desc(fc), desc(max.grp.size))
to %>% print(n=50)
fo = file.path(dirw, '02.hs.tsv')
write_tsv(to, fo)

tz = ev_go %>% filter(nid==!!nid) %>%
    select(enrich_grp) %>% unnest() %>%
    filter(ctag == 'li2013', n >=10, net_size==50000, pval < .01) %>%
    select(grp, n, fc, reg.gid=max.reg.gid, max.reg.size) %>%
    arrange(desc(fc)) %>%
    left_join(hs, by=c('grp'='qid')) %>%
    left_join(t_gl, by=c('reg.gid'='gid')) %>% rename(gchrom=chrom,gpos=pos)
tz %>% select(-grp,-max.reg.size) %>% print(n=50)

tnl = ev_tf %>% filter(nid==!!nid) %>% pull(tn)
tx0 = tnl[[1]] %>%
    count(reg.gid) %>% rename(n.tgt=n)
gid = 'Zm00001d023987' #RBP
gid = 'Zm00001d028842' #p1
gid = 'Zm00001d028851'
tx0 %>% filter(reg.gid == gid)

#tx1 = val.br$tf %>% mutate(ctag='br')
#tx2 = val.bm.spe$tf %>% mutate(ctag='bm.spe')
#tx3 = val.bm.spc$tf %>% mutate(ctag='bm.spc')
tp = rbind(tx1,tx2,tx3) %>% filter(n.tgt>=10) %>% inner_join(t_gl, by=c('reg.gid'='gid'))
#
p1 = ggplot() +
    geom_rect(data=hs, aes(xmin=start,xmax=end,ymin=0,ymax=1), fill='steelblue') +
    geom_point(data=hs, aes(x=qpos,y=.5),shape=3,col='steelblue') +
    geom_point(data=tp, aes(x=pos,y=.5,size=n.tgt),shape=4) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    facet_grid(chrom~.) +
    otheme(xtitle=T,xtext=T,xgrid=T)
fo = file.path(dirw, '25.pdf')
ggsave(p1, file=fo, width=8,height=6)

chrom='B10'
t_hs %>% filter(chrom==!!chrom)
tp %>% filter(chrom==!!chrom) %>% print(n=30)
hs='hs65'
hs.tgts = t_hs0 %>% filter(mid==hs, type=='trans') %>% pull(gid)
reg.gid = 'Zm00001d024894'
br.tgts = val.br$tf.tgt %>% filter(n.net>=2) %>% filter(reg.gid == !!reg.gid) %>% pull(tgt.gid)
net.tgts = tnl[[1]] %>% filter(reg.gid == !!reg.gid) %>% pull(tgt.gid)

sum(hs.tgts %in% net.tgts)
#}}}

fun_ann = gs$fun_ann %>% distinct(ctag, grp, note)
gotag = 'CornCyc'
gotag = 'GO_arabidopsis'
gotag = 'GO_uniprot.plants'
gotag = 'li2013'

nid='n13a'
te = ev_go %>% filter(nid==!!nid) %>% pull(enrich_term)
te = te[[1]]
te %>% filter(net_size==50000, ctag==!!gotag, pval<.05,n>=10) %>%
    inner_join(fun_ann, by=c('ctag','grp')) %>%
    distinct(n,fc,pval,note,grp) %>% arrange(desc(fc)) %>% print(n=20)

grp='qtl2862'
gs$fun_ann %>% filter(ctag==!!gotag, grp==!!grp)


res$hs %>% filter(qid==!!grp)
