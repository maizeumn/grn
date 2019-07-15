source("functions.R")
require(pROC)
diri = '~/projects/rnaseq'
dirw = file.path(dird, '14_eval_sum')
gopt = 'rf'

#{{{ tables
#{{{ GRN sum
tp = t_cfg %>% select(net_type,nid,study,note,sample_size,col) %>%
    mutate(study=str_replace(str_to_title(study), '(\\d+)$', ' et al. \\1')) %>%
    arrange(net_type, nid) %>%
    select(-nid)
x = tp %>% select(-col) %>%
    mutate(net_type = cell_spec(net_type, color=tp$col),
           study = cell_spec(study, color=tp$col),
           note = cell_spec(note, color=tp$col),
           sample_size = cell_spec(sample_size, color=tp$col)) %>%
    rename(`network type`=net_type, N=sample_size) %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T,
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left') %>%
    column_spec(1, bold=T) %>%
    collapse_rows(columns=c(1,2), latex_hline="major", valign="top")
fo = file.path(dirw, '01.t1.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t1.rds')
saveRDS(x, file=fo)
#}}}

#{{{ TF mutant
fi = '~/projects/barn/data/01.cfg.xlsx'
tp = read_xlsx(fi, 'mutants') %>%
    mutate(accession=str_sub(accession, 1, 11)) %>%
    replace_na(list(author='')) %>%
    select(`TF alias`=gene_alias, `TF name`=gene_name, `TF ID`=gene_id,
           Study=author, Accession=accession, Tissue=tissue, N=n)
x = tp %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirw, '01.t2.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t2.rds')
saveRDS(x, file=fo)
#}}}

#{{{ natural variation DEGs
fi = '~/projects/master.xlsx'
th = read_xlsx(fi, 'barn')
fi = file.path(dirw, '../06_deg/all.rds')
t_ds = readRDS(fi)

tp = t_ds %>% count(yid, cond, group1, group2, DE) %>% spread(DE, n) %>%
    rename(`non-DE` = non_DE) %>%
    mutate(cond = str_replace(cond, "_", " ")) %>%
    filter(!str_detect(group2,'x'))
tp = th %>% select(yid, author, study) %>%
#    mutate(author=str_replace(str_to_title(author), '(\\d+)$', ' et al. \\1')) %>%
    inner_join(tp, by='yid') %>% select(-yid) %>%
    unite('contrast', group1, group2, sep=' vs ') %>%
    rename(condition = cond)
x = tp %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 8, position='left') %>%
    collapse_rows(columns=c(1,2), latex_hline="major", valign="top")
fo = file.path(dirw, '01.t3.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t3.rds')
saveRDS(x, file=fo)
#}}}

#|{{{ trans hotspots
fi = file.path(dird, '14_eval_sum/33.hs.tsv')
tp = read_tsv(fi) %>%
    replace_na(list(reg.note='',txt='')) %>%
    mutate(reg.note = str_replace(reg.note, " *\\[.*\\]", '')) %>%
    mutate(reg.note = str_replace(reg.note, ";.*$", '')) %>%
    select(ID=reg.gid, `Support eQTL study`=qtags, `Support GRN`=studies, `TF Annotation`=reg.note, `Target enrichment`=txt)
x = tp %>%
    kable(format='latex', escape=F, longtable=T, booktabs=T,
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=T, font_size=8, position='left') %>%
    #column_spec(1:3, width = "1.85cm") %>%
    column_spec(5, width = "5cm")
fo = file.path(dirw, '01.t9.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t9.rds')
saveRDS(x, file=fo)
#}}}
#}}}

#{{{ general stats
fv = sprintf("%s/%s.100k.rds", dirr, gopt)
ev = readRDS(fv)
eopt = 'tf'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_tf = readRDS(fv)

#{{{ #mse stats
tp = ev %>%
    select(nid, mse) %>% unnest() %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
tps = tp %>% group_by(nid) %>%
    summarise(n=n(), n_hq = sum(mse >= .8),
              q25=quantile(mse,.25),
              q50=quantile(mse,.5),
              q75=quantile(mse,.75)) %>% ungroup() %>%
    mutate(lgd2 = sprintf("N=%d", n)) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))

p1 = ggplot(tp) +
    geom_violin(aes(x=lgd, y=mse)) +
    geom_errorbar(data=tps, aes(x=lgd, ymin=q25, ymax=q75), width=.2) +
    geom_point(data=tps, aes(x=lgd, y=q50)) +
    coord_flip() +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'MSE') +
    scale_color_aaas(name = 'Network type') +
    otheme(legend.pos = 'top.left', legend.dir = 'v', legend.title = T,
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fo = file.path(dirw, '02.1.mse.pdf')
ggsave(p1, file=fo, width=6, height=8)
#}}}

#{{{ #r2 eval stats
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
fo = file.path(dirw, '02.1.stat.eval.pdf')
p %>% ggexport(filename = fo, width = 10, height = 8)
#}}}

#{{{ #topology stats
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
    filter(nid %in% t_cfg$nid) %>%
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

th0 = t_cfg %>% select(nid, lgd, col)
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
p = heatmap_hc(tp, top=.4, bottom=4.5, ratio=4)
fo = file.path(dirw, '03.hc.pdf')
p %>% ggexport(filename=fo, width=10, height=8)
#}}}

#{{{ tSNE
require(Rtsne)
tt = tu %>% spread(nid, score) %>% replace(., is.na(.), 0)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 2000)

tp = as_tibble(tsne$Y) %>%
    add_column(nid = colnames(tt)[-1]) %>%
    inner_join(t_cfg, by = 'nid')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=lgd), size=2.5) +
    geom_point(aes(x=V1, y=V2, color=net_type, shape=net_type), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_color_aaas() +
    scale_shape_manual(values = c(0:7)) +
    otheme(legend.pos='bottom.right', legend.dir='v', legend.title=F, legend.border=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
fp = file.path(dirw, "03.tsne.pdf")
ggsave(p_tsne, filename = fp, width=6, height=6)
#}}}
#}}}

#{{{ TF & target stats
diro = file.path(dird, '17_degree')
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

# get TF & target stats
fi = sprintf("%s/%s.500k.rds", dirr, gopt)
tn = readRDS(fi) %>% mutate(tn2 = map(tn, bin_network, bins=10)) %>%
    select(nid, tn2) %>% unnest()

deg.tgt = tn %>% group_by(nid, tgt.gid) %>%
    summarise(r1 = sum(score >= 1), r2 = sum(score >= 2),
              r3 = sum(score >= 3), r4 = sum(score >= 4)) %>%
    ungroup() %>%
    gather(score, deg, -nid, -tgt.gid) %>%
    mutate(score = as.integer(str_replace(score, 'r', '')))
deg.reg = tn %>% group_by(nid, reg.gid) %>%
    summarise(r1 = sum(score >= 1), r2 = sum(score >= 2),
              r3 = sum(score >= 3), r4 = sum(score >= 4)) %>%
    ungroup() %>%
    gather(score, deg, -nid, -reg.gid) %>%
    mutate(score = as.integer(str_replace(score, 'r', '')))

fo = file.path(diro, 'degree.rds')
saveRDS(list(tgt=deg.tgt, reg=deg.reg), file=fo)
#}}}

#{{{ # [obsolete] Y1H eval
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

#{{{ TF KO / TFBS AUROC
eopt = 'tf'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_tf = readRDS(fv)
#{{{ read ko & networks
fd = file.path(dird, '07_mutants', 'degs.rds')
ds = readRDS(fd) %>% filter(gene_alias != 'P1')
dss = ds %>% unnest() %>% group_by(gene_alias, Tissue) %>%
    summarise(n_tot=n(), n_de=sum(padj<.01), prop_de=n_de/n_tot) %>%
    ungroup() %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags = dss$ctag
#
fv = sprintf("%s/%s.1m.rds", dirr, gopt)
ev = readRDS(fv)
ko = gs$ko %>% mutate(kid=1:n()) %>%
    select(kid,gene_id,gene_alias,Tissue,ds) %>%
    unnest() %>%
    mutate(response=ifelse(padj < .01, 1, 0)) %>%
    rename(reg.gid = gene_id, tgt.gid = gid) %>%
    select(-disp,-log2fc)
rids_ko = gs$ko %>% distinct(gene_id) %>% pull(gene_id)
#}}}
#
#{{{ functions
complete_tn <- function(tn, tids, rids) {
    #{{{
    rids = rids[rids %in% tn$reg.gid]
    crossing(reg.gid = rids, tgt.gid = tids) %>% as_tibble() %>%
        left_join(tn, by=c('reg.gid','tgt.gid')) %>%
        replace_na(list(score=0))
    #}}}
}
get_auroc <- function(ti, fpr=1) {
    #{{{
    if(max(ti$response) == 0)
        NA
    else if (min(ti$response) == 1)
        NA
    else
        roc(ti$response, ti$score, partial.auc=c(1,1-fpr), levels=c(0,1))$auc
    #}}}
}
get_spc_pval <- function(ti)
    -log10(cor.test(ti$score, ti$padj, method='kendall')$p.value)
get_wil_pval <- function(ti) {
    #{{{
    scores1 = ti %>% filter(response==1) %>% pull(score)
    scores2 = ti %>% filter(response==0) %>% pull(score)
    if(length(scores1) > 0 & length(scores2) > 0)
        -log10(wilcox.test(scores1, scores2, alternative='greater')$p.value)
    else
        NA
    #}}}
}
#}}}

#{{{ ko auroc + pval
tp0 = ev %>% select(nid,tids,tn) %>%
    mutate(res = map2(tn, tids, complete_tn, rids = rids_ko)) %>%
    select(nid, res) %>% unnest() %>%
    select(nid, reg.gid, tgt.gid, score) %>%
    inner_join(ko, by=c('reg.gid','tgt.gid')) %>%
    group_by(nid, kid, gene_alias, Tissue) %>%
    nest()
tp1 = tp0 %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval))

tp = tp1 %>% #mutate(lab = str_remove(sprintf("%.04f", auroc), '^0+')) %>%
    mutate(lab = number(auroc, accuracy=.001)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
#    mutate(tcol = ifelse(auroc > .0, 'black', 'white')) %>%
    inner_join(dss, by=c('gene_alias','Tissue')) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auroc)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), color='black', hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name = 'Area Under Curve (AUC)') +
    scale_fill_gradientn(name = 'AUC (FPR<=0.1)', colors=cols100) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/05.tf.auc.pdf', dirw)
ggsave(p1, file = fp, width = 8, height = 7)
#}}}
#{{{ pval
tp = tp1 %>% #mutate(lab = str_remove(sprintf("%.04f", auroc), '^0+')) %>%
    mutate(lab = number(spc.pval, accuracy=1)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
#    mutate(tcol = ifelse(auroc > .0, 'black', 'white')) %>%
    inner_join(dss, by=c('gene_alias','Tissue')) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=spc.pval)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), color='black', hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name = 'Area Under Curve (AUC)') +
    scale_fill_gradientn(name = 'Wilcox test p-value (-log10)', colors=cols100) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/05.tf.pval.pdf', dirw)
ggsave(p1, file = fp, width = 8, height = 7)
#}}}

#{{{ TFBS
rids_tfbs = gs$tfbs %>% distinct(reg.gid) %>% pull(reg.gid)
ctags_tfbs = gs$tfbs %>% distinct(ctag) %>% pull(ctag)
tfbs = gs$tfbs %>% mutate(response=1)
tps = tfbs %>% count(ctag) %>% mutate(ctagl=sprintf("%s [%s]", ctag, number(n)))
tp0 = ev %>% select(nid,tids,tn) %>%
    #mutate(res = map2(tn, tids, complete_tn, rids = rids_tfbs)) %>%
    select(nid, tn) %>% unnest() %>%
    filter(reg.gid %in% rids_tfbs) %>%
    crossing(ctag = ctags_tfbs) %>%
    select(nid, ctag, reg.gid, tgt.gid, score) %>%
    left_join(tfbs, by=c('ctag','reg.gid','tgt.gid')) %>%
    replace_na(list(response = 0)) %>%
    group_by(nid, ctag) %>%
    nest()
tp1 = tp0 %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval))

tp = tp1 %>% #mutate(lab = str_remove(sprintf("%.04f", auroc), '^0+')) %>%
    mutate(lab = number(auroc, accuracy=.001)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auroc)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), color='black', hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0)), breaks=tps$ctag, labels=tps$ctagl) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(name = 'AUC (FPR<=0.1)', colors=cols100) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/06.tfbs.auc.pdf', dirw)
ggsave(p1, file = fp, width = 4.5, height = 7)
#}}}
#{{{ pval
tp = tp1 %>% #mutate(lab = str_remove(sprintf("%.04f", auroc), '^0+')) %>%
    mutate(lab = number(spc.pval, accuracy=1)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
#    mutate(tcol = ifelse(auroc > .0, 'black', 'white')) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=spc.pval)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), color='black', hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0)), breaks=tps$ctag, labels=tps$ctagl) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name = 'Area Under Curve (AUC)') +
    scale_fill_gradientn(name = 'Wilcox test p-value (-log10)', colors=cols100) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/06.tfbs.pval.pdf', dirw)
ggsave(p1, file = fp, width = 4.5, height = 7)
#}}}

#{{{ # [obsolete] - selected roc/prc plot
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
#{{{ # [obsolete] - use TF mutant RNA-Seq to validate
tv = ev_tf %>% select(nid, ko) %>% unnest() %>% filter(!is.na(pval)) %>%
    inner_join(dss, by=c('gene_alias','Tissue'))
tp = tv %>%
    mutate(lab = ifelse(pval<.05, number(-log10(pval),accuracy=2), '')) %>%
    mutate(pval=-log10(pval)) %>%
    select(nid,ctag,pval,lab) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
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
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fo = sprintf('%s/05.tf.pval.pdf', dirw)
ggsave(p1, file = fo, width = 8, height = 7)
#}}}
#{{{ # [obsolete] aupr/auroc bar-plot
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
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
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
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/05.tf.auc.pdf', dirw)
ggsave(p1, file = fp, width = 11, height = 7)
#}}}
#{{{ # [obsolete] TFBS
tv = ev_tf %>% select(nid, tfbs) %>% unnest() %>% filter(!is.na(pval))
tp = tv %>%
    mutate(lab = ifelse(pval<.05, number(-log10(pval),accuracy=2), '')) %>%
    mutate(pval=-log10(pval)) %>%
    select(nid,ctag,pval,lab) %>%
    #mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
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
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
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
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels = rev(t_cfg$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auc.norm)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="Area Under Curve (AUC)", begin = .4) +
    scale_fill_gradientn(name = 'Area Under Curve (AUC)', colors = cols100) +
    facet_wrap(~type, nrow = 1) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fp = sprintf('%s/06.tfbs.auc.pdf', dirw)
ggsave(p1, file = fp, width = 6, height = 6)
#}}}
#}}}

#{{{ GO/CornCyc
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
ctags = c("li2013","liu2017","wang2018")
fo = file.path(dirw, "08.go.eQTL.pdf")
ctags = c("GO_HC", "GO_arabidopsis","CornCyc")
fo = file.path(dirw, "08.go.pdf")
tp = tp1 %>% left_join(tp2, by=c('nid','net_size','ctag')) %>%
    mutate(sig = ifelse(pval<.05, 1, 0)) %>%
    filter(ctag %in% ctags) %>%
    mutate(fc = log2(fc)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd))) %>%
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
    theme(axis.text.y = element_text(color=rev(t_cfg$col))) +
    theme(legend.box = "horizontal") +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.3)) +
    guides(color = guide_legend("network size:", nrow=1, order=1),
           shape = guide_legend("network size:", nrow=1, order=1),
           alpha = guide_legend("signifance:", nrow=1, order=2))
    #theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}

#{{{ GO heatmap
gotag = 'GO_arabidopsis'
gotag = 'GO_uniprot.plants'
goname = gs$fun_ann %>% filter(ctag==gotag) %>% distinct(grp,note) %>%
    mutate(note = str_sub(note, 1, 60))
tx = ev_go %>% select(nid,enrich_grp) %>%
    unnest() %>%
    filter(net_size==5e4,ctag==gotag) %>%
    select(-net_size, -ctag) %>%
    filter(n>=50, pval<.01, fc>=2) %>%
    mutate(fc = log2(fc)) %>%
    inner_join(t_cfg, by = 'nid')

tx %>% select(-net_type) %>% inner_join(t_cfg, by = 'nid') %>%
    count(grp,net_type) %>% spread(net_type,n) %>%
    select(-ril) %>% inner_join(goname, by='grp') %>%
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
s_ncfg = t_cfg %>% mutate(nid=factor(nid, levels=s_nids)) %>% arrange(nid)
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

#{{{ # [obsolete] case study
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

#{{{ # [need rework] cross-network consistency
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

#{{{ natural variation evaluation
#{{{ read
nv = readRDS('~/projects/grn/data/06_deg/all.rds')
eopt = 'nv'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_nv = readRDS(fv)
#
fi = '~/projects/master.xlsx'
t_cfg_rn = read_xlsx(fi, 'barn') %>% select(yid, author)
#
nvs = nv %>% group_by(yid,cond,group1,group2) %>%
    summarise(popSize = n(), hitInPop = sum(DE!='non_DE'), bg.de = popSize, hitInPop) %>%
    ungroup() %>%
    select(yid,cond,group1,group2, popSize, hitInPop, bg.de)
#}}}

#{{{ nc01 case study
tp = ev_nv %>% filter(nid == 'nc01') %>%
    filter(!str_detect(group2, 'x')) %>%
    filter((cond=='root1' & group2 %in% c('Mo17','Oh43')) |
        (cond %in% c('leaf3_cold','leaf3_heat') & group2 %in% c('Mo17','PH207') |
        (cond %in% c('coleoptile','husk','spikelet')))) %>%
    group_by(score, yid,cond,group1,group2,reg.DE) %>%
    summarise(prop.de = sum(n[tgt.DE!='non_DE' & consis])/sum(n), n=sum(n)) %>%
    ungroup() %>% filter(n>=50) %>%
    inner_join(t_cfg_rn, by='yid') %>%
    mutate(pan = sprintf("%s %s:\n%s vs %s", str_to_title(author), cond, group1, group2))
tps = tp %>% distinct(yid,cond,group1,group2,pan) %>%
    inner_join(nvs, by=c('yid','cond','group1','group2'))
p1 = ggplot(tp) +
    geom_point(aes(x=score,y=prop.de,color=reg.DE, shape=reg.DE), size=2) +
    geom_line(aes(x=score,y=prop.de,color=reg.DE), size=.5) +
    geom_hline(data=tps, aes(yintercept=bg.de), linetype='dashed') +
    scale_x_continuous(name = 'TF-target rank') +
    scale_y_continuous(name='Proportion targets that are DE',limits=c(0,1), expand=expand_scale(mult=c(0,0))) +
    facet_wrap(~pan, ncol=3) +
    scale_color_aaas() +
    scale_shape_manual(values=0:5) +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
           legend.pos = 'top.center.out', legend.dir = 'h') +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.8))
fo = file.path(dirw, '11.nv.1.pdf')
ggsave(p1, file=fo, width=6, height=6)

tp = ev_nv %>% filter(nid == 'nc01') %>%
    filter(!str_detect(group2, 'x')) %>%
    group_by(score, yid,cond,group1,group2,reg.DE) %>%
    summarise(prop.de = sum(n[tgt.DE!='non_DE' & consis])/sum(n), n=sum(n)) %>%
    ungroup() %>% filter(n>=50) %>%
    inner_join(t_cfg_rn, by='yid') %>%
    mutate(pan = sprintf("%s %s:\n%s vs %s", str_to_title(author), cond, group1, group2))
tps = tp %>% distinct(cond,group1,group2,pan) %>%
    inner_join(nvs, by=c('yid','cond','group1','group2'))
p1 = ggplot(tp) +
    geom_point(aes(x=score,y=prop.de,color=reg.DE, shape=reg.DE), size=2) +
    geom_line(aes(x=score,y=prop.de,color=reg.DE), size=.5) +
    geom_hline(data=tps, aes(yintercept=bg.de), linetype='dashed') +
    scale_x_continuous(name = 'TF-target rank') +
    scale_y_continuous(name='Proportion targets that are DE',limits=c(0,1), expand=expand_scale(mult=c(0,0))) +
    facet_wrap(~pan, ncol=5) +
    scale_color_aaas() +
    scale_shape_manual(values=0:5) +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T,ytick=T,
           legend.pos = 'top.center.out', legend.dir = 'h') +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.8))
fo = file.path(dirw, '11.nv.1s.pdf')
ggsave(p1, file=fo, width=8, height=11)
#}}}

#{{{ DE enrichment
tp = ev_nv %>%
    filter(!str_detect(group2, 'x')) %>%
    filter(score >= 3, reg.DE == 'SPE') %>%
    group_by(nid,yid,cond,group1,group2,reg.DE) %>%
    summarise(sampleSize=sum(n), hitInSample = sum(n[tgt.DE!='non_DE' & consis])) %>%
    ungroup() %>% filter(sampleSize>=40) %>%
    inner_join(nvs, by=c('yid','cond','group1','group2')) %>%
    mutate(pval = pmap_dbl(list(hitInSample-1,hitInPop,popSize-hitInPop,sampleSize), phyper, lower.tail=F)) %>%
    mutate(pval = ifelse(pval==0, min(pval[pval!=0]), pval)) %>%
    mutate(logp = -log10(pval)) %>%
    mutate(lab = str_remove(sprintf("%.0f", logp), '^0+')) %>%
    mutate(lab = ifelse(pval<.05, lab, '')) %>%
    inner_join(t_cfg_rn, by='yid') %>%
    mutate(author = str_to_title(author)) %>%
    mutate(xlab0 = sprintf("%s %s: %s vs %s", author, cond, group1, group2)) %>%
    mutate(xlab1 = sprintf("%s: %s vs %s", cond, group1, group2)) %>%
    mutate(xlab2 = sprintf("%s vs %s", group1, group2)) %>%
    inner_join(t_cfg) %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
tpx = tp %>% distinct(author, cond, xlab0, xlab1, xlab2, group2) %>%
    arrange(xlab0) %>% mutate(x = 1:n()) %>%
    mutate(xcol=ifelse(group2=='Mo17',pal_uchicago()(3)[2], pal_uchicago()(3)[1]))
tpx1 = tpx %>% group_by(author,cond) %>% summarise(xs=min(x),xe=max(x),x=mean(x)) %>% ungroup()
tpx2 = tpx %>% group_by(author) %>% summarise(x=mean(x)) %>% ungroup()
tpy = tp %>% distinct(lgd, col) %>% arrange(lgd)
leg = 'Enrichment in Proportion Target DE (-log10 phyper):'
p1 = ggplot(tp) +
    geom_tile(aes(x=xlab0, y=lgd, fill = logp)) +
    geom_text(aes(x=xlab0, y=lgd, label = lab), size=2, lineheight=.8, color='black') +
    scale_x_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_y_discrete(name = '', expand = expand_scale(mult=c(0,0))) +
    scale_fill_gradientn(name=leg, colors = pal_gradient(reverse=T)) +
    #facet_grid(.~author, space='free_x', scales='free_x', switch='x') +
    #theme_classic(base_size = 14, base_family = 'mono') +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           margin = c(1, .1, .1, .1),
           xtitle=F, xtext=F, ytick=T, ytext=T) +
    #theme(panel.grid.minor.x = element_blank()) +
    #theme(panel.spacing.x = unit(0,"line")) +
    #theme(strip.placement='outside', strip.background.x=element_blank(), strip.text.x=element_text(angle = 90)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    #theme(axis.text.x = element_text(color = tpx$xcol, size=7, angle=90, vjust=1, hjust=1)) +
    theme(axis.text.y = element_text(color = tpy$col)) +
    expand_limits(y = -10, data = tibble(homeworld='pz')) +
    geom_rect(data=tibble(homeworld='pz'), ymin=-10,ymax=0,xmin=0,xmax=nrow(tpx)+2, fill='white') +
    geom_text(data=tpx, aes(x=x,label=xlab2), color=tpx$xcol, y=0, hjust=1, angle=90, size=2.5) +
    geom_segment(data=tpx1, aes(x=xs,xend=xe,y=-3.8,yend=-3.8), size=.5) +
    geom_text(data=tpx1, aes(x=x,label=cond), y=-4, hjust=1, vjust=1, angle=30, size=2.5) +
    geom_text(data=tpx2, aes(x=x,label=author), y=-8, hjust=.5, vjust=1, angle=30, size=3) +
    theme(panel.border = element_blank()) +
    guides(direction = 'horizontal', fill=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/11.nv.3.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 9)
#}}}
#}}}

#{{{ # [obsolete] evaluate using briggs data
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
#{{{ # [obsolete] eval using biomap data
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
#{{{ # [obsolete] check & plot high-conf TF-target pairs
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

#{{{ trans-eQTL hotspot
#{{{ read
eopt = 'go'
fv = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)
ev_go = readRDS(fv)
fun_ann_note = gs$fun_ann %>% distinct(ctag,grp,note)
t_gl = gcfg$gene %>% mutate(pos=(start+end)/2) %>% select(gid,chrom,pos)
#
qtags = c('li2013','liu2017','wang2018')
hs = tibble(qtag=qtags) %>%
    mutate(fi=sprintf('~/projects/genome/data2/%s/10.rds', qtag)) %>%
    mutate(data=map(fi, readRDS)) %>%
    mutate(data=map(data, 'hs')) %>%
    select(qtag, data) %>% unnest() %>%
    select(qtag,qid,qchrom,qpos,n.tgt)
#
fv = sprintf("%s/%s.50k.rds", dirr, gopt)
ev = readRDS(fv)
#
t_pick = tibble(reg.gid=c('Zm00001d026147','Zm00001d046405','Zm00001d053124','Zm00001d030028'),
    gname=c('R1', 'hb87','myb61','myc7'),
    grp = c('PWY-5125', 'CHLOROPHYLL-SYN', 'CALVIN-PWY','PWY-321'),
    #c('PWY-6787','GO:0010027','GO:0010027')
    path = c('anthocyanin biosynthesis',
      "3,8-divinyl-chlorophyllide a biosynthesis I (aerobic, light-dependent)",
      "Calvin-Benson-Bassham cycle",'cutin biosynthesis'))
t_pick = t_pick[c(1,2,3),]
#}}}

#{{{ prepare
nids_hc = c("nc01",'n13a','n15d')
nids_hc = c("nc01",'n18e','n99a','n18g','n13a','n15c','n15d','n18c')
tz0 = ev_go %>% filter(nid %in% nids_hc) %>%
    select(nid, enrich_reg) %>% unnest() %>%
    filter(ctag %in% qtags, n >= 10, net_size==5e4, pval < .05) %>%
    select(ctag, nid, reg.gid, n, fc, grp=max.grp, max.grp.size) %>%
    inner_join(t_cfg, by = 'nid')
tz = tz0 %>%
    group_by(ctag, reg.gid, grp) %>%
    summarise(n = n(), studies=list(study),
              max.grp.size = max(max.grp.size)) %>%
    ungroup() %>%
    filter(n>=2) %>% #arrange(ctag,grp,desc(fc)) %>%
    left_join(t_gl, by=c('reg.gid'='gid')) %>%
    rename(gchrom=chrom,gpos=pos) %>%
    left_join(hs, by=c('ctag'='qtag','grp'='qid')) %>%
    mutate(hit=ifelse(gchrom==qchrom & abs(gpos-qpos)<=5e7,'hit','non-hit'))
length(tz$hit)
length(unique(tz$reg.gid))
tz %>% filter(hit == 'hit') %>%
    select(ctag,reg.gid,n,max.grp.size,studies,gchrom,qchrom) %>% print(n=40)
length(tz$hit[tz$hit=='hit'])
length(unique(tz$reg.gid[tz$hit=='hit']))
#}}}

#{{{ plot
tx = gcfg$chrom
ti = tz
gpos = flattern_gcoord(ti %>% select(chrom=gchrom, pos=gpos), gcfg$chrom)
qpos = flattern_gcoord(ti %>% select(chrom=qchrom, pos=qpos), gcfg$chrom)
tp = ti %>% mutate(gpos=!!gpos, qpos=!!qpos)
p = ggplot(tp) +
    geom_point(aes(x=qpos, y=gpos, color=max.grp.size, shape=hit)) +
    geom_vline(xintercept = tx$start, alpha=.1) +
    geom_vline(xintercept = tx$end, alpha=.1) +
    geom_hline(yintercept = tx$start, alpha=.1) +
    geom_hline(yintercept = tx$end, alpha=.1) +
    scale_x_continuous(name='trans-eQTL hotsplot position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_y_continuous(name='co-regulating TF position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_color_viridis(name='N_targets:',option = "plasma") +
    scale_shape_manual(name='co-localization:',values=c(16,4),labels=c('yes','no')) +
    facet_wrap(~ctag,nrow=1) +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, legend.title=T,
           legend.pos='top.center.out', legend.dir='h', margin=c(1,.2,.2,.2)) +
    theme(legend.box = "horizontal")
    #theme(legend.background = element_rect(fill='white'))
fp = sprintf("%s/32.hs.pdf", dirw)
ggsave(p, filename = fp, width = 10, height = 4.2)

ti = tz %>% filter(hit=='hit')
gpos = flattern_gcoord(ti %>% select(chrom=gchrom, pos=gpos), gcfg$chrom)
qpos = flattern_gcoord(ti %>% select(chrom=qchrom, pos=qpos), gcfg$chrom)
tp = ti %>% mutate(gpos=!!gpos, qpos=!!qpos)
tps = t_pick %>% inner_join(tp, by = 'reg.gid')
p0 = ggplot(tp) +
    geom_point(aes(x=qpos, y=gpos, color=ctag, shape=hit)) +
    geom_text_repel(data=tps, aes(qpos,gpos,label=gname), nudge_x=-2e8, direction='y', segment.size=.3, size=2.5) +
    geom_vline(xintercept = tx$start, alpha=.1) +
    geom_vline(xintercept = tx$end, alpha=.1) +
    geom_hline(yintercept = tx$start, alpha=.1) +
    geom_hline(yintercept = tx$end, alpha=.1) +
    #geom_abline(intercept = 0, slope = 1, alpha=.1) +
    scale_x_continuous(name='trans-eQTL hotsplot position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_y_continuous(name='co-regulating TF position', breaks=tx$pos, labels=tx$chrom, expand=c(0,0)) +
    scale_color_d3(name = 'eQTL study') +
    scale_shape_manual(values=c(16,4)) +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, legend.title=T,
           legend.pos='bottom.right', legend.dir='v') +
    theme(legend.background = element_rect(fill='white')) +
    guides(shape = F)
#}}}

#{{{ save as table
reg.gids = tz %>% distinct(reg.gid) %>% pull(reg.gid)
tn0 = tz0 %>% select(nid,reg.gid,ctag,grp) %>% filter(reg.gid %in% reg.gids)
tn = ev %>% select(nid,tn) %>% unnest() %>%
    inner_join(tn0, by=c('nid','reg.gid')) %>%
    select(nid,reg.gid,tgt.gid, ctag,grp) %>%
    inner_join(gs$fun_ann, by=c('ctag','grp','tgt.gid'='gid')) %>%
    select(-note) %>%
    left_join(gcfg$gene[,c('gid','note2')], by=c('reg.gid'='gid')) %>%
    rename(reg.note=note2) %>%
    mutate(reg.note = str_replace(reg.note, ' *\\[.*\\]', ''))
#
t_fun = gs$fun_ann %>% distinct(ctag, grp, note)
tv0 = ev_go %>% select(nid, enrich_reg) %>% unnest() %>%
    inner_join(unique(tn0[c('nid','reg.gid')]), by=c('nid','reg.gid')) %>%
    filter(net_size==5e4) %>%
    filter(ctag %in% c("CornCyc",'GO_arabidopsis','GO_uniprot.plants')) %>%
    filter(fc > 5) %>% #arrange(reg.gid, nid) %>% print(n=80)
    select(nid, reg.gid,ctag,grp=max.grp, fc)
tv = tv0 %>%
    inner_join(t_fun, by=c('ctag','grp')) %>%
    group_by(reg.gid,ctag,note) %>%
    summarise(fc = max(fc)) %>% ungroup() %>%
    mutate(txt = sprintf("(%.0f) %s", fc, note)) %>%
    arrange(reg.gid, desc(fc)) %>%
    group_by(reg.gid) %>%
    summarise(txt = str_c(txt, collapse="; ")) %>%
    ungroup()
    #filter(reg.gid == 'Zm00001d046405') %>% print(n=50)
#
to = tz %>%
    group_by(reg.gid) %>%
    summarise(n_qtag = n(), studies=str_c(unlist(studies), collapse=', '),
              max.grp.size = max(max.grp.size),
              qtags = paste(ctag, collapse=',')) %>% ungroup() %>%
    filter(max.grp.size >= 1) %>%
    arrange(desc(n_qtag), desc(max.grp.size)) %>%
    left_join(gcfg$gene[,c('gid','note2')], by=c('reg.gid'='gid')) %>%
    rename(reg.note = note2) %>%
    left_join(tv, by='reg.gid') %>% replace_na(list(txt='')) %>%
    arrange(reg.gid)
to %>% print(n=50)

fo = file.path(dirw, '33.hs.tsv')
write_tsv(to, fo)
#}}}

#{{{ # [obsolete] individual loci test
ztag = 'GO_uniprot.plants'
ztag = 'CornCyc'
tfid = 'Zm00001d026147'
tgid = 'Zm00001d017077'
tgid = 'Zm00001d026141'
tz2 = ev_go %>% filter(nid %in% nids_hc) %>%
    select(nid, enrich_reg) %>% unnest() %>%
    filter(ctag == ztag, n >= 5, net_size==5e4, pval < .05) %>%
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
#}}}

#{{{ pathway plot
#{{{ read trans-eQTL - requires tv0 and tz
tn0 = ev %>% #filter(nid %in% nids_hc) %>%
    select(nid,tn) %>% unnest() %>%
    select(nid,reg.gid,tgt.gid)
tv0 %>% filter(reg.gid %in% t_pick$reg.gid) %>% arrange(reg.gid,desc(fc)) %>% print(n=40)
#x = to %>% filter(reg.gid %in% t_pick$reg.gid) %>% print(width=Inf)
#x$txt
#}}}
#{{{ read
require(ggnetwork)
require(network)
fi = '~/projects/genome/data/Zmays_B73/61_functional/07.corncyc.rds'
cc = readRDS(fi)
prepare_nodes <- function(n0) { #from_id, to_id, slab, elab, stype, etype
    #{{{
    require(Rgraphviz)
    zn1 = ftM2graphNEL(as.matrix(n0[,c('from_id','to_id')]), edgemode="directed")
    attrs <- list(graph=list(rankdir="LR"))
    zn2 = layoutGraph(zn1, attrs=attrs)
    zn = nodeRenderInfo(zn2) %>% as_tibble() %>%
        select(x=nodeX, y=nodeY) %>% mutate(node=nodes(zn2))
    #
    v1 = n0 %>% select(node=from_id, type=stype, lab=slab)
    v2 = n0 %>% select(node=to_id, type=etype, lab=elab)
    v = rbind(v1,v2) %>% distinct(node, type, lab) %>%
        inner_join(zn, by='node')
    v
    #}}}
}
#}}}

#i = 2
plot_pathway <- function(i) {
#{{{
grp = t_pick$grp[i]; reg.gid = t_pick$reg.gid[i]; path=t_pick$path[i]
q_ctag = tz$ctag[tz$reg.gid==reg.gid]; q_grp = tz$grp[tz$reg.gid==reg.gid]
#{{{ get GRN and eQTL support
tps = gs$fun_ann %>% select(grp, tgt.gid=gid) %>%
    filter(grp == !!grp) %>% mutate(reg.gid = !!reg.gid) %>%
    left_join(tn0, by=c('reg.gid','tgt.gid')) %>%
    group_by(tgt.gid) %>%
    summarise(n_net = sum(!is.na(nid)), nids = str_c(nid, collapse=',')) %>%
    ungroup() %>% mutate(grn = n_net > 0) %>%
    left_join(gcfg$gene[,c('gid','note2')], by=c('tgt.gid'='gid')) %>%
    rename(note=note2) %>%
    mutate(note = str_replace(note, ' *\\[.*\\]', ''))
tq = gs$fun_ann %>% filter(ctag==q_ctag, grp==q_grp)
tps = tps %>% mutate(eqtl = tgt.gid %in% tq$gid)
tps %>% count(grn, eqtl)
gnote = tps$note; names(gnote) = tps$tgt.gid
#}}}
#
concat_gene_desc <- function(gid, gnote)
    ifelse(gid %in% names(gnote), str_c(gid, gnote[gid], sep=' '), gid)
exc_ids = c('2-oxoglutarate','UDP-alpha-D-glucose','succinate',
            "hydrogen peroxide",'phosphate','Mg2+',
            "AMP","diphosphate",'beta-D-fructose 1,6-bisphosphate')
n0 = cc %>% filter(pathway==path) %>% select(pathway,net) %>% unnest() %>%
    filter(!snode %in% exc_ids, !enode %in% exc_ids) %>%
    distinct(pathway,snode,enode,stype,etype) %>%
    mutate(slab = map_chr(snode, concat_gene_desc, gnote=gnote)) %>%
    mutate(elab = map_chr(enode, concat_gene_desc, gnote=gnote)) %>%
    mutate(slab = str_wrap(slab, width=25)) %>%
    mutate(elab = str_wrap(elab, width=25)) %>%
    rename(from_id=snode,to_id=enode) %>%
    mutate(type = str_c(stype,etype,sep='')) %>%
    select(from_id, to_id, everything())
#
n1 = network(n0, matrix.type='edgelist')
vlevels = network.vertex.names(n1)
v = prepare_nodes(n0) %>%
    left_join(tps[,c('tgt.gid','grn')], by=c('node'='tgt.gid')) %>%
    mutate(type = ifelse(type=='g' & !is.na(grn) & grn, 'g2', type)) %>%
    mutate(node=factor(node,levels=vlevels)) %>% arrange(node)
n1 %v% 'type' = v$type
n1 %v% 'lab' = v$lab
n2 = ggnetwork(n1, layout=as.matrix(v[,c('x','y')]))
cols3 = c(brewer.pal(8,'Paired')[5:6],pal_lancet()(8)[8])
cols3 = c(pal_npg()(4)[2:1],pal_uchicago()(4)[1])
leg.pos = ifelse(i==1, 'top.right', 'none')
brks=c('s','g','g2')
labs=c('substrate / product', 'gene / enzyme', 'gene / enzyme [GRN supported]')
p = ggplot(n2, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_edges(aes(),size=.3,arrow=arrow(angle=15,length=unit(5,'pt'),type='open')) +
    geom_nodes(aes(color=type,shape=type), size=4) +
    geom_nodetext_repel(aes(label=lab), size=2) +
    scale_color_manual(values=cols3,breaks=brks,labels=labs) +
    scale_shape_manual(values=c(16,16,2),breaks=brks,labels=labs) +
    otheme(legend.pos=leg.pos, margin=c(.1,.1,.1,.1)) +
    theme(panel.border = element_blank()) +
    guides(linetype=F)
p
#}}}
}

p1 = plot_pathway(1)
p2 = plot_pathway(2)
p3 = plot_pathway(3)
fo= file.path(dirw, '35.1.pdf')
ggarrange(p0, p1, p2, p3,
    nrow = 2, ncol = 2, labels = LETTERS[1:4], heights = c(2,2)) %>%
    ggexport(filename = fo, width=10, height=10)
#}}}



