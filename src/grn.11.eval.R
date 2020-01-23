source("functions.R")
gopts = c("rf",'et','xgb')
colmap = pal_aaas()(3)
names(colmap) = gopts
dirw = file.path(dird, '13_eval')

#{{{ network clustering
ti = tibble(gopt = gopts) %>%
    mutate(fi = sprintf("%s/%s.50k.rds", dirr, gopt)) %>%
    mutate(data = map(fi, readRDS)) %>%
    unnest(cols=data) %>% select(-fi) %>%
    filter(nid %in% t_cfg$nid) %>%
    mutate(nid = sprintf("%s_%s", nid, str_to_upper(gopt))) %>%
    mutate(lgd = sprintf("%s_%s", lgd, str_to_upper(gopt)))

tu = ti %>% select(nid, tn) %>% unnest() %>%
    select(nid, reg.gid, tgt.gid, score) %>%
    group_by(nid) %>% slice(1:1e5) %>% ungroup() %>%
    mutate(reg.tgt = str_c(reg.gid, tgt.gid, sep = '-')) %>%
    select(-reg.gid, -tgt.gid)
reg.tgts = tu %>% count(reg.tgt) %>% filter(n>=2) %>% pull(reg.tgt)
length(reg.tgts)
tu = tu %>% filter(reg.tgt %in% reg.tgts)
#
tuw = tu %>% spread(reg.tgt, score)
mat = as.matrix(tuw[,-1])
rownames(mat) = tuw$nid
mat[is.na(mat)] = 0
dist_mat = dist(mat, method = 'binary')

#{{{ only hc-tree
tp = as.matrix(dist_mat)
colnames(tp) = tuw$nid
tp = as_tibble(tp) %>% mutate(aid = tuw$nid) %>%
    gather(bid, dis, -aid) %>%
    mutate(dis = ifelse(dis <= 0, 0, dis)) %>%
    mutate(dis = ifelse(aid==bid, NA, dis)) %>%
    rename(dist = dis)
#
require(cluster)
require(ape)
require(ggtree)
hc = hclust(as.dist(dist_mat), method = "ward.D")
oidx = hc$order
tree = as.phylo(hc)
#
tpt = ti %>% select(taxa=nid, everything())
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=lgd, col=gopt), size=2.5) +
    scale_x_continuous(expand=expand_scale(mult=c(0,.3))) +
    scale_color_aaas() +
    theme(plot.margin = margin(-1,.5,-1,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '02.gopt.1.hc.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 11)
#}}}

#{{{ hc + heatmap
th0 = ti %>% select(nid, lgd, gopt)
#dis = daisy(tuw[,-1], metric = 'gower')
tp = as.matrix(dist_mat)
colnames(tp) = tuw$nid
tp = as_tibble(tp) %>% mutate(aid = tuw$nid) %>%
    gather(bid, dis, -aid) %>%
    #mutate(dis = ifelse(dis <= .6, .6, dis)) %>%
    mutate(dis = ifelse(aid==bid, NA, dis)) %>%
    rename(dist = dis) %>%
    inner_join(th0, by=c('aid'='nid')) %>% rename(algd=lgd, acol=gopt) %>%
    inner_join(th0, by=c('bid'='nid')) %>% rename(blgd=lgd, bcol=gopt) %>%
    mutate(acol = colmap[acol], bcol = colmap[bcol])
#
p = heatmap_hc(tp, top=0, bottom=4.6, r.top=1.2, text.size=1,ratio=4)
fo = file.path(dirw, '02.gopt.2.heat.pdf')
p %>% ggexport(filename = fo, width = 14, height = 12)
#}}}

#{{{ tSNE
require(Rtsne)
tt = tu %>% spread(nid, score) %>% replace(., is.na(.), 0)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=9,
              pca = T, max_iter = 2000)

tp = as_tibble(tsne$Y) %>%
    add_column(nid = colnames(tt)[-1]) %>%
    inner_join(ti[,c('gopt','nid','lgd','col','net_type')], by = 'nid') %>%
    mutate(lgd0 = str_replace(lgd,'_[A-Z]{2,3}$',''))
tps = tp %>% mutate(lgd0 = ifelse(gopt=='rf', lgd0, ''))
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1, y=V2)) +
    geom_text_repel(data=tps,aes(x=V1,y=V2,label=lgd0), size=3) +
    geom_mark_ellipse(aes(fill=lgd0),
        expand=unit(2,'mm'), alpha=0, size = .2,
        con.type='none',label.fontsize=7,label.minwidth=unit(0,'mm'),
        label.buffer=unit(0,'mm'),label.margin = margin(0,0,0,0,"mm")) +
    geom_point(aes( color=gopt, shape=gopt), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_color_aaas() +
    scale_shape(solid=F) +
    otheme(legend.pos='top.left', legend.dir='v', legend.title=F,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(fill=F)
fp = file.path(dirw, "02.tsne.pdf")
ggsave(p_tsne, filename = fp, width=10, height=10)
#}}}
#}}}

#{{{ read eval
ev_ko = read_eval_ko()
ev_bs = read_eval_bs()
ev_go = read_eval_go()

evk = ev_ko %>% filter(net_size==1e7) %>%  mutate(xlab=ctag)
evc = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^REF'))
evd = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^(Galli2018)|(Ricci2019)$'))
evdu = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^(Galli2018)|(Ricci2019)_umr$'))
evde = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^(Galli2018)|(Ricci2019)_acrE$'))
evdl = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^(Galli2018)|(Ricci2019)_acrL$'))
evbc = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^cisbp'))
evbp = ev_bs %>% filter(net_size==1e7) %>% filter(str_detect(ctag, '^plantregmap'))
#}}}

#{{{ plot chip / dap / ko / bs
ev = evc; eopt = 'cp'; wid=6; hei=10
ev = evk; eopt = 'ko'; wid=12; hei=10
ev = evd; eopt = 'dp'; wid=15; hei=10

tp1 = ev %>% rename(score=score1, lab=lab1)
tp2 = ev %>% rename(score=score2, lab=lab2)
tp3 = ev %>% rename(score=score3, lab=lab3)
tp4 = ev %>% rename(score=score4, lab=lab4)
p1 = plot_tile(tp1, t_cfg, lgd.opt=1, faceting=T)
p2 = plot_tile(tp2, t_cfg, lgd.opt=2, faceting=T)
p3 = plot_tile(tp3, t_cfg, lgd.opt=3, faceting=T)
p4 = plot_tile(tp4, t_cfg, lgd.opt=4, faceting=T)
fo = sprintf('%s/11.%s.pdf', dirw, eopt)
ggarrange(p1, p2,
    nrow = 2, ncol = 1, labels = LETTERS[1:2], heights = c(2,2)) %>%
    ggexport(filename = fo, width = wid, height = hei)
fo = sprintf('%s/11.%s2.pdf', dirw, eopt)
ggarrange(p3, p4,
    nrow = 2, ncol = 1, labels = LETTERS[1:2], heights = c(2,2)) %>%
    ggexport(filename = fo, width = wid, height = hei)
#}}}

#{{{ GO evaluation
#{{{ enrichment [heatmap]
ctags = c("GO_HC","GO_arabidopsis_P","GO_uniprot.plants_P","CornCyc")
tp = ev_go %>% select(gopt, nid, enrich) %>% unnest(enrich) %>%
    mutate(gopt = str_to_upper(gopt)) %>%
    filter(ctag %in% ctags, score == 10) %>%
    group_by(gopt, ctag) %>%
    #mutate(fcn = scale(fc)) %>% ungroup() %>%
    mutate(fcn = fc) %>% ungroup() %>%
    mutate(lab = sprintf("%.01f", fc)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
tps = tp %>% distinct(lgd, col) %>% arrange(lgd)
swit = max(tp$fcn) / 2
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=fcn)) +
    geom_tile() +
    geom_text(aes(x=ctag, y=lgd, label=lab, color=fc>swit), hjust=.5, size=2.5) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name = 'Fold Enrichment', direction=-1) +
    scale_fill_gradientn(name='Fold Enrichment', colors=cols100v) +
    scale_color_manual(values=c('black','white')) +
    facet_grid(.~gopt) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=tps$col)) +
    guides(color = F)
fp = sprintf('%s/15.go.heat.pdf', dirw)
ggsave(p1, file = fp, width = 6, height = 6)
#}}}

#{{{ # enrichment [obsolete]
ctags = c("GO_HC", "GO_arabidopsis","CornCyc")
fo = file.path(dirw, "15.go.pdf")
tp1 = evg %>% select(nid, gopt, lgd, enrich) %>% unnest() %>% filter(score==10)
tp2 = evg %>% select(nid, gopt, enrich_grp) %>% unnest() %>%
    filter(score == 10) %>%
    filter(n >= 10) %>%
    group_by(nid, gopt, score, ctag) %>%
    summarise(n_grp=length(grp), n_grp_sig=sum(pval<.05)) %>%
    ungroup() %>%
    mutate(sigtxt = str_c(n_grp_sig,n_grp,sep='/')) %>%
    select(nid,gopt,net_size,ctag,sigtxt)
tp = tp1 %>% left_join(tp2, by=c('nid','gopt','net_size','ctag')) %>%
    mutate(sig = ifelse(pval<.05, 1, 0)) %>%
    filter(ctag %in% ctags) %>%
    mutate(fc = log2(fc)) %>%
    mutate(gopt = str_to_upper(gopt)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
tpi = tp %>% filter(sig == 0)
tps = tp %>% distinct(lgd, col) %>% arrange(lgd)
#
p1 = ggplot(tp, aes(lgd, fc)) +
    geom_point(aes(color=gopt, shape=gopt), size=2) +
    geom_point(data=tpi, aes(x=lgd,y=fc), shape=4, size=2) +
    geom_text_repel(aes(label=sigtxt), size=2) +
    #geom_hline(yintercept = 1, alpha= .5, linetype='dotted') +
    scale_x_discrete(expand = expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name = 'log2 Fold Enrichment', expand = c(.05,0)) +
    coord_flip() +
    facet_grid(.~ctag, scale='free') +
    scale_color_npg() +
    scale_shape_manual(values=c(0:5)) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           strip.size = 8,
           xtitle=T, xtext=T, ytext=T, ygrid=T, xtick=T, ytick=T) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col))) +
    theme(legend.box = "horizontal") +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.3)) +
    guides(color = guide_legend("inference approach:", nrow=1, order=1),
           shape = guide_legend("inference approach:", nrow=1, order=1))
    #theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}
#}}}





#{{{ # [obsolete] eval TF KO / TFBS AUROC
ev = tibble(gopt=gopts, fv=sprintf("%s/%s.1m.rds", dirr, gopts)) %>%
    mutate(data = map(fv, readRDS)) %>% select(-fv) %>% unnest()

#{{{ read ko
nrow_positive <- function(ti) sum(ti$response == 1)
ko = read_ko() %>% filter(gene_alias != 'P1', gene_alias != 'fl3') %>%
    mutate(kid=1:n()) %>%
    select(kid,gene_id,gene_alias,Tissue,ds) %>%
    unnest() %>%
    mutate(response=ifelse(padj < .01, 1, 0)) %>%
    rename(reg.gid = gene_id, tgt.gid = gid) %>%
    group_by(kid, gene_alias, Tissue) %>% nest() %>% rename(res = data) %>%
    mutate(n_tot = map_dbl(res, nrow)) %>%
    mutate(n_de = map_dbl(res, nrow_positive)) %>%
    mutate(prop_de=n_de/n_tot) %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags_ko = ko$ctag
rids_ko = ko %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
ko_u = ko %>% unnest()
#}}}
#{{{ read tf (ko direct)
tf0 = read_ko_direct() %>% select(gene_alias,tissue,ctag=lab,data)
tfs = tf0 %>% select(gene_alias, tissue, ctag)
tmp = ko %>% rename(tissue=Tissue) %>%
    mutate(gene_alias = str_to_upper(gene_alias)) %>%
    mutate(gene_alias=ifelse(gene_alias=='BZIP22','bZIP22',gene_alias)) %>%
    mutate(tissue=ifelse(tissue=='ear_1mm', 'ear',tissue)) %>%
    filter(gene_alias %in% tf0$gene_alias) %>%
    filter(tissue != 'SAM', tissue != 'ear_2mm') %>%
    unnest() %>%
    select(gene_alias,tissue, reg.gid, tgt.gid)
tf = tf0 %>% select(-ctag) %>% unnest() %>%
    mutate(response=1) %>%
    right_join(tmp, by=c('gene_alias','tissue','reg.gid','tgt.gid')) %>%
    replace_na(list(response=0)) %>%
    inner_join(tfs, by=c('gene_alias','tissue')) %>%
    group_by(gene_alias,tissue,ctag) %>% nest() %>%
    rename(res=data)
ctags_tf = tf$ctag
rids_tf = tf %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
tf_u = tf %>% unnest()
#}}}
#{{{ read y1h
rids = gs$y1h %>% filter(reg.gid %in% gs$tf_ids) %>% distinct(reg.gid) %>% pull(reg.gid)
tmp = gs$y1h %>% distinct(tgt.gid) %>% crossing(reg.gid = rids)
nrow_positive <- function(ti) sum(ti$response == 1)
yh = gs$y1h %>% mutate(response = 1) %>%
    right_join(tmp, by=c('reg.gid','tgt.gid')) %>%
    mutate(ctag = 'Y1H') %>%
    replace_na(list(response=0)) %>%
    group_by(ctag) %>% nest() %>%
    rename(res = data) %>% arrange(ctag) %>%
    mutate(n = map_dbl(res, nrow_positive)) %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = sprintf("%s [%s]", ctag, number(n))) %>%
    mutate(ctag = factor(ctag, levels=ctag))
ctags_yh = yh$ctag
rids_yh = yh %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
yh_u = yh %>% unnest()

#}}}
#{{{ read tfbs
gids = gcfg$gene %>% filter(ttype=='mRNA') %>% pull(gid)
tmp = gs$tfbs %>% distinct(ctag, reg.gid) %>% crossing(tgt.gid = gids)
nrow_positive <- function(ti) sum(ti$response == 1)
bs = gs$tfbs %>% mutate(response = 1) %>%
    right_join(tmp, by=c('ctag','reg.gid','tgt.gid')) %>%
    replace_na(list(response=0)) %>%
    group_by(ctag) %>% nest() %>%
    rename(res = data) %>% arrange(ctag) %>%
    mutate(n = map_dbl(res, nrow_positive)) %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = sprintf("%s [%s]", ctag, number(n))) %>%
    mutate(ctag = factor(ctag, levels=ctag))
ctags_bs = bs$ctag
rids_bs = bs %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
bs_u = bs %>% unnest()
#}}}

#{{{ functions
complete_tn <- function(tn, tids, rids) {
    #{{{
    if(sum(rids %in% tn$reg.gid)==0)
        tibble()
    else {
        rids = rids[rids %in% tn$reg.gid]
        crossing(reg.gid = rids, tgt.gid = tids) %>% as_tibble() %>%
            left_join(tn, by=c('reg.gid','tgt.gid')) %>%
            replace_na(list(score=0))
    }
    #}}}
}
complete_tn2 <- function(tn, rids, tids) {
    #{{{
    if(sum(tids %in% tn$tgt.gid)==0)
        tibble()
    else {
        tids = tids[tids %in% tn$tgt.gid]
        crossing(reg.gid = rids, tgt.gid = tids) %>% as_tibble() %>%
            left_join(tn, by=c('reg.gid','tgt.gid')) %>%
            replace_na(list(score=0))
    }
    #}}}
}
filter_tn <- function(tn, rids) {
    #{{{
    if(sum(rids %in% tn$reg.gid)==0)
        tibble()
    else
        tn %>% filter(reg.gid %in% rids)
    #}}}
}
get_auroc <- function(ti, fpr=1) {
    #{{{
    if(max(ti$response) == 0 || max(ti$score) == 0)
        NA
    else if (min(ti$response) == 1)
        NA
    else
        roc(ti$response, ti$score, partial.auc=c(1,1-fpr), levels=c(0,1))$auc
    #}}}
}
get_spc_pval <- function(ti)
    cor.test(ti$score, ti$padj, method='kendall')$p.value
get_wil_pval <- function(ti) {
    #{{{
    scores1 = ti %>% filter(response==1) %>% pull(score)
    scores2 = ti %>% filter(response==0) %>% pull(score)
    if(length(scores1) > 0 & length(scores2) > 0)
        wilcox.test(scores1, scores2, alternative='greater')$p.value
    else
        NA
    #}}}
}
#}}}

rids = tf$gid[13]
rids = gs$bs %>% filter(str_detect(ctag, 'Ricci')) %>% distinct(reg.gid) %>% pull(reg.gid)

x1 = gs$bs %>% filter(ctag == 'Ricci2019', reg.gid %in% rids) %>%
    mutate(tf = reg.gid) %>%
    group_by(ctag, tf) %>% nest()
tne = ev %>% select(gopt,nid,rids,tids,tn) %>% filter(nid=='nc01', gopt=='rf')
rids = unlist(tne$rids); tids = unlist(tne$tids); tn1 = tne$tn[[1]]

ev2 = x1 %>%
    mutate(enc = map(data, eval_bs1, tn=tn1, rids=rids, tids=tids)) %>%
    mutate(auroc=map_dbl(enc, 'auroc')) %>%
    mutate(pval=map_dbl(enc, 'pval')) %>%
    select(ctag, tf, auroc, pval)

eval_bs1 <- function(t_bs, tn, rids, tids, fpr=.1, net_size=1e7) {
    #{{{
    tb0 = t_bs %>% mutate(response = 1)
    rids0 = rids[rids %in% tb0$reg.gid]
    tt = complete_tn2(tn, rids0, tids) %>%
        filter(reg.gid %in% rids0) %>%
        filter(reg.gid != tgt.gid) %>%
        mutate(score = as.numeric(score)) %>%
        mutate(score = ifelse(is.na(score), 0, score)) %>%
        arrange(desc(score)) %>%
        filter(row_number() <= net_size)
    if(max(tt$score) == 0) tt$score[1] = 0.1
    to = tt %>%
        left_join(tb0, by = c('reg.gid','tgt.gid')) %>%
        replace_na(list(response=0))
    # auroc at FPR
    auroc = NA
    if(max(to$response) != 0 && max(to$score) != 0 && min(to$response) != 1)
        auroc = roc(to$response, to$score, partial.auc=c(1,1-fpr), levels=c(0,1))$auc
    # wilcox p-value
    pval = NA
    scores1 = to %>% filter(response==1) %>% pull(score)
    scores2 = to %>% filter(response==0) %>% pull(score)
    if(length(scores1) > 0 & length(scores2) > 0)
        pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
    list(auroc=auroc, pval=pval)
    #}}}
}

#{{{ auroc + pval
tp_tf = ev %>% select(gopt,nid,tids,tn) %>%
    mutate(res = map2(tn, tids, complete_tn, rids = rids_tf)) %>%
    select(gopt, nid, res) %>% unnest() %>%
    select(gopt, nid, reg.gid, tgt.gid, score) %>%
    inner_join(tf_u, by=c('reg.gid','tgt.gid')) %>%
    group_by(gopt, nid, gene_alias, tissue, ctag) %>%
    nest() %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval)) %>% select(-data)

tp_yh = ev %>% select(gopt,nid,rids,tn) %>%
    mutate(res = map2(tn, rids, complete_tn2, tids = unique(gs$y1h$tgt.gid))) %>%
    select(gopt, nid, res) %>% unnest() %>%
    select(gopt, nid, reg.gid, tgt.gid, score) %>%
    inner_join(yh_u, by=c('reg.gid','tgt.gid')) %>%
    group_by(gopt, nid, ctag) %>%
    nest() %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval)) %>% select(-data)

tp_ko = ev %>% select(gopt,nid,tids,tn) %>%
    mutate(res = map2(tn, tids, complete_tn, rids = rids_ko)) %>%
    select(gopt, nid, res) %>% unnest() %>%
    select(gopt, nid, reg.gid, tgt.gid, score) %>%
    inner_join(ko_u, by=c('reg.gid','tgt.gid')) %>% rename(tissue=Tissue) %>%
    group_by(gopt, nid, kid, gene_alias, tissue, ctag) %>%
    nest() %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval)) %>% select(-data)
#
tp_bs = ev %>% select(gopt,nid,tids,tn) %>%
    mutate(res = map(tn, filter_tn, rids = rids_bs)) %>%
    select(gopt, nid, res) %>% unnest() %>%
    select(gopt, nid, reg.gid, tgt.gid, score) %>%
    inner_join(bs_u, by=c('reg.gid','tgt.gid')) %>%
    group_by(gopt, nid, ctag) %>%
    nest() %>%
    mutate(auroc = map_dbl(data, get_auroc, fpr=.1)) %>%
    mutate(spc.pval = map_dbl(data, get_wil_pval)) %>% select(-data)
#
res = list(tf=tp_tf, ko=tp_ko, bs=tp_bs)
fo = file.path(dird, '13_eval', '05.tf.ko.bs.rds')
saveRDS(res, file=fo)
#}}}
#}}}
#{{{ # [obsolete] read meta info
#{{{ read ko
nrow_positive <- function(ti) sum(ti$response == 1)
ko = read_ko() %>% filter(gene_alias != 'P1', gene_alias != 'fl3') %>%
    mutate(kid=1:n()) %>%
    select(kid,gene_id,gene_alias,Tissue,ds) %>%
    unnest() %>%
    mutate(response=ifelse(padj < .01, 1, 0)) %>%
    rename(reg.gid = gene_id, tgt.gid = gid) %>%
    group_by(kid, gene_alias, Tissue) %>% nest() %>% rename(res = data) %>%
    mutate(n_tot = map_dbl(res, nrow)) %>%
    mutate(n_de = map_dbl(res, nrow_positive)) %>%
    mutate(prop_de=n_de/n_tot) %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags_ko = ko$ctag
rids_ko = ko %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
ko_u = ko %>% unnest()
#}}}
#{{{ read tf (ko direct)
tf0 = read_ko_direct() %>% select(gene_alias,tissue,ctag=lab,data)
tfs = tf0 %>% select(gene_alias, tissue, ctag)
tmp = ko %>% rename(tissue=Tissue) %>%
    mutate(gene_alias = str_to_upper(gene_alias)) %>%
    mutate(gene_alias=ifelse(gene_alias=='BZIP22','bZIP22',gene_alias)) %>%
    mutate(tissue=ifelse(tissue=='ear_1mm', 'ear',tissue)) %>%
    filter(gene_alias %in% tf0$gene_alias) %>%
    filter(tissue != 'SAM', tissue != 'ear_2mm') %>%
    unnest() %>%
    select(gene_alias,tissue, reg.gid, tgt.gid)
tf = tf0 %>% select(-ctag) %>% unnest() %>%
    mutate(response=1) %>%
    right_join(tmp, by=c('gene_alias','tissue','reg.gid','tgt.gid')) %>%
    replace_na(list(response=0)) %>%
    inner_join(tfs, by=c('gene_alias','tissue')) %>%
    group_by(gene_alias,tissue,ctag) %>% nest() %>%
    rename(res=data)
ctags_tf = tf$ctag
rids_tf = tf %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
tf_u = tf %>% unnest()
#}}}
#{{{ read y1h
rids = gs$y1h %>% filter(reg.gid %in% gs$tf_ids) %>% distinct(reg.gid) %>% pull(reg.gid)
tmp = gs$y1h %>% distinct(tgt.gid) %>% crossing(reg.gid = rids)
nrow_positive <- function(ti) sum(ti$response == 1)
yh = gs$y1h %>% mutate(response = 1) %>%
    right_join(tmp, by=c('reg.gid','tgt.gid')) %>%
    mutate(ctag = 'Y1H') %>%
    replace_na(list(response=0)) %>%
    group_by(ctag) %>% nest() %>%
    rename(res = data) %>% arrange(ctag) %>%
    mutate(n = map_dbl(res, nrow_positive)) %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = sprintf("%s [%s]", ctag, number(n))) %>%
    mutate(ctag = factor(ctag, levels=ctag))
ctags_yh = yh$ctag
rids_yh = yh %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
yh_u = yh %>% unnest()

#}}}
#{{{ read tfbs
gids = gcfg$gene %>% filter(ttype=='mRNA') %>% pull(gid)
tmp = gs$tfbs %>% distinct(ctag, reg.gid) %>% crossing(tgt.gid = gids)
nrow_positive <- function(ti) sum(ti$response == 1)
bs = gs$tfbs %>% mutate(response = 1) %>%
    right_join(tmp, by=c('ctag','reg.gid','tgt.gid')) %>%
    replace_na(list(response=0)) %>%
    group_by(ctag) %>% nest() %>%
    rename(res = data) %>% arrange(ctag) %>%
    mutate(n = map_dbl(res, nrow_positive)) %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = sprintf("%s [%s]", ctag, number(n))) %>%
    mutate(ctag = factor(ctag, levels=ctag))
ctags_bs = bs$ctag
rids_bs = bs %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
bs_u = bs %>% unnest()
#}}}
#}}}


#{{{ # test bx gene clusters
tf = read_tf_info()
tfbs_regs = gs$bs
gids = gcfg$gene %>% filter(ttype=='mRNA') %>% pull(gid)
tmp = tfbs_regs %>% distinct(ctag, reg.gid) %>% crossing(tgt.gid=gids)
tf %>% select(tf, gid) %>% print(n=40)

ev = tibble(gopt=gopts, fv=sprintf("%s/%s.100k.rds", dirr, gopts)) %>%
    mutate(data = map(fv, readRDS)) %>% select(-fv) %>% unnest()

gids_bx = tg %>% filter(pathway=='dimboa', role=='tgt') %>% pull(gid)

tfbs_regs %>% filter(ctag=='cisbp',reg.gid==tg$gid[5])
tfbs_regs %>% filter(ctag=='cisbp',reg.gid==tg$gid[18], tgt.gid %in% gids_bx)

nrow_positive <- function(ti) sum(ti$response == 1)
bs = tfbs_regs %>% rename(response = score) %>%
    right_join(tmp, by=c('ctag','reg.gid','tgt.gid')) %>%
    replace_na(list(response=0)) %>%
    mutate(ctag=str_c(ctag, reg.gid, sep=' ')) %>%
    group_by(ctag) %>% nest() %>%
    rename(res = data) %>% arrange(ctag) %>%
    mutate(n = map_dbl(res, nrow_positive)) %>%
    filter(n >= 10) %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = sprintf("%s [%s]", ctag, number(n))) %>%
    mutate(ctag = factor(ctag, levels=ctag))
ctags_bs = bs$ctag
rids_bs = bs %>% unnest() %>% distinct(reg.gid) %>% pull(reg.gid)
bs_u = bs %>% unnest()
#}}}




