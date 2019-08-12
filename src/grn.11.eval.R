source("functions.R")
gopts = c("rf",'et','xgb')
colmap = pal_aaas()(3)
names(colmap) = gopts
dirw = file.path(dird, '13_eval')

#{{{ network clustering
ti = tibble(gopt = gopts) %>%
    mutate(fi = sprintf("%s/%s.50k.rds", dirr, gopt)) %>%
    mutate(data = map(fi, readRDS)) %>%
    unnest() %>% select(-fi) %>%
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
    theme(plot.margin = margin(-1,.5,-1,.5, 'lines'))
#
fo = file.path(dirw, '02.gopt.1.hc.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 10)
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
p = heatmap_hc(tp, top=0, bottom=4.8, r.top=1, text.size=1.5,ratio=4)
fo = file.path(dirw, '02.gopt.2.heat.pdf')
p %>% ggexport(filename = fo, width = 15, height = 12)
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

#{{{ eval TF KO / TFBS AUROC [takes long]
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
#
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
#
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

ev = read_tf_ko_bs()
#{{{ plot tf / ko / bs
eopt = 'bs'; wid=5; hei=10
eopt = 'tf'; wid=6; hei=10
eopt = 'ko'; wid=12; hei=10
tp1 = ev %>% filter(eopt==!!eopt, key=='auroc')
tp2 = ev %>% filter(eopt==!!eopt, key=='spc.pval')
p1 = plot_tile(tp1, lgd.opt=1, faceting=T)
p2 = plot_tile(tp2, lgd.opt=2, faceting=T)
fo = sprintf('%s/11.%s.pdf', dirw, eopt)
ggarrange(p1, p2,
    nrow = 2, ncol = 1, labels = LETTERS[1:2], heights = c(2,2)) %>%
    ggexport(filename = fo, width = wid, height = hei)
#}}}


#{{{ GO evaluation
evg = tibble(gopt = gopts) %>%
    mutate(fi = sprintf("%s/%s.go.rds", dirr, gopt)) %>%
    mutate(data = map(fi, readRDS)) %>% select(-fi) %>% unnest()

#{{{ enrichment [heatmap]
ctags = c("GO_HC","GO_arabidopsis","GO_uniprot.plants","CornCyc")
tp = evg %>% select(gopt, nid, enrich) %>% unnest() %>%
    mutate(gopt = str_to_upper(gopt)) %>%
    filter(ctag %in% ctags, score == 10) %>%
    group_by(gopt, ctag) %>%
    #mutate(fcn = scale(fc)) %>% ungroup() %>%
    mutate(fcn = fc) %>% ungroup() %>%
    mutate(lab = sprintf("%.01f", fc)) %>%
    mutate(lab = str_remove(lab, '^0+')) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by='nid') %>%
    mutate(lgd = factor(lgd, levels=t_cfg$lgd))
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
    theme(axis.text.y = element_text(color=rev(tps$col))) +
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



