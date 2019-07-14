source("functions.R")
require(PRROC)
dirw = file.path(dird, '14_eval_sum')
gopts = c("rf",'et','xgb')
colmap = pal_aaas()(3)
names(colmap) = gopts

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
fp = file.path(dirw, "03.tsne.pdf")
ggsave(p_tsne, filename = fp, width=10, height=10)
#}}}
#}}}

#{{{ knockout pval
fd = file.path(dird, '07_mutants', 'degs.rds')
ds = readRDS(fd) %>% filter(gene_alias != 'P1')
dss = ds %>% unnest() %>% group_by(gene_alias, Tissue) %>%
    summarise(n_tot=n(), n_de=sum(padj<.01), prop_de=n_de/n_tot) %>%
    ungroup() %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags = dss$ctag

eopt = 'tf'
tvk = tibble(gopt = gopts) %>%
    mutate(fi = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)) %>%
    mutate(data = map(fi, readRDS)) %>%
    unnest() %>% select(nid, gopt, ko) %>% unnest() %>%
    filter(!is.na(pval)) %>%
    inner_join(t_cfg, by = 'nid')

tv0 = tvk %>% inner_join(dss, by=c('gene_alias','Tissue'))
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
tp = tv0 %>%
    mutate(lab = ifelse(pval<.05, number(-log10(pval),accuracy=2), '')) %>%
    mutate(pval=-log10(pval)) %>%
    select(gopt,nid,ctag,pval,lab) %>%
    mutate(gopt = str_to_upper(gopt)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
pval.max = max(tp$pval)
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=pval)) +
    geom_tile() +
    geom_text(aes(label=lab, col=abs(pval-pval.max/2)<pval.max/4), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="-log10(P-value)", begin = .4) +
    scale_color_manual(values = c('white','black')) +
    scale_fill_gradientn(name = '-log10(P-value)', colors = cols100) +
    facet_grid(.~gopt) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(t_cfg$col)))
fo = file.path(dirw, '02.gopt.3.ko.pval.pdf')
ggsave(p1, file = fo, width = 15, height = 8)
#}}}

#{{{ GO evaluation
eopt = 'go'
tvg = tibble(gopt = gopts) %>%
    mutate(fi = sprintf("%s/%s.%s.rds", dirr, gopt, eopt)) %>%
    mutate(data = map(fi, readRDS)) %>%
    unnest() %>% select(-fi) %>%
    inner_join(t_cfg, by = 'nid')

#{{{ enrichment
net_size = 5e4
ctags = c("GO_HC", "GO_arabidopsis","CornCyc")
fo = file.path(dirw, "02.gopt.6.go.pdf")
tp1 = tvg %>% select(nid, gopt, lgd, enrich) %>% unnest() %>%
    filter(net_size == !!net_size)
tp2 = tvg %>% select(nid, gopt, enrich_grp) %>% unnest() %>%
    filter(net_size == !!net_size) %>%
    filter(n >= 10) %>%
    group_by(nid, gopt, net_size, ctag) %>%
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


#{{{ ## eval knockout mutant RNA-Seq at once
ev_tf = readRDS(fi_tf)

fd = file.path(dird, '07_mutants', 'degs.rds')
ds = readRDS(fd) %>% filter(gene_alias != 'P1')
dss = ds %>% unnest() %>% group_by(gene_alias, Tissue) %>%
    summarise(n_tot=n(), n_de=sum(padj<.01), prop_de=n_de/n_tot) %>%
    ungroup() %>%
    mutate(ctag=sprintf("%s [%s] [%s] [%s]", gene_alias, Tissue, number(n_de), percent(prop_de)))
ctags = dss$ctag

eval_tf_1 <- function(gene_id,gene_alias,Tissue,t_ds, tids,tn) {
    #{{{
    rid = gene_id
    gids = tids
    #
    if(! rid %in% tn$reg.gid) {
        list(auroc=NA, auprc=NA, pval=NA)
    } else {
        tr = tn %>% filter(reg.gid == rid, reg.gid != tgt.gid) %>%
            select(gid = tgt.gid, score) %>%
            mutate(score = as.numeric(score)) %>%
            replace_na(list(score=0))
        if(max(tr$score) == 0) tr$score[1] = 0.1
        tt = t_ds %>% mutate(weight=padj < .01) %>% select(gid, weight)
        to = tt %>% filter(gid %in% gids) %>%
            left_join(tr, by = 'gid') %>%
            replace_na(list(score=0))
        resR = roc.curve(scores.class0=to$score, weights.class0=to$weight)
        resP = pr.curve(scores.class0=to$score, weights.class0=to$weight)
        scores1 = to %>% filter(weight) %>% pull(score)
        scores2 = to %>% filter(!weight) %>% pull(score)
        res = wilcox.test(scores1, scores2, alternative='greater')
        pval = res$p.value
        auroc = resR$auc
        auprc = resP$auc.integral
        list(auroc=auroc, auprc=auprc, pval=pval)
    }
    #}}}
}
eval_tf <- function(tids,tn, ds) {
    #{{{
    ds %>%
        mutate(res = pmap(list(gene_id,gene_alias,Tissue,ds), eval_tf_1,
                          tids=!!tids,tn=!!tn)) %>%
        mutate(auroc=map_dbl(res,'auroc'), auprc=map_dbl(res,'auprc'),
               pval=map_dbl(res,'pval')) %>%
        select(-ds,-res)
    #}}}
}

tv = ev_tf %>% #filter(nid=='n18a') %>%
    mutate(r=map2(tids,tn, eval_tf, ds=ds)) %>%
    select(nid, r) %>% unnest() %>%
    filter(!is.na(auroc))
#}}}

