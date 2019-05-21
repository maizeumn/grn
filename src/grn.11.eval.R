#{{{
source("functions.R")
require(PRROC)
dirw = file.path(dird, '14_eval_sum')
tsyn = read_syn(gcfg)
#
gs = read_gs()
fi_tf = file.path(dirr, '01.tf.rds')
fi_br = file.path(dirr, '01.br.rds')
fi_go = file.path(dirr, '01.go.rds')
fi_bm = file.path(dirr, '01.bm.rds')
#}}}

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

#{{{ aupr/auroc bar-plot
tv0 = tv %>% inner_join(dss, by=c('gene_alias','Tissue'))
levs = c("AUROC","AUPR")
cols100 = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
tps = th %>% filter(nid %in% tv0$nid)
tp = tv0 %>%
    select(nid,ctag,AUROC=auroc,AUPR=auprc) %>%
    gather(type, auc, -nid, -ctag) %>%
    group_by(type) %>%
    mutate(auc.norm = as.numeric(scale(auc))) %>% ungroup() %>%
    mutate(lab = str_remove(sprintf("%.03f", auc), '^0+')) %>%
    mutate(type = factor(type, levels = levs)) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(tps$lgd)))
#
#tpx = gs$tnk %>% filter(ctag %in% ctags) %>%
    #count(ctag) %>% mutate(lab=sprintf("%s (%d)",ctag,n)) %>%
    #mutate(ctag=factor(ctag, levels=ctags)) %>% arrange(ctag)
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=auc.norm)) +
    geom_tile() +
    geom_text(aes(label=lab), hjust=.5, size=2) +
    #scale_x_discrete(breaks=tpx$ctag, labels=tpx$lab, expand=expand_scale(mult=c(0,0))) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="Area Under Curve (AUC)", begin = .4) +
    scale_fill_gradientn(name = 'Area Under Curve (AUC)', colors = cols100) +
    facet_wrap(~type, nrow = 1) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=7)) +
    theme(axis.text.y = element_text(color=rev(tps$col)))
fo = file.path(dirw, '05.tf.auc.pdf')
ggsave(p1, file = fo, width = 10, height = 8)
#}}}

#{{{ pval plot
tv0 = tv %>% inner_join(dss, by=c('gene_alias','Tissue'))
cols100 = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
tps = th %>% filter(nid %in% tv0$nid)
tp = tv0 %>% mutate(lab = ifelse(pval < .05, scientific(pval,digits=2), 'ns')) %>%
    select(nid,ctag,pval,lab) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th, by = 'nid') %>%
    mutate(lgd = factor(lgd, levels=rev(tps$lgd)))
p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=-log(pval))) +
    geom_tile() +
    geom_text(aes(label=lab), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(expand=c(0,0)) +
    #scale_fill_viridis(name="Area Under Curve (AUC)", begin = .4) +
    scale_fill_gradientn(name = '-log(pvalue)', colors = cols100) +
    otheme(strip.size=8, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=7)) +
    theme(axis.text.y = element_text(color=rev(tps$col)))
fo = file.path(dirw, '05.tf.pval.pdf')
ggsave(p1, file = fo, width = 11, height = 8)
#}}}


