#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("net", nargs=1, help="GRN file (*.rds)")
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument("--opt", default='tf',
               help="evaluation option [default: '%(default)s']")
p$add_argument("--permut", type='integer', default=100,
               help="number permutations to evaluate significance [default: %(default)s]")
args <- p$parse_args()

f_net = args$net
f_out = args$out
opt = args$opt
permut = args$permut
if( file.access(f_net) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_net))

source("~/projects/grn/src/functions.R")
y = readRDS(f_net)
rids=y$rids; tids=y$tids; tn=y$tn

eval_netstat <- function(net_size=1e4, tn) {
#{{{
    tn = tn %>% filter(row_number() <= net_size)
    n.reg = length(unique(tn$reg.gid))
    n.tgt = length(unique(tn$tgt.gid))
    deg.reg = tn %>% count(reg.gid)
    deg.tgt = tn %>% count(tgt.gid)
    list(n.reg = n.reg, n.tgt = n.tgt, deg.reg = deg.reg, deg.tgt = deg.tgt)
#}}}
}
eval_tf <- function(ctag, tnk, rids, tids, tn, net_size=Inf) {
    #{{{
    tnk0 = tnk %>% filter(ctag == !!ctag) %>% select(-ctag) %>%
        mutate(weight = T)
    rids0 = rids[rids %in% tnk0$reg.gid]
    tt = tn %>% filter(reg.gid %in% rids0) %>%
        filter(reg.gid != tgt.gid) %>%
        mutate(score = as.numeric(score)) %>%
        mutate(score = ifelse(is.na(score), 0, score)) %>%
        arrange(desc(score)) %>%
        filter(row_number() <= net_size)
    if(max(tt$score) == 0) tt$score[1] = 0.1
    to = tt %>%
        left_join(tnk0, by = c('reg.gid','tgt.gid')) %>%
        replace_na(list(weight=F))
    resR = roc.curve(scores.class0=to$score, weights.class0=to$weight)
    resP = pr.curve(scores.class0=to$score, weights.class0=to$weight)
    scores1 = to %>% filter(weight) %>% pull(score)
    scores2 = to %>% filter(!weight) %>% pull(score)
    pval = NA
    if(length(scores1) > 0 & length(scores2) > 0)
        pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
    #roc = resR$curve %>% as_tibble() %>%
        #transmute(TPR = V1, FPR = V2, score = V3)
    #prc = resP$curve %>% as_tibble() %>%
        #transmute(recall = V1, precision = V2, score = V3)
    auroc = resR$auc
    auprc = resP$auc.integral
    list(auroc=auroc, auprc=auprc, pval=pval)
    #}}}
}
eval_y1h <- function(net_size=1e4, tn, reg.gids, tgt.gids) {
#{{{
    tn %>% filter(row_number() <= net_size) %>%
        filter(reg.gid %in% reg.gids, tgt.gid %in% tgt.gids) %>%
        select(reg.gid, tgt.gid, score)
#}}}
}
eval_ko <- function(gene_id,gene_alias,Tissue,t_ds, tids,tn) {
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
        pval = NA
        if(length(scores1) > 0 & length(scores2) > 0)
            pval = wilcox.test(scores1, scores2, alternative='greater')$p.value
        auroc = resR$auc
        auprc = resP$auc.integral
        list(auroc=auroc, auprc=auprc, pval=pval)
    }
    #}}}
}

#' evaluate network stats, TF/target pairs, Y1H overlap
eval_gs <- function(f_net, gs, y1h, net_sizes=c(1e4,5e4,1e5,5e5)) {
    #{{{
    tf = gs$tf; tfbs = gs$tfbs; ko = gs$ko
    y = readRDS(f_net)
    rids=y$rids; tids=y$tids; tn=y$tn
    tfstat = tf %>% filter(reg.gid %in% rids) %>%
        distinct(ctag) %>%
        mutate(res=map(ctag, eval_tf, tnk=tf, rids=rids, tids=tids, tn=tn)) %>%
        mutate(auroc=map_dbl(res,'auroc'), auprc=map_dbl(res,'auprc'),
               pval=map_dbl(res, 'pval')) %>%
        select(-res)
    tfbsstat = tfbs %>% filter(reg.gid %in% rids) %>%
        distinct(ctag) %>%
        mutate(res=map(ctag, eval_tf, tnk=tfbs, rids=rids, tids=tids, tn=tn)) %>%
        mutate(auroc=map_dbl(res,'auroc'), auprc=map_dbl(res,'auprc'),
               pval=map_dbl(res, 'pval')) %>%
        select(-res)
    kostat = ko %>%
        mutate(res = pmap(list(gene_id,gene_alias,Tissue,ds), eval_ko,
                          tids=!!tids,tn=!!tn)) %>%
        mutate(auroc=map_dbl(res,'auroc'), auprc=map_dbl(res,'auprc'),
               pval=map_dbl(res,'pval')) %>%
        select(-ds,-res)
    # network properties
    nstat = tibble(net_size = net_sizes) %>%
        mutate(res = map(net_size, eval_netstat, tn)) %>%
        mutate(n.reg = map_int(res, 'n.reg'), n.tgt = map_int(res, 'n.tgt'),
            deg.reg = map(res, 'deg.reg'), deg.tgt = map(res, 'deg.tgt')) %>%
        select(-res)
    ystat = tibble(net_size = net_sizes) %>%
        mutate(res = map(net_size, eval_y1h, tn, y1h$reg.gids, y1h$tgt.gids))
    list(tfstat=tfstat, tfbsstat=tfbsstat, kostat=kostat, nstat=nstat, ystat=ystat)
    #}}}
}

eval_fun_ann_1 <- function(net_size=1e4, perm=0, tn, fun_ann) {
    #{{{
    tn = tn %>% filter(row_number() <= net_size)
    if(perm != 0) {
        set.seed(perm)
        tn = tn %>% mutate(tgt.gid = sample(tgt.gid))
    }
    if(perm %% 10 == 0) {
        cat('net_size:', net_size, 'permut:', perm, '\n')
    }
    tn2 = tn %>%
        inner_join(fun_ann, by = c("tgt.gid" = 'gid')) %>%
        count(ctag, grp, reg.gid)
    enc1 = tn2 %>%
        arrange(ctag, grp, desc(n), reg.gid) %>%
        group_by(ctag, grp) %>%
        summarise(coreg = sum(n*(n-1)/2),
                  max.reg.gid=reg.gid[1], max.reg.size=n[1], n=sum(n)) %>%
        ungroup()
    enc2 = tn2 %>%
        arrange(ctag, reg.gid, desc(n), grp) %>%
        group_by(ctag, reg.gid) %>%
        summarise(coreg = sum(n*(n-1)/2),
                  max.grp=grp[1], max.grp.size=n[1], n=sum(n)) %>%
        ungroup()
    list(enc.grp = enc1, enc.reg = enc2)
    #}}}
}
eval_fun_ann_1a <- function(net_size=1e4, perm=0, ctag='GO_HC', tn, fun_ann) {
    #{{{
    tn = tn %>% filter(row_number() <= net_size)
    fun_ann = fun_ann %>% filter(ctag == !!ctag) %>% select(-ctag)
    if(perm != 0) {
        set.seed(perm)
        tn = tn %>% mutate(tgt.gid = sample(tgt.gid))
    }
    if(perm %% 10 == 0 & ctag == 'GO_HC') {
        cat('net_size:', net_size, 'permut:', perm, '\n')
    }
    tn %>%
        inner_join(fun_ann, by = c("tgt.gid" = 'gid')) %>%
        count(reg.gid)
    tn2 = tn %>%
        inner_join(fun_ann, by = c("tgt.gid" = 'gid')) %>%
        count(grp, reg.gid)
    enc1 = tn2 %>%
        arrange(grp, desc(n), reg.gid) %>%
        group_by(grp) %>%
        summarise(coreg = sum(n*(n-1)/2),
                  max.reg.gid=reg.gid[1], max.reg.size=n[1], n=sum(n)) %>%
        ungroup()
    enc2 = tn2 %>%
        arrange(reg.gid, desc(n), grp) %>%
        group_by(reg.gid) %>%
        summarise(coreg = sum(n*(n-1)/2),
                  max.grp=grp[1], max.grp.size=n[1], n=sum(n)) %>%
        ungroup()
    list(enc.grp = enc1, enc.reg = enc2)
    #}}}
}

#' evaluate functional enrichment of GO CornCyc Y1H
eval_fun_ann <- function(f_net, gs, n_permut=permut, net_sizes=c(1e4,5e4,1e5,5e5)) {
    #{{{
    fun_ann = gs$fun_ann
    #ctags = fun_ann %>% distinct(ctag) %>% pull(ctag)
    y = readRDS(f_net)
    rids=y$rids; tids=y$tids; tn=y$tn
    #ptm <- proc.time()
    #x = proc.time() - ptm
    #cat(x[3], '\n')
    res = crossing(net_size=net_sizes, perm=0:n_permut) %>%
        mutate(data = map2(net_size, perm, eval_fun_ann_1, tn=tn, fun_ann=fun_ann))
    res1 = res %>% mutate(data=map(data, 'enc.grp')) %>% unnest()
    res2 = res %>% mutate(data=map(data, 'enc.reg')) %>% unnest()
    res10 = res1 %>% filter(perm==0) %>% select(-perm,-coreg)
    res20 = res2 %>% filter(perm==0) %>% select(-perm,-coreg)
    enrich_grp = res1 %>%
        group_by(net_size, ctag, grp) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup() %>%
        inner_join(res10, by=c('net_size','ctag','grp'))
    enrich_reg = res2 %>%
        inner_join(res20[,1:3], by=c('net_size','ctag','reg.gid')) %>%
        group_by(net_size, ctag, reg.gid) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup() %>%
        inner_join(res20, by=c('net_size','ctag','reg.gid'))
    enrich = res1 %>%
        group_by(net_size, ctag, perm) %>%
        summarise(n_grp=length(grp), coreg=sum(coreg)) %>% ungroup() %>%
        group_by(net_size, ctag, n_grp) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup()
    list(enrich=enrich, enrich_grp=enrich_grp, enrich_reg=enrich_reg)
    #}}}
}

#' evaluate briggs dataset
eval_briggs <- function(f_net, br, net_size=5e4) {
    #{{{
    y = readRDS(f_net)
    rids=y$rids; tids=y$tids; tn=y$tn
    res = tibble()
    for (tissue in br$tissues) {
        t_de = br$de %>% filter(Tissue == tissue, gid %in% tids) %>% select(-Tissue)
        gids = t_de$gid
        rids1 = rids[rids %in% gids]
        tids1 = tids[tids %in% gids]
        tn1 = tn %>% filter(reg.gid %in% rids1, tgt.gid %in% tids1) %>%
            filter(row_number() <= net_size) %>%
            inner_join(t_de, by = c('reg.gid' = 'gid')) %>%
            rename(reg.DE = DE, reg.DEdir = DEdir) %>%
            inner_join(t_de, by = c('tgt.gid' = 'gid')) %>%
            rename(tgt.DE = DE, tgt.DEdir = DEdir) %>%
            mutate(tissue = !!tissue) %>%
            select(tissue, everything())
        res = rbind(res, tn1)
    }
    res
    #}}}
}

#' evaluate biomap dataset
eval_biomap <- function(f_net, bm, net_size=5e4) {
    #{{{
    gts_exl = c("B73","Mo17")
    gts_exl = c()
    tm = bm %>%
        filter(! Genotype %in% gts_exl) %>%
        #mutate(CPM = asinh(CPM)) %>%
        arrange(gid, Tissue, Genotype) %>%
        group_by(gid, Tissue) %>%
        summarise(v = list(CPM)) %>% ungroup()
    y = readRDS(f_net)
    rids=y$rids; tids=y$tids; tn=y$tn
    if(!'pcc' %in% colnames(tn)) tn = tn %>% mutate(pcc = 1)
    tn1 = tn %>%
        filter(row_number() <= net_size) %>%
        mutate(p.drc = ifelse(pcc < 0, -1, 1)) %>%
        mutate(simu = F) %>%
        select(simu,reg.gid,tgt.gid,p.drc)
    # permutation
    set.seed(001)
    tn2 = tn1 %>% mutate(simu = T, tgt.gid = sample(tgt.gid))
    tp = tn1 %>%
        bind_rows(tn2) %>%
        inner_join(tm, by = c('reg.gid'='gid')) %>%
        rename(reg.v = v) %>%
        inner_join(tm, by = c('tgt.gid'='gid','Tissue'='Tissue')) %>%
        rename(tgt.v = v) %>%
        mutate(bm.pval.spc = pmap_dbl(list(reg.v, tgt.v, p.drc), eval_bm_spc),
               bm.pval.spe = pmap_dbl(list(reg.v, tgt.v, p.drc), eval_bm_spe)) %>%
        select(-reg.v, -tgt.v)
    tp
    #}}}
}
eval_bm_spc <- function(reg.v, tgt.v, p.drc) {
    #{{{
    alt.hyp = ifelse(p.drc == -1, 'less', 'greater')
    res = cor.test(reg.v, tgt.v, alt.hyp, 'spearman')
    pval = res$p.value
    pval
    #}}}
}
eval_bm_spe <- function(reg.v, tgt.v, p.drc) {
    #{{{
    alt.hyp = ifelse(p.drc == -1, 'greater', 'less')
    idxs1 = which(reg.v <= .1)
    idxs2 = which(reg.v >= 1)
    pval = NA
    if(length(idxs1)>=3 & length(idxs2)>=3) {
        #res = t.test(tgt.v[idxs1], tgt.v[idxs2], alternative=alt.hyp)
        #res = t.test(tgt.v[idxs1], tgt.v[idxs2])
        res = wilcox.test(tgt.v[idxs1], tgt.v[idxs2], alternative=alt.hyp)
        pval = res$p.value
    }
    pval
    #}}}
}

if (opt == 'tf') {
    require(PRROC)
    gs = read_gs()
    fi = file.path(dird, '08_y1h', '01.rds')
    y1h = readRDS(fi)
    res = eval_gs(f_net, gs, y1h)
} else if (opt == 'go') {
    gs = read_gs()
    res = eval_fun_ann(f_net, gs)
} else if (opt == 'br') {
    br = read_briggs()
    res = eval_briggs(f_net, br)
} else if (opt == 'bm') {
    bm = read_biomap(opt='inbred')
    res = eval_biomap(f_net, bm)
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(res, file = f_out)


