#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("net", nargs=1, help="GRN file (*.rds)")
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument('-t', "--opt", default='tf',
               help="evaluation option [default: '%(default)s']")
p$add_argument('-m', "--permut", type='integer', default=100,
               help="number permutations to evaluate significance [default: %(default)s]")
p$add_argument('-p', "--thread", type='integer', default=1,
               help="number threads [default: %(default)s]")
args <- p$parse_args()

f_net = args$net
f_out = args$out
opt = args$opt
permut = args$permut
thread = args$thread
if( file.access(f_net) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_net))

source("~/projects/grn/src/functions.R")
require(future)
require(furrr)
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

eval_fun_ann_1 <- function(net_size=1e4, perm=0, tn, fun_ann) {
    #{{{
    tn = tn %>% filter(row_number() <= net_size)
    if(perm != 0) {
        set.seed(perm)
        tn = tn %>% mutate(tgt.gid = sample(tgt.gid))
    }
    if(perm %% 10 == 0)
        cat('net_size:', net_size, 'permut:', perm, '\n')
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

if(thread > 1) {
    plan(multiprocess, workers=thread)
    options(future.globals.maxSize=10e9)
}

#gs was already read
if (opt == 'bs') {
    #{{{ evaluate network stats, TF/target pairs, Y1H overlap
    require(PRROC)
    fi = file.path(dird, '08_y1h', '01.rds')
    y1h = readRDS(fi)
    res = eval_gs(f_net, gs, y1h)
    net_sizes=c(1e4,5e4,1e5,5e5)
    tf = gs$tf; tfbs = gs$tfbs; ko = gs$ko
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
} else if (opt == 'go') {
    #{{{
    fun_ann = gs$fun_ann
    n_permut=permut; net_size=1e6; nbin=10
    res = tn %>% slice(1:net_size) %>%
        mutate(score = as.integer(cut_number(score,nbin))) %>%
        select(reg.gid,tgt.gid,score) %>%
        group_by(score) %>% nest() %>% rename(tn=data) %>%
        crossing(perm=0:n_permut) %>%
        mutate(data = future_map2(perm, tn, eval_go_1, fun_ann=fun_ann, .progress=T)) %>%
        #mutate(data = map2(perm, tn, eval_go_1, fun_ann=fun_ann)) %>%
        mutate(enc.grp = map(data, 'enc.grp')) %>%
        mutate(enc.reg = map(data, 'enc.reg'))
    res1 = res %>% select(-tn, -data, -enc.reg) %>% unnest()
    res2 = res %>% select(-tn, -data, -enc.grp) %>% unnest()
    res10 = res1 %>% filter(perm==0) %>% select(-perm,-coreg)
    res20 = res2 %>% filter(perm==0) %>% select(-perm,-coreg)
    enrich_grp = res1 %>%
        group_by(score, ctag, grp) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup() %>%
        inner_join(res10, by=c('score','ctag','grp'))
    enrich_reg = res2 %>%
        inner_join(res20[,1:3], by=c('score','ctag','reg.gid')) %>%
        group_by(score, ctag, reg.gid) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup() %>%
        inner_join(res20, by=c('score','ctag','reg.gid'))
    enrich = res1 %>%
        group_by(score, ctag, perm) %>%
        summarise(n_grp=length(grp), coreg=sum(coreg)) %>% ungroup() %>%
        group_by(score, ctag, n_grp) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup()
    res = list(enrich=enrich, enrich_grp=enrich_grp, enrich_reg=enrich_reg)
    #}}}
} else if (opt == 'nv') {
    #{{{ evaluate natural variation datasets
    net_size=1e6; nbin=10
    nv = readRDS('~/projects/grn/data/06_deg/all.rds')
    tn0 = tn %>% slice(1:net_size) %>%
        mutate(score = as.integer(cut_number(score,nbin))) %>%
        mutate(pcc_sign=ifelse(pcc < 0, '-', '+')) %>%
        select(reg.gid,tgt.gid,score,pcc_sign)
    nv0 = nv %>% select(yid,cond,group1,group2,gid,DE,DEdir)
    res = tn0 %>% inner_join(nv0, by=c('reg.gid'='gid')) %>%
        rename(reg.DE=DE, reg.DEdir=DEdir) %>%
        inner_join(nv0, by=c('tgt.gid'='gid','yid','cond','group1','group2')) %>%
        rename(tgt.DE=DE, tgt.DEdir=DEdir) %>%
        mutate(DE_sign=ifelse(reg.DEdir==tgt.DEdir, '+', '-')) %>%
        mutate(consis = pcc_sign == DE_sign) %>%
        count(yid,cond,group1,group2,reg.DE,tgt.DE,score, consis)
    #}}}
} else if (opt == 'bm') {
    bm = read_biomap(opt='inbred')
    res = eval_biomap(f_net, bm)
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(res, file = f_out)


