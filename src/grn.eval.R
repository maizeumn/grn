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
require(furrr)
require(purrr)
y = readRDS(f_net)
rids=y$rids; tids=y$tids; tn=y$tn

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
eval_tf1 <- function(ti, net_size, tn, rids, tids, fpr=.1) {
    #{{{
    rids0 = rids[rids %in% ti$reg.gid]
    tn0 = tn %>% filter(row_number() <= net_size)
    tt = complete_tn(tn0, tids, rids0)
    if(nrow(tt) == 0) return(list(auroc=NA,pval=NA,auroc0=NA,auprc=NA))
    tt = tt %>% filter(reg.gid %in% rids0) %>%
        filter(reg.gid != tgt.gid) %>%
        mutate(score = as.numeric(score)) %>%
        mutate(score = ifelse(is.na(score), 0, score)) %>%
        arrange(desc(score))
    if(max(tt$score) == 0) tt$score[1] = 0.1
    to = tt %>%
        left_join(ti, by = c('reg.gid','tgt.gid')) %>%
        replace_na(list(response=0))
    # auroc & auprc
    resR = roc.curve(scores.class0=to$score, weights.class0=to$response)
    resP = pr.curve(scores.class0=to$score, weights.class0=to$response)
    auroc0 = resR$auc
    auprc = resP$auc.integral
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
    list(auroc=auroc, pval=pval, auroc0=auroc0, auprc=auprc)
    #}}}
}
eval_go_1 <- function(perm=0, tn, fun_ann) {
    #{{{
    if(perm != 0) {
        set.seed(perm)
        tn = tn %>% mutate(tgt.gid = sample(tgt.gid))
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

if(thread > 1) {
    plan(multiprocess, workers=thread)
    options(future.globals.maxSize=10e9)
}

if (opt %in% c('bs','ko')) {
    #{{{ ChIP-Seq / DAP-Seq / PWM-TFBS | knock-out mutant
    require(PRROC)
    require(pROC)
    fpr = .1; net_sizes = c(1e5, 1e6, 1e7)
    x1 = gs[[opt]] %>% filter(reg.gid %in% rids)
    if(nrow(x1) == 0) {
        res = NULL
    } else {
        if(opt == 'bs') {
            x2 = x1 %>%
                mutate(tf = reg.gid, response = 1) %>%
                group_by(ctag, tf) %>% nest()
        } else if(opt == 'ko') {
            x2 = x1 %>%
                unnest(ds) %>% rename(tgt.gid=gid) %>%
                mutate(response = ifelse(padj < .01, 1, 0)) %>%
                select(-padj, -log2fc) %>%
                group_by(yid,author,tf,tissue) %>% nest()
        }
        res = x2 %>% crossing(net_size = net_sizes) %>%
            mutate(r = map2(data, net_size, eval_tf1, tn=tn, rids=rids, tids=tids, fpr=fpr)) %>%
            mutate(auroc=map_dbl(r, 'auroc')) %>%
            mutate(pval=map_dbl(r, 'pval')) %>%
            mutate(auroc0=map_dbl(r, 'auroc0')) %>%
            mutate(auprc=map_dbl(r, 'auprc')) %>%
            select(-data, -r)
    }
    #}}}
} else if (opt == 'go') {
    #{{{
    fun_ann = gs$fun_ann %>% filter(!str_detect(ctag, "^GO_.*_[CF]$"))
    n_permut=permut; net_size=1e6; nbin=10
    res = tn %>% slice(1:net_size) %>%
        mutate(score = as.integer(cut_number(score,nbin))) %>%
        select(reg.gid,tgt.gid,score) %>%
        group_by(score) %>% nest() %>% rename(tn=data) %>%
        crossing(perm=0:n_permut) %>%
        mutate(data = future_map2(perm, tn, eval_go_1, fun_ann=fun_ann, .progress=T)) %>%
        #mutate(data = map2(perm, tn, eval_go_1, fun_ann=fun_ann)) %>%
        mutate(enc.grp = map(data, 'enc.grp')) %>%
        mutate(enc.reg = map(data, 'enc.reg')) %>%
        select(-tn, -data)
    res1 = res %>% select(-enc.reg) %>% unnest(cols = c(enc.grp))
    res2 = res %>% select(-enc.grp) %>% unnest(cols = c(enc.reg))
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
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(res, file = f_out)


