#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("net", nargs=1, help="GRN file (*.rda)")
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
x = load(f_net)

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
eval_tf <- function(rid, tf, tids, reg.mat) {
    #{{{
    score = as.numeric(reg.mat[rid,])
    score[is.na(score)] = 0
    true.tgts = tf %>% filter(reg.gid == rid) %>% pull(tgt.gid)
    categ = rep(0, length(tids))
    names(categ) = tids
    categ[true.tgts] = 1
    if(max(score) == 0) score[1] = 0.1
    resR = roc.curve(scores.class0 = score, weights.class0 = categ, curve = T)
    resP = pr.curve(scores.class0 = score, weights.class0 = categ, curve = T)
    roc = resR$curve %>% as_tibble() %>%
        transmute(TPR = V1, FPR = V2, score = V3)
    prc = resP$curve %>% as_tibble() %>%
        transmute(recall = V1, precision = V2, score = V3)
    auroc = resR$auc
    auprc = resP$auc.integral
    list(roc=roc, prc=prc, auroc=auroc, auprc=auprc)
    #}}}
}
eval_y1h <- function(net_size=1e4, tn, reg.gids, tgt.gids) {
#{{{
    tn %>% filter(row_number() <= net_size) %>%
        filter(reg.gid %in% reg.gids, tgt.gid %in% tgt.gids) %>%
        select(reg.gid, tgt.gid, score)
#}}}
}

#' evaluate network stats, TF/target pairs, Y1H overlap
eval_gs <- function(f_net, gs, y1h, net_sizes=c(1e4,5e4,1e5,5e5)) {
    #{{{
    tf = gs$tf; tfs = gs$tfs
    y = load(f_net)
    tfstat = tfs %>% filter(reg.gid %in% rids) %>%
        mutate(res = map(reg.gid, eval_tf, tf=tf, tids=tids, reg.mat=reg.mat)) %>%
        mutate(auroc = map_dbl(res, 'auroc'), auprc = map_dbl(res, 'auprc'),
            roc = map(res, 'roc'), prc = map(res, 'prc')) %>%
        select(-res)
    # network properties
    nstat = tibble(net_size = net_sizes) %>%
        mutate(res = map(net_size, eval_netstat, tn)) %>%
        mutate(n.reg = map_int(res, 'n.reg'), n.tgt = map_int(res, 'n.tgt'),
            deg.reg = map(res, 'deg.reg'), deg.tgt = map(res, 'deg.tgt')) %>%
        select(-res)
    ystat = tibble(net_size = net_sizes) %>%
        mutate(res = map(net_size, eval_y1h, tn, y1h$reg.gids, y1h$tgt.gids))
    list(tfstat=tfstat, nstat=nstat, ystat=ystat)
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
    enc = tn %>%
        inner_join(fun_ann, by = c("tgt.gid" = 'gid')) %>%
        count(ctag, grp, reg.gid) %>%
        group_by(ctag, grp) %>%
        summarise(coreg = sum(n*(n-1)/2), n = sum(n)) %>%
        ungroup()
    enc
    #}}}
}

#' evaluate functional enrichment of GO CornCyc Y1H
eval_fun_ann <- function(f_net, gs, n_permut=permut, net_sizes=c(1e4,5e4,1e5,5e5)) {
    #{{{
    fun_ann = gs$fun_ann
    y = load(f_net)
    #ptm <- proc.time()
    #x = proc.time() - ptm
    #cat(x[3], '\n')
    res = crossing(net_size=net_sizes, perm=0:n_permut) %>%
        mutate(res = map2(net_size, perm, eval_fun_ann_1, tn=tn, fun_ann=fun_ann)) %>%
        unnest() %>%
        filter(n >= 2)
    enrich_term = res %>%
        group_by(net_size, ctag, grp, n) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup()
    enrich = res %>%
        group_by(net_size, ctag, perm) %>%
        summarise(n_grp=length(grp), coreg=sum(coreg)) %>% ungroup() %>%
        group_by(net_size, ctag, n_grp) %>%
        summarise(fc = coreg[perm==0]/mean(coreg[perm!=0]),
            pval = sum(coreg[perm!=0]>coreg[perm==0]) / sum(perm!=0)) %>%
        ungroup()
    list(enrich=enrich, enrich_term=enrich_term)
    #}}}
}

#' evaluate briggs dataset
eval_briggs <- function(f_net, br, net_size=5e4) {
    #{{{
    y = load(f_net)
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
    x = load(f_net)
    tp = tn %>%
        filter(row_number() <= net_size) %>%
        inner_join(bm, by = c('reg.gid' = 'gid')) %>%
        mutate(r.cpm = asinh(CPM)) %>% select(-CPM) %>%
        inner_join(bm, by = c('tgt.gid'='gid','Genotype'='Genotype','Tissue'='Tissue')) %>%
        mutate(t.cpm = asinh(CPM)) %>%
        group_by(reg.gid, tgt.gid, Tissue) %>%
        summarise(pcc = cor(r.cpm, t.cpm)) %>% ungroup() %>%
        filter(!is.na(pcc))
    tp
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


