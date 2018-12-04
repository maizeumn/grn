#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("net", nargs=1, help="GRN file (*.rda)")
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument("--opt", default='tf',
               help="evaluation option [default: '%(default)s']")
args <- p$parse_args()

f_net = args$net
f_out = args$out
opt = args$opt
if( file.access(f_net) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_net))

source("~/projects/grn/src/functions.R")
x = load(f_net)

eval_gs <- function(f_net, gs) {
    #{{{
    tf = gs$tf; tfs = gs$tfs
    y = load(f_net)
    ctags_v = tfs %>% filter(reg.gid %in% rids) %>% pull(ctag)
    t_pr = tibble(); t_roc = tibble(); t_aupr = tibble(); t_auroc = tibble()
    for (ctag in ctags_v) {
        rid = tfs %>% filter(ctag == !!ctag) %>% pull(reg.gid)
        score = as.numeric(reg.mat[rid,])
        score[is.na(score)] = 0
        true.tgts = tf %>% filter(reg.gid == rid) %>% pull(tgt.gid)
        categ = rep(0, length(tids))
        names(categ) = tids
        categ[true.tgts] = 1
        if(max(score) == 0) score[1] = 0.1
        pr = pr.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        roc = roc.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        t_pr1 = pr$curve %>% as_tibble() %>%
            transmute(ctag = ctag,
                      recall = V1, precision = V2, score = V3)
        t_roc1 = roc$curve %>% as_tibble() %>%
            transmute(ctag = ctag,
                      TPR = V1, FPR = V2, score = V3)
        t_pr = rbind(t_pr, t_pr1)
        t_roc = rbind(t_roc, t_roc1)
        aupr = pr$auc.integral
        auroc = roc$auc
        t_aupr1 = tibble(ctag = ctag, aupr = aupr)
        t_auroc1 = tibble(ctag = ctag, auroc = auroc)
        t_aupr = rbind(t_aupr, t_aupr1)
        t_auroc = rbind(t_auroc, t_auroc1)
    }
    list('pr'=t_pr, 'roc'=t_roc, 'aupr'=t_aupr, 'auroc'=t_auroc)
    #}}}
}
eval_go_corncyc <- function(f_net, gs, net_sizes = c(1e4,5e4,1e5,5e5)) {
    #{{{
    y = load(f_net)
    set.seed(1026)
    tn1 = tn %>% mutate(permut='observed') %>% select(permut, reg.gid, tgt.gid)
    tn2 = tn %>% mutate(tgt.gid = sample(tgt.gid)) %>%
        mutate(permut='random') %>% select(permut,reg.gid,tgt.gid)
    tnc = tn1 %>% bind_rows(tn2)
    res = tibble()
    for (net_size in net_sizes) {
        tn1 = tnc %>%
            group_by(permut) %>%
            filter(row_number() <= net_size) %>% ungroup() %>%
            inner_join(gs$grp_f, by = c("tgt.gid" = 'gid')) %>%
            rename(grp_tag = ctag) %>%
            count(permut, grp_tag, grp, reg.gid) %>%
            rename(ng = n) %>%
            group_by(permut, grp_tag, grp) %>%
            summarise(n.reg = n(), n.tgt = sum(ng),
                      pairs.total = n.tgt * (n.tgt-1) / 2,
                      pairs.coreg = sum( ng*(ng-1)/2 ),
                      rich = pairs.coreg/pairs.total) %>%
            ungroup() %>%
            mutate(net_size = !!net_size) %>%
            select(net_size,everything()) %>%
            arrange(desc(rich))
        res = rbind(res, tn1)
    }
    res
    #}}}
}
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
    res = eval_gs(f_net, gs)
} else if (opt == 'go') {
    gs = read_gs()
    res = eval_go_corncyc(f_net, gs)
} else if (opt == 'briggs') {
    br = read_briggs()
    res = eval_briggs(f_net, br)
} else if (opt == 'biomap') {
    bm = read_biomap(opt='inbred')
    res = eval_biomap(f_net, bm)
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(res, file = f_out)


