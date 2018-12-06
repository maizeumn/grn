#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument("--opt", default='tf',
               help="evaluation option [default: '%(default)s']")
args <- p$parse_args()

f_out = args$out
opt = args$opt

source("~/projects/grn/src/functions.R")
dirw = file.path(dird, '14_eval_sum')

if (opt == 'tf') {
    ev = th %>% mutate(fi = sprintf("%s/13_eval/%s_tf.rds", dird, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        mutate(roc = map(res, "roc"), pr = map(res, "pr"),
               auroc = map(res, "auroc"), aupr = map(res, "aupr"),
               nstat = map(res, 'nstat'), ystat = map(res, 'ystat')) %>%
        select(nid,roc,pr,auroc,aupr,nstat,ystat)
} else if (opt == 'go') {
    ev = th %>% mutate(fi = sprintf("%s/13_eval/%s_go.rds", dird, nid)) %>%
        mutate(res = map(fi, readRDS))
} else if (opt == 'br') {
    ev = th %>% mutate(fi = sprintf("%s/13_eval/%s_br.rds", dird, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        select(nid,res) %>% unnest()
} else if (opt == 'bm') {
    ev = th %>% mutate(fi = sprintf("%s/13_eval/%s_bm.rds", dird, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        select(nid,res) %>% unnest()
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(ev, file = f_out)


