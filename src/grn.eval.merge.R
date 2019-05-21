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
dirw = file.path(dird, 'cache')
diri = file.path(dirw, '17_eval')

if (opt == 'tf') {
    ev = th %>% mutate(fi = sprintf("%s/%s_tf.rds", diri, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        mutate(tfstat = map(res, 'tfstat'),
               nstat = map(res, 'nstat'),
               ystat = map(res, 'ystat')) %>%
        mutate(tfstat = map(tfstat, select, ctag, auroc, auprc)) %>%
        select(nid,tfstat,nstat,ystat)
} else if (opt == 'go') {
    ev = th %>% mutate(fi = sprintf("%s/%s_go.rds", diri, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        mutate(enrich = map(res, "enrich"),
               enrich_grp = map(res, "enrich_grp"),
               enrich_reg = map(res, "enrich_reg")) %>%
        select(nid,enrich,enrich_grp,enrich_reg)
} else if (opt == 'br') {
    ev = th %>% mutate(fi = sprintf("%s/%s_br.rds", diri, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        select(nid,res) %>% unnest()
    evr = ev %>%
        mutate(p.drc = ifelse(pcc < 0, -1, 1)) %>%
        mutate(b.drc = ifelse(reg.DEdir==tgt.DEdir, 1, -1)) %>%
        select(nid,tissue,reg.gid,tgt.gid,p.drc,b.drc,reg.DE,tgt.DE)
    ev = evr
} else if (opt == 'bm') {
    ev = th %>% mutate(fi = sprintf("%s/%s_bm.rds", diri, nid)) %>%
        mutate(res = map(fi, readRDS)) %>%
        select(nid,res) %>% unnest()
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(ev, file = f_out)


