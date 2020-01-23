#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'GRN evaluation utilities')
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument("--eopt", default='tf', help="evaluation option [default: '%(default)s']")
p$add_argument("--gopt", default='rf', help="GRN option [default: '%(default)s']")
args <- p$parse_args()

f_out = args$out
gopt = args$gopt
eopt = args$eopt

source("~/projects/grn/src/functions.R")
dirw = file.path(dird, 'cache')
diri = file.path(dirw, '17_eval')

ev0 = t_cfg %>%
    mutate(fi = sprintf("%s/%s.%s.%s.rds", diri, eopt, gopt, nid)) %>%
    mutate(res = map(fi, readRDS))
if (eopt %in% c('bs','ko')) {
    ev = ev0 %>% unnest(res)
} else if (eopt == 'tf') {
    ev = ev0 %>%
        mutate(tf = map(res, 'tfstat'),
               tfbs = map(res, 'tfbsstat'),
               ko = map(res, 'kostat'),
               nstat = map(res, 'nstat'),
               ystat = map(res, 'ystat')) %>%
        select(nid,tf,tfbs,ko,nstat,ystat)
} else if (eopt == 'go') {
    ev = ev0 %>%
        mutate(enrich = map(res, "enrich"),
               enrich_grp = map(res, "enrich_grp"),
               enrich_reg = map(res, "enrich_reg")) %>%
        select(nid,enrich,enrich_grp,enrich_reg)
} else if (eopt == 'nv') {
    ev = ev0 %>%
        select(nid,res) %>% unnest()
} else if (eopt == 'br') {
    ev = ev0 %>%
        select(nid,res) %>% unnest()
    evr = ev %>%
        mutate(p.drc = ifelse(pcc < 0, -1, 1)) %>%
        mutate(b.drc = ifelse(reg.DEdir==tgt.DEdir, 1, -1)) %>%
        select(nid,tissue,reg.gid,tgt.gid,p.drc,b.drc,reg.DE,tgt.DE)
    ev = evr
} else if (eopt == 'bm') {
    ev = ev0 %>%
        select(nid,res) %>% unnest(res)
} else {
    stop(sprintf("unknown option: %s\n", opt))
}
saveRDS(ev, file = f_out)


