#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'merge GRNs')
p$add_argument("out", nargs=1, help="output file (*.rds)")
p$add_argument("--gopt", default='rf',
               help="GRN option [default: '%(default)s']")
args <- p$parse_args()

f_out = args$out
gopt = args$gopt
stopifnot(gopt %in% c('rf','et','xgb'))

source("~/projects/grn/src/functions.R")
dirw = file.path(dird, 'cache')
diri = file.path(dirw, '15_grn_rds')

tn = t_cfg %>% mutate(fi = sprintf("%s/%s.%s.rds", diri, gopt, nid)) %>%
    mutate(res = map(fi, readRDS)) %>%
    mutate(rids = map(res, 'rids'),
           tids = map(res, 'tids'),
           mse = map(res, 'mse'),
           tn = map(res, 'tn')) %>%
    select(-fi, -res)

saveRDS(tn, file = f_out)

tn1 = tn %>% mutate(tn = map(tn, slice, 1:5e4))
tn2 = tn %>% mutate(tn = map(tn, slice, 1:1e5))
tn3 = tn %>% mutate(tn = map(tn, slice, 1:5e5))
tn4 = tn %>% mutate(tn = map(tn, slice, 1:1e6))

fo1 = str_replace(f_out, ".rds$", ".50k.rds")
fo2 = str_replace(f_out, ".rds$", ".100k.rds")
fo3 = str_replace(f_out, ".rds$", ".500k.rds")
fo4 = str_replace(f_out, ".rds$", ".1m.rds")
saveRDS(tn1, file = fo1)
saveRDS(tn2, file = fo2)
saveRDS(tn3, file = fo3)
saveRDS(tn4, file = fo4)



