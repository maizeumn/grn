#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'convert GRN format')
parser$add_argument("exp", nargs=1, help="expression matrix (*.tsv)")
parser$add_argument("net", nargs=1, help="GRN result file (*.pkl)")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
parser$add_argument("--net_size", type="integer", default=2e6,
                    help="num. edges in net to keep [default %(default)s]")
args <- parser$parse_args()

f_exp = args$exp
f_net = args$net
f_out = args$out
net_size = args$net_size
if( file.access(f_net) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_net))

require(tidyverse)
require(reticulate)
x = py_load_object(normalizePath(f_net))
rids = x[[1]]
tids = x[[2]]
#tn = as_tibble(py_to_r(x[[3]]))
tn = as_tibble(x[[3]])
mse = as.numeric(x[[4]])
#stopifnot(nrow(tn) == length(rids)*length(tids))
tn = tn %>% rename(reg.gid = rid, tgt.gid = tid) %>%
    arrange(desc(score)) %>%
    filter(row_number() <= net_size)
t_mse = tibble(gid = tids, mse = mse)

tx = read_tsv(f_exp) %>%
    gather(sid, cpm, -gid) %>%
    mutate(cpm = asinh(cpm)) %>%
    group_by(gid) %>%
    summarise(v = list(cpm)) %>% ungroup()

tn = tn %>% inner_join(tx, by = c('reg.gid'='gid')) %>%
    rename(reg.v = v) %>%
    inner_join(tx, by = c('tgt.gid'='gid')) %>%
    rename(tgt.v = v) %>%
    mutate(res1 = map2(reg.v, tgt.v, cor.test, 'two.sided', 'pearson'),
        res2 = map2(reg.v, tgt.v, cor.test, 'two.sided', 'spearman')) %>%
    mutate(pcc = map_dbl(res1, 'estimate'),
        pval.p = map_dbl(res1, 'p.value'),
        spc = map_dbl(res2, 'estimate'),
        pval.s = map_dbl(res2, 'p.value')) %>%
    select(reg.gid, tgt.gid, score, pcc, pval.p, spc, pval.s)

res = list(rids=rids, tids=tids, mse=t_mse, tn=tn)
saveRDS(res, file=f_out)

