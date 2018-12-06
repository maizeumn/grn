#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Run GENIE3')
parser$add_argument("exp", nargs=1, help="Expression file (*.pkl)")
parser$add_argument("net", nargs=1, help="GRN (genie3) file (*.pkl)")
parser$add_argument("out", nargs=1, help="output file (*.rda)")
parser$add_argument("--net_size", type="integer", default=1e6,
                    help="num. edges in net to keep [default %(default)s]")
args <- parser$parse_args()

f_exp = args$exp
f_net = args$net
f_out = args$out
net_size = args$net_size
if( file.access(f_exp) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_exp))
if( file.access(f_net) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_net))

require(tidyverse)
require(reticulate)
x = py_load_object(normalizePath(f_exp))
VIM = py_load_object(normalizePath(f_net))
tids = unique(x[[2]]); rids = unique(x[[3]])
stopifnot(dim(VIM)[1] == length(tids))
rownames(VIM) = tids
colnames(VIM) = tids
trash = VIM[tids[!tids %in% rids],]
stopifnot(sum(trash) == 0)
reg.mat = VIM[rids,]
tn = as_tibble(reg.mat) %>% mutate(reg.gid = rids) %>%
    gather(tgt.gid, score, -reg.gid) %>%
    filter(reg.gid != tgt.gid) %>%
    arrange(desc(score)) %>%
    filter(row_number() <= net_size)
save(rids, tids, reg.mat, tn, file = f_out)

