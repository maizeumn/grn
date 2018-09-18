#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Run GENIE3')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                     help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                     dest="verbose", help="Print little output")
parser$add_argument("fi", nargs=1, help="Input (expression matrix) file")
parser$add_argument("fo", nargs=1, help="Output file")
parser$add_argument("--tf", default="~/data/genome/Zmays_v4/TF/10.tsv",
                    help = "TF IDs file [default %(default)s]")
parser$add_argument("--cpm", action = 'store_true', default = F,
                    help = "use CPM [default %(default)s (i.e., use FPKM)]")
parser$add_argument("-p", "--thread", type="integer", default=1,
                    help="Num. threads to use [default %(default)s]")
args <- parser$parse_args()

if( args$verbose ) { 
    write("writing some verbose output to standard error...\n", stderr()) 
}

fi = args$fi
fo = args$fo
ff = args$tf
thread = args$thread
use_cpm = args$cpm
if( file.access(fi) == -1 ) {
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))
}

require(tidyverse)
require(doParallel)
require(doRNG)
require(GENIE3)

tf = read_tsv(ff)
tf_ids = tf$gid

x = load(fi)

n_cond = length(unique(t_exp$condition))
gids = t_exp %>% group_by(gid) %>%
    summarise(n_cond_exp = sum(CPM >= 1)) %>%
    filter(n_cond_exp >= n_cond * .1) %>%
    pull(gid)
t_flt = t_exp %>% filter(gid %in% gids)

if(use_cpm) {
    t_flt = t_flt %>% mutate(exp.val = asinh(CPM))
} else {
    t_flt = t_flt %>% mutate(exp.val = asinh(FPKM))
}

et_b = t_flt %>% 
    select(condition, gid, exp.val) %>%
    spread(condition, exp.val)
em_b = as.matrix(et_b[,-1])
rownames(em_b) = et_b$gid

gids = rownames(em_b)
rids = tf_ids[tf_ids %in% gids]
length(rids)
tids = gids
cat(sprintf("%d rids, %d tids\n", length(rids), length(tids)))

weightMat <- GENIE3(em_b, 
                    regulators = rids,
                    treeMethod = "RF",
                    K = 'sqrt',
                    nTrees = 1000,
                    nCores = thread,
                    verbose = T)
dim(weightMat)
#links <- getLinkList(weightMat, reportMax=10000000)
#dim(links)
tn = as_tibble(weightMat) %>%
    mutate(reg.gid = rownames(weightMat)) %>%
    gather(tgt.gid, weight, -reg.gid) %>%
    filter(weight > 0, reg.gid != tgt.gid) %>%
    arrange(desc(weight))

save(rids, tids, tn, file = fo)

