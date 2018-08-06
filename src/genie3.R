#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Run GENIE3')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                     help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                     dest="verbose", help="Print little output")
parser$add_argument("fi", nargs=1, help="Input (expression matrix) file")
parser$add_argument("fo", nargs=1, help="Output file")
parser$add_argument("--tf", default="~/data/genome/Zmays_v4/TF/11.tsv",
                    help = "TF IDs file [default %(default)s]")
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

et_b = t_exp %>% spread(condition, CPM)
em_b = as.matrix(et_b[,-1])
rownames(em_b) = et_b$gid

gids = rownames(em_b)
rids = tf_ids[tf_ids %in% gids]
length(rids)
tids = gids

weightMat <- GENIE3(em_b, 
                    regulators = rids,
                    treeMethod = "RF",
                    K = 'sqrt',
                    nTrees = 1000,
                    nCores = thread,
                    verbose = T)
dim(weightMat)
links <- getLinkList(weightMat, reportMax=1000000)
dim(links)
tn = as_tibble(links) %>% 
    transmute(reg.gid = regulatoryGene,
              tgt.gid = targetGene,
              weight = weight)

save(rids, tids, tn, file = fo)

