#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Merge GRN model eval output')
parser$add_argument("fi", nargs='+', help = "model eval file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="model.eval.tsv",
                    help = "output file [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(tidyverse)

tv = tibble(fi=fis) %>%
    mutate(fname = basename(fi)) %>%
    mutate(nid = str_replace(fname, "\\.[a-z]+$", "")) %>%
    mutate(data = map(fi, read_tsv, col_names=c('gid','nid_b','score'))) %>%
	select(nid, data) %>% unnest()

saveRDS(tv, file=fo)

