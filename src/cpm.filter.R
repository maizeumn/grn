#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

p <- ArgumentParser(description = 'Filter expression matrix for network analysis')
p$add_argument("fi", nargs=1, help="Input file with expression matrix (cpm.rds)")
p$add_argument("fo", nargs=1, help="output file (*.tsv)")
p$add_argument("--subid", default='z',
               help="subset ID used to filter expression matrix [default: %(default)s]")
p$add_argument("--num_sam_on", type='integer', default=0,
               help="a gene needs to be expressed in at least XX samples [default: %(default)s]")
p$add_argument("--pct_sam_on", type='double', default=0,
               help="a gene needs to be expressed in at least XX%% of samples [default: %(default)s]")
p$add_argument("--min_cpm", type='double', default=1,
               help="CPM threshold to determine expression on/off [default: '%(default)s']")
p$add_argument("--min_var_p", type='double', default=0,
               help="genes with variance below this percentile will be removed [default: '%(default)s']")
p$add_argument("--no_asinh", action='store_true',
               help="use raw CPM/FPKM values (no asinh transformation) [default: '%(default)s']")
p$add_argument("--use_fpkm", action='store_true',
               help="use FPKM instead of CPM [default: '%(default)s']")
args <- p$parse_args()

fi = args$fi
fo = args$fo
subid = args$subid
if( file.access(fi) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi))

source("~/projects/grn/src/functions.R")
x = readRDS(fi)
tm = x$tm_m

if(subid != 'z') {
    if(subid %in% x$th_m$Tissue) {
        sids = x$th_m %>% filter(Tissue == subid) %>% pull(SampleID)
    } else if(subid %in% x$th_m$Genotype) {
        sids = x$th_m %>% filter(Genotype == subid) %>% pull(SampleID)
    } else {
        cat(sprintf("subid [%s] not found\n", subid))
        return(1)
    }
    tm = tm %>% filter(SampleID %in% sids)
}

if(args$use_fpkm) {
    ti = tm %>% select(gid, SampleID, val=FPKM)
} else {
    ti = tm %>% select(gid, SampleID, val=CPM)
}
cat(sprintf("%d genes in %d samples read\n", length(unique(ti$gid)), length(unique(ti$SampleID))))

if(!args$no_asinh)
    ti = ti %>% mutate(val = asinh(val))

tis = ti %>% group_by(gid) %>%
    summarise(nsam_on = sum(val >= args$min_cpm),
              psam_on = nsam_on/n(),
              val_sd = sd(val)) %>%
    ungroup()
gids = tis %>%
    filter(nsam_on >= args$num_sam_on,
           psam_on >= args$pct_sam_on,
           val_sd >= 0) %>% pull(gid)

tis2 = tis %>% filter(gid %in% gids)
min_sd = quantile(tis2$val_sd, args$min_var_p)
gids = tis2 %>% filter(val_sd >= as.numeric(min_sd)) %>% pull(gid)

to = ti %>% filter(gid %in% gids) %>%
    spread(SampleID, val)

cat(sprintf("%d genes passed filtering\n", nrow(to)))

write_tsv(to, fo)


