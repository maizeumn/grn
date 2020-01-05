source("functions.R")
dirw = file.path(dird, '03_tfbs_chipseq')
tf = read_tf_info() %>% mutate(name=ifelse(is.na(name), tf, name))


#{{{ merge NF chipseq/dapseq results
read_chipseq_peaks <- function(yid, min_frip=.05, tf_info=tf, diri = '~/projects/nf/data/out') {
    #{{{
    fi = sprintf("%s/%s.rds", diri, yid)
    res = readRDS(fi)
    grps = res$frip %>% filter(frip >= min_frip) %>% select(-frip)
    to = res$peaks %>% inner_join(grps, by=c('group','rep')) %>%
        select(group, rep, chrom, start, end, srd, score,tgt.gid=gid,dist_tss)
    tf1 = tf %>% filter(yid==!!yid) %>% select(tf, reg.gid=gid, ctag=reference)
    if(yid == 'ca19a4')
        to1 = to %>% rename(tf = group, note = rep)
    else if (yid == 'cp12b')
        to1 = to %>% mutate(group = str_to_upper(group)) %>%
            rename(tf=group, note=rep)
    else if (yid == 'cp14g')
        to1 = to %>% mutate(tf = 'RA1', note=str_replace(group, '_control','')) %>%
            select(-group,-rep)
    else if (yid == 'cp15a')
        to1 = to %>% mutate(tf='FEA4') %>% select(-group) %>% rename(note=rep)
    else if (yid == 'cp15b' | yid == 'cp18a')
        to1 = to %>% mutate(tf='O2') %>% select(-group) %>% rename(note=rep)
    to1 %>% inner_join(tf1, by='tf') %>%
            select(ctag, tf, note, reg.gid, everything())
    #}}}
}

to1 = read_chipseq_peaks('ca19a4')
to2 = read_chipseq_peaks('cp12b')
to3 = read_chipseq_peaks('cp14g')
to4 = read_chipseq_peaks('cp15a')
to5 = read_chipseq_peaks('cp18a')

to = rbind(to1, to2, to3, to4, to5)
fo = file.path(dirw, "tfbs.nf.tsv")
write_tsv(to, fo)
#}}}


#{{{ Ricci2019 DAP-Seq 27 TFs
fi = file.path(dirw, 'Ricci2019.dapseq.txt')
ti = read_tsv(fi, col_names='fname') %>%
    mutate(fi = sprintf("%s/Ricci2019/%s", dirw, fname)) %>%
    mutate(tf = str_replace(fname,'\\.dap\\.bed.gz','')) %>%
    mutate(tf = str_replace(tf,'^GSM.*DAP_','')) %>%
    inner_join(tf, by='tf') %>% select(-tf) %>% rename(tf=gid)

to = ti %>%
    mutate(data=map(fi, read_tsv, col_names=c('chrom','start','end'))) %>%
    select(tf,data) %>% unnest() %>%
    select(chrom,start,end,tf)
fo = file.path(dirw, 'Ricci2019_tmp/01.raw.bed')
write_tsv(to, fo, col_names=F)
#}}}


