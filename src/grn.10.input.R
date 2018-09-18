#{{{ head
source("grn.fun.r")
f_cfg = file.path(dird, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
#th = t_cfg %>% mutate(fgn = sprintf("%s/12_output/%s.rda", dirw, nid)) 
#}}}

#{{{ prepare GRN input
diri = '~/projects/maize.expression/data/09_output'
mids = t_cfg %>% distinct(mid) %>% pull(mid)
for (mid in mids) {
    tc1 = t_cfg %>% filter(mid == !!mid)
    fi = sprintf("%s/%s.rda", diri, mid)
    stopifnot(file.exists(fi))
    x = load(fi)
    x
    th = tl %>% select(-paired) %>%
        replace_na(list(Tissue = '', Genotype = '', Treatment = '')) %>%
        mutate(condition = sprintf("%s.%s.%s", Tissue, Genotype, Treatment)) 
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(condition, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>%
        ungroup() 
}

    t_exp2 = t_exp %>% 
        mutate(condition = sprintf("%s_%s", tag, condition),
               tag = '5_tissues')
    t_exp = t_exp %>% bind_rows(t_exp2)
    #}}}
} else if(study == 'kremling2018') {
    #{{{
    th = th %>% select(SampleID, Genotype, Tissue) %>%
        mutate(tag = Tissue)
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Genotype, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    t_exp2 = t_exp %>% 
        mutate(condition = sprintf("%s_%s", tag, condition),
               tag = '7_tissues')
    t_exp = t_exp %>% bind_rows(t_exp2)
    #}}}
} else if(study == 'briggs') {
    #{{{
    th = th %>% select(SampleID, Genotype, Tissue) %>%
        mutate(tag = Genotype)
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'dev41' | study == 'dev64') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else {
    cat('unknown study', study, "\n")
}
nrow(t_exp)/ngene
tt = t_exp
#
tags = unique(tt$tag)
for (tag in tags) {
    tc1 = t_cfg %>% filter(study == !!study, tag == !!tag)
    if(nrow(tc1) != 1)
        stop(sprintf("not 1 row: %s %s\n", study, tag))
    nid = tc1$nid[1]
    t_exp = tt %>% filter(tag == !!tag) %>% select(-tag)
    fo = sprintf("%s/%s.rda", diro, nid)
    save(t_exp, file = fo)
}
#}}}


