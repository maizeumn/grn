#{{{ head
source("grn.fun.r")
f_cfg = file.path(dird, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
#th = t_cfg %>% mutate(fgn = sprintf("%s/12_output/%s.rda", dirw, nid)) 
#}}}

#{{{ prepare genie3 input
diri = '~/projects/maize.expression/data'
#study = 'li2013'
#study = 'hirsch2014'
#study = 'leiboff2015'
#study = 'jin2016'
#study = 'stelpflug2016'
study = 'walley2016'
#study = 'lin2017'
#study = 'kremling2018'
#study = 'briggs'
#study = 'dev41'
#study = 'dev64'
study
diri1 = file.path(diri, study, 'data')
if(study %in% c("dev41","dev64")) diri1 = file.path(diri, 'data', study)
fi = file.path(diri1, '20.rc.norm.rda')
x = load(fi)
fh1 = file.path(diri1, '01.reads.tsv')
fh2 = file.path(diri1, '02.reads.corrected.tsv')
fh = ifelse(file.exists(fh2), fh2, fh1) 
th = read_tsv(fh)
diro = file.path(dird, '11_input')
ngene = 46117
#
if(study == 'li2013') {
    #{{{ li2013
    th = th %>% select(SampleID, Genotype, Replicate) %>%
        mutate(Replicate = sprintf("SAM%d", Replicate))
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = Replicate,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'hirsch2014') {
    #{{{ hirsch2014
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'seedling_503')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'leiboff2015') {
    #{{{ leiboff2015
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'SAM_380')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Genotype, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'jin2016') {
    #{{{
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'kernel_368')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'stelpflug2016') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73_18')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'walley2016') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73_23')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'lin2017') {
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


