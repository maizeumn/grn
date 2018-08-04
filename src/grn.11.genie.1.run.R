#!/usr/bin/env Rscript
require(tidyverse)
require(doParallel)
require(doRNG)
require(GENIE3)

#{{{ head
diri = '~/projects/maize.expression'
dirw = '~/projects/maize.grn/analysis/11_genie3'
# read TFs
ff = '~/data/genome/Zmays_v4/TF/11.tsv'
tf = read_tsv(ff)
tf_ids = tf$gid
#}}}

#{{{ 
study = 'briggs'
fi = file.path(diri, study, 'data/20.rc.norm.RData')
x = load(fi)
fh = file.path(diri, study, 'data/01.reads.tsv')
th = read_tsv(fh)
th = th %>% filter(SampleID %in% tl$SampleID) %>%
    select(SampleID, Tissue, Genotype)

tm2 = tm %>% select(gid, SampleID, CPM) %>% 
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, gid) %>%
    summarise(CPM = mean(CPM)) %>% 
    ungroup()
gts = c("B73", "Mo17", "B73xMo17")

tn = tibble()
for (gt in gts) {
    cat(gt, "\n")
    et_b = tm2 %>% filter(Genotype == gt) %>% select(-Genotype) %>%
        spread(Tissue, CPM)
    em_b = as.matrix(et_b[,-1])
    rownames(em_b) = et_b$gid

    gids = rownames(em_b)
    rids = tf_ids[tf_ids %in% gids]
    length(rids)

    weightMat <- GENIE3(em_b, 
                        regulators = rids,
                        treeMethod = "RF",
                        K = 'sqrt',
                        nTrees = 1000,
                        nCores = 8)
    dim(weightMat)
    links <- getLinkList(weightMat, reportMax=1000000)
    dim(links)
    tn1 = as_tibble(links) %>% 
        transmute(ctag = study,
                  tag = gt,
                  reg.gid = regulatoryGene,
                  tgt.gid = targetGene,
                  weight = weight)
    tn = rbind(tn, tn1)
}

fo = sprintf("%s/%s.rda", dirw, study)
save(tn, file = fo)
#}}}

#{{{

#}}}
