#{{{ head
source("grn.fun.r")
#source("enrich.R")
fi = file.path(dird, '05.previous.grns/10.RData')
x = load(fi)
fi = file.path(dird, '07.known.tf/10.RData')
x = load(fi)
# read RNA-Seq data
dirw = file.path("~/projects/briggs/data", "49.coop")
fi = file.path(dirw, "01.master.RData")
x = load(fi)
fd = file.path(dirw, "03.sharing.RData")
x = load(fd)
# TFs
ff = '~/data/genome/Zmays_v4/TF/11.tsv'
tf = read_tsv(ff)
tf_ids = tf$gid
#}}}


#{{{ 
require(doParallel)
require(doRNG)
fi = file.path('~/projects/briggs/data', "41.qc/10.rep.merged.RData")
x = load(fi)
nrow(t_rep_merged)/69

et_b = t_rep_merged %>% filter(Genotype == 'B73') %>%
    select(Tissue, gid, CPM) %>%
    spread(Tissue, CPM)
em_b = as.matrix(et_b[,-1])
rownames(em_b) = et_b$gid

gids = rownames(em_b)
rids = tf_ids[tf_ids %in% gids]

weightMat <- GENIE3(em_b, 
                    regulators = rids,
                    treeMethod = "RF", 
                    K=7,
                    nTrees=50,
                    nCores = 4)
dim(weightMat)
links <- getLinkList(weightMat, reportMax=1000000)
dim(links)

fo = file.path(dirw, '10.genie3.RData')
save(weightMat, links, file = fo)
#}}}

#{{{

#}}}
