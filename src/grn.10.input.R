#{{{ head
source("grn.fun.r")
require(reticulate)
f_cfg = file.path(dird, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
fi = file.path(dird, '09.gs.rda')
x = load(fi)
x
#}}}

#{{{ prepare GRN input
write_genie3_input <- function(mid, t_cfg, diro, use_cpm = T) {
    #{{{
    cat(mid, '\n')
    tc1 = t_cfg %>% filter(mid == !!mid)
    fi = sprintf("%s/%s.rda", diri, mid)
    stopifnot(file.exists(fi))
    x = load(fi)
    x
    th0 = th; tm0 = tm
    for (i in 1:nrow(tc1)) {
        nid = tc1$nid[i]; tag = tc1$tag[i]
        if(mid %in% c('me17a', 'me18a') & !is.na(tag)) {
            th = th0 %>% filter(Tissue == tag)
            t_exp = tm0 %>% filter(SampleID %in% th$SampleID)
        } else if(mid %in% c('me99b') & !is.na(tag)) {
            th = th0 %>% filter(Genotype == tag)
            t_exp = tm0 %>% filter(SampleID %in% th$SampleID)
        } else {
            th = th0
            t_exp = tm0
        }
        fo = sprintf("%s/%s.pkl", diro, nid)
        write_genie3_input1(t_exp, th, tf_ids, fo, use_cpm)
    }
    #}}}
}
write_genie3_input1 <- function(t_exp, th, tf_ids, fo, use_cpm = T) {
    #{{{
    nsam = length(unique(t_exp$SampleID))
    gids = t_exp %>% group_by(gid) %>%
        summarise(nsam_exp = sum(CPM >= 1)) %>%
        filter(nsam_exp >= nsam * .1) %>%
        pull(gid)
    t_flt = t_exp %>% filter(gid %in% gids)
    if(use_cpm) {
        t_flt = t_flt %>% mutate(exp.val = asinh(CPM))
    } else {
        t_flt = t_flt %>% mutate(exp.val = asinh(FPKM))
    }
    et_b = t_flt %>% 
        select(SampleID, gid, exp.val) %>%
        spread(SampleID, exp.val)
    tids = et_b %>% pull(gid)
    em_b = t(as.matrix(et_b[,-1]))
    rids = tf_ids[tf_ids %in% gids]
    fo = normalizePath(fo, mustWork = F)
    x = list('em' = em_b, 'gene_names' = tids, 'regulators' = rids)
    x = list(em_b, tids, rids)
    py_save_object(x, fo) 
    #}}}
}

diri = '~/projects/maize.expression/data/15_output'
diro = file.path(dird, '11_input')
mids = t_cfg %>% distinct(mid) %>% pull(mid)
mids = c("me13c")
#sapply(mids, write_genie3_input, t_cfg = t_cfg, diro = diro)

#diro = file.path(dird, '11_input_fpkm')
#sapply(mids, write_genie3_input, t_cfg = t_cfg, diro = diro, use_cpm = F)
#}}}


