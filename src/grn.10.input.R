#{{{ head
source("grn.fun.r")
require(reticulate)
f_cfg = file.path(dird, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
#th = t_cfg %>% mutate(fgn = sprintf("%s/12_output/%s.rda", dirw, nid)) 
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
    th0 = tl %>% select(-paired) %>%
        replace_na(list(Tissue = '', Genotype = '', Treatment = '')) %>%
        mutate(condition = sprintf("%s.%s.%s", Tissue, Genotype, Treatment))
    th1 = th0 %>% select(SampleID, condition)
    if(!str_detect(mid, 'me[ct]')) {
        tm0 = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
            inner_join(th1, by = 'SampleID') %>%
            group_by(condition, gid) %>%
            summarise(CPM = mean(CPM), FPKM = mean(FPKM)) %>%
            ungroup()
    }
    for (i in 1:nrow(tc1)) {
        nid = tc1$nid[i]; tag = tc1$tag[i]
        if(mid %in% c('me17a', 'me18a') & !is.na(tag)) {
            th = th0 %>% filter(Tissue == tag)
            t_exp = tm0 %>% filter(condition %in% th$condition)
        } else if(mid %in% c('me99b') & !is.na(tag)) {
            th = th0 %>% filter(Genotype == tag)
            t_exp = tm0 %>% filter(condition %in% th$condition)
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

sapply(mids, write_genie3_input, t_cfg = t_cfg, diro = diro)
#}}}


