source("functions.R")
require(reticulate)
fi = file.path(dird, '09.gs.rds')
gs = readRDS(fi)
tf_ids = gs$tf_ids

write_me <- function(mid, t_cfg, diro, diri='~/projects/rnaseq/data', use_cpm=T) {
    #{{{
    fh1 = sprintf("%s/05_read_list/%s.tsv", diri, mid)
    fh2 = sprintf("%s/05_read_list/%s.c.tsv", diri, mid)
    fh = ifelse(file.exists(fh2), fh2, fh1)
    stopifnot(file.exists(fh))
    th = read_tsv(fh)
    #
    fi = sprintf("%s/08_raw_output/%s/cpm.rds", diri, mid)
    stopifnot(file.exists(fi))
    res = readRDS(fi)
    tm = res$tm
    #
    if(mid %in% c('me99a','me99c')) {
        ths = th %>% distinct(Tissue, Genotype, Treatment, inbred) %>%
            mutate(nSampleID = sprintf("%s_%d", mid, 1:length(Tissue)))
        th = th %>% inner_join(ths, by = c("Tissue","Genotype","Treatment",'inbred'))
    } else {
        ths = th %>% distinct(Tissue, Genotype, Treatment) %>%
            mutate(nSampleID = sprintf("%s_%d", mid, 1:length(Tissue)))
        th = th %>% inner_join(ths, by = c("Tissue","Genotype","Treatment"))
    }
    t_map = th %>% select(SampleID, nSampleID)
    th = ths %>% mutate(SampleID=nSampleID) %>% select(-nSampleID)
    #
    tm = tm %>% inner_join(t_map, by = 'SampleID') %>%
        mutate(SampleID = nSampleID) %>%
        group_by(gid, SampleID) %>%
        summarise(CPM = mean(CPM), FPKM = mean(FPKM)) %>%
        ungroup()
    #
    tc = t_cfg %>% filter(mid == !!mid)
    for (i in 1:nrow(tc)) {
        nid = tc$nid[i]; tag = tc$note[i]
        #{{{ get th1 and tm1
        if(mid %in% c('me17a','me18a','me99c') &
           ! nid %in% c('n17a','n18a','n99c')) {
            th1 = th %>% filter(tolower(Tissue) == tolower(tag))
        } else if(mid %in% c('me99b') & !nid %in% c('n99b')) {
            th$Genotype[th$Genotype=='B73xMo17'] = 'BxM'
            th1 = th %>% filter(Genotype == tag)
        } else if(mid == 'me99a' & str_detect(nid, "_[1-5]$")) {
            th1 = th %>% filter(Tissue == tag)
        } else if(mid == 'me99a' & nid == 'n99a_6') {
            th1 = th %>% filter(Tissue == 'seedling', inbred == T)
        } else if(mid == 'me99a' & nid == 'n99a_7') {
            th1 = th %>% filter(Tissue == 'seedling', inbred == F)
        } else {
            th1 = th
        }
        tm1 = tm %>% filter(SampleID %in% th1$SampleID)
        #}}}
        cat(sprintf('%s %d\n', nid, nrow(th1)))
        fo = sprintf("%s/%s", diro, nid)
        write_me1(tm1, th1, tf_ids, fo, use_cpm)
    }
    T
    #}}}
}
write_me1 <- function(t_exp, th, tf_ids, fo, use_cpm = T) {
    #{{{
    if(use_cpm) {
        t_exp = t_exp %>% mutate(exp.val = asinh(CPM))
    } else {
        t_exp = t_exp %>% mutate(exp.val = asinh(FPKM))
    }
    et = t_exp %>%
        select(SampleID, gid, exp.val) %>%
        spread(SampleID, exp.val)
    gids = et %>% pull(gid)
    fo1 = sprintf("%s.tsv", fo)
    em = as_tibble(t(et[,-1]))
    write_tsv(em, fo1, col_names = F)
    fo2 = sprintf("%s.gid.txt", fo)
    write(gids, fo2)
    #}}}
}
write_me2 <- function(t_exp, th, tf_ids, fo, use_cpm = T) {
    #{{{
    if(use_cpm) {
        t_exp = t_exp %>% mutate(exp.val = asinh(CPM))
    } else {
        t_exp = t_exp %>% mutate(exp.val = asinh(FPKM))
    }
    nsam = length(unique(t_exp$SampleID))
    gids = t_exp %>% group_by(gid) %>%
        summarise(nsam_exp = sum(CPM >= 1)) %>%
        filter(nsam_exp >= nsam * .1) %>%
        pull(gid)
    t_flt = t_exp %>% filter(gid %in% gids)
    et = t_exp %>%
        select(SampleID, gid, exp.val) %>%
        spread(SampleID, exp.val)
    et_f = t_flt %>%
        select(SampleID, gid, exp.val) %>%
        spread(SampleID, exp.val)
    tids = et %>% pull(gid)
    tids_f = et_f %>% pull(gid)
    em = t(as.matrix(et[,-1]))
    em_f = t(as.matrix(et_f[,-1]))
    rids = unique(tf_ids[tf_ids %in% gids])
    fo = normalizePath(fo, mustWork = F)
    x = list(em_f, rids, tids_f, tids, em)
    py_save_object(x, fo)
    #}}}
}

diro = file.path(dird, '11_exp_mat')
mids = t_cfg %>% filter(str_detect(mid,'^me')) %>% distinct(mid) %>% pull(mid)
#mids = c('me17a','me18a',sprintf('me99%s',letters[1:3]))
map_int(mids, write_me, t_cfg = t_cfg, diro = diro)


