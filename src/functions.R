require(rmaize)
dirg = '~/data/genome'
dirp = '~/projects/grn'
dird = file.path(dirp, 'data')
#dirr = file.path(dirp, 'Rmd')
#
f_cfg = file.path(dird, '10.dataset.xlsx')
t_cfg = read_xlsx(f_cfg) %>% fill(mid, study)
studies = t_cfg %>% distinct(study) %>% pull(study)
net_types = c("tissue","genotype","tissue*genotype",'ril','liftover')
net_cols = pal_aaas()(length(net_types))
names(net_cols) = net_types
th = t_cfg %>% select(-mid) %>%
    filter(net_type %in% net_types) %>%
    mutate(net_type = factor(net_type, levels = net_types)) %>%
    arrange(net_type, nid, sample_size) %>%
    mutate(col = net_cols[net_type]) %>%
    mutate(txt = sprintf("%s %s [%d]", study,note,sample_size))
nids_geno = th %>% filter(net_type == 'genotype') %>% pull(nid)
nids_dev = th %>% filter(net_type == 'tissue') %>% pull(nid)

load_maize_dataset <- function(id = 'me99b', opt = "exp") {
    #{{{
    if(opt == 'exp') {
        diri = '~/projects/rnaseq/data/08_raw_output'
        fi = sprintf("%s/%s.rda", diri, id)
        stopifnot(file.exists(fi))
        x = load(fi)
        env(th = th, tm = tm)
    } else if (opt == 'grn') {
        diri = '~/projects/grn/data/12_output'
        fi = sprintf("%s/%s.rda", diri, id)
        stopifnot(file.exists(fi))
        x = load(fi)
        env(rids = rids, tids = tids, reg.mat = reg.mat, tn = tn)
    } else {
        stop(sprintf("unknown dataset opt: %s", opt))
    }
    #}}}
}
read_gs <- function(fi='~/projects/grn/data/09.gs.rds') {
    readRDS(fi)
}
read_briggs <- function(fi="~/projects/briggs/data/49_coop/01.master.rda") {
    #{{{ read briggs data
    x = load(fi)
    tissues = unique(tm$Tissue)
    de = tm %>% filter(silent == 0) %>%
        mutate(DE= ifelse(pDE == 'non_DE', 'non_DE',
            ifelse(abs(log2mb) < 1, 'DE1-2',
            ifelse(abs(log2mb) < 2, 'DE2-4',
            ifelse(SPE == 'non_SPE', 'DE4+', 'SPE'))))) %>%
        mutate(DEdir = ifelse(log2mb<0, 'B>M', 'B<M')) %>%
        select(Tissue, gid, DE, DEdir)
    de %>% count(Tissue, DE, DEdir)
    des = de %>% group_by(Tissue) %>%
        summarise(propDE = sum(DE!='non_DE')/n()) %>% ungroup()
    list(tissues=tissues, de=de, des=des)
    #}}}
}
read_biomap <- function(opt='all') {
    #{{{ read biomap data
    fi = '~/projects/biomap/data/41_qc/10.rc.ase.rds'
    bm = readRDS(fi)
    if (opt == 'inbred') {
        sids = bm$th %>%
            separate(Genotype, c('pa1','pa2'), sep='x', fill='right', extra='merge') %>%
            filter(is.na(pa2)) %>% pull(SampleID)
        bm$th = bm$th %>% filter(SampleID %in% sids)
        bm$tm = bm$tm %>% filter(SampleID %in% sids)
        bm = bm$tm %>% inner_join(bm$th, by = 'SampleID') %>%
            select(Tissue,Genotype,gid,CPM)
    }
    bm
    #}}}
}
radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}



