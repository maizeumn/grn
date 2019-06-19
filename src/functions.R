require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
dirg = '~/data/genome'
dirp = '~/projects/grn'
dird = file.path(dirp, 'data')
dirr = file.path(dird, 'raw')
gcfg = read_genome_conf()
#
f_cfg = file.path(dird, '10.dataset.xlsx')
t_cfg = read_xlsx(f_cfg) %>% fill(mid, study)
studies = t_cfg %>% distinct(study) %>% pull(study)
net_types = c("tissue","genotype","tissue*genotype",'ril','liftover')
net_cols = pal_aaas()(length(net_types))
names(net_cols) = net_types
nids_meta = c("n17a","n18a","n99a","n99b","n99c")
th = t_cfg %>% select(-mid) %>%
    filter(net_type %in% net_types) %>%
    mutate(net_type = ifelse(nid %in% nids_meta, 'genotype', net_type)) %>%
    mutate(net_type = factor(net_type, levels = net_types)) %>%
    arrange(net_type, nid, sample_size) %>%
    mutate(net_type = as.character(net_type)) %>%
    mutate(net_type = ifelse(nid %in% nids_meta, 'tissue*genotype', net_type)) %>%
    mutate(net_type = factor(net_type, levels = net_types)) %>%
    mutate(col = net_cols[net_type]) %>%
    mutate(lgd = sprintf("%s %s [%d]", study,note,sample_size)) #%>%
    #select(-study,-note,-sample_size)
th1 = th %>% filter(!str_detect(nid, "^n99[a]")) %>%
    select(nid, net_type, sample_size, col, lgd)
th2 = th %>% filter(!str_detect(nid, "^n99[a]")) %>%
    filter(!str_detect(nid, "^n18g")) %>%
    select(nid, net_type, sample_size, col, lgd)
nids_geno = th %>% filter(net_type == 'genotype') %>% pull(nid)
nids_dev = th %>% filter(net_type == 'tissue') %>% pull(nid)
nids22 = t_cfg %>%
    filter(!str_detect(nid, '^n((17a)|(18a)|(99a)|(99b)|(99c))_')) %>%
    pull(nid)
nids25 = t_cfg %>%
    filter(!str_detect(nid, '^n((17a)|(18a)|(99a)|(99c))_')) %>% pull(nid)
tsyn = read_syn(gcfg)
gs = readRDS('~/projects/grn/data/09.gs.rds')

read_briggs <- function(fi="~/projects/briggs/data/49_coop/01.master.rda") {
    #{{{ read briggs data
    x = load(fi)
    tissues23 = c("seedlingleaf_11DAS", "blade_v12", "flagleaf_0DAP",
            "husk_0DAP", "sheath_v12", "auricle_v12",
            "floret_0DAP", "tasselstem_0DAP",
            "internode_v12",
            "root_0DAP", "seedlingroot_11DAS", "radicle_root",
            "coleoptile_tip",
            "silk_0DAP",
            "tassel_v12", "spikelets_0DAP", "ear_v14",
            "seedlingmeristem_11DAS",
            "embryo_27DAP", "embryo_imbibedseed",
            "endosperm_27DAP", "kernel_14DAP", "endosperm_14DAP")
    de = tm %>% filter(silent == 0) %>%
        mutate(DE= ifelse(pDE == 'non_DE', 'non_DE',
            ifelse(abs(log2mb) < 1, 'DE1-2',
            ifelse(abs(log2mb) < 2, 'DE2-4',
            ifelse(SPE == 'non_SPE', 'DE4+', 'SPE'))))) %>%
        mutate(DEdir = ifelse(log2mb<0, 'B>M', 'B<M')) %>%
        mutate(Tissue = as.character(Tissue)) %>%
        select(Tissue, gid, DE, DEdir)
    de %>% count(Tissue, DE, DEdir)
    des = de %>% group_by(Tissue) %>%
        summarise(propDE = sum(DE!='non_DE')/n()) %>% ungroup()
    list(tissues=tissues23, de=de, des=des)
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



