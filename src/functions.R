require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
require(ggforce)
require(knitr)
require(kableExtra)
options(knitr.table.format = "latex")
require(pROC)
require(PRROC)
dirg = '~/data/genome'
dirp = '~/projects/grn'
dird = file.path(dirp, 'data')
dirr = file.path(dird, 'raw')
gcfg = read_genome_conf()
#
f_cfg = file.path(dird, '10.dataset.xlsx')
t_cfg = read_xlsx(f_cfg) %>% fill(mid, study)
studies = t_cfg %>% distinct(study) %>% pull(study)
net_types = c("Tissue","Genotype","Tissue*Genotype",'RIL','Meta')
net_cols = pal_aaas()(length(net_types))
names(net_cols) = net_types
nids_meta = c("n17a","n18a",'n18g',"n99a")
t_cfg = t_cfg %>% select(-mid) %>%
    mutate(net_type = factor(net_type, levels = net_types)) %>%
    mutate(col = net_cols[net_type]) %>%
    mutate(lgd = sprintf("%s %s [%d]", str_to_title(study),note,sample_size)) #%>%
    #select(-study,-note,-sample_size)
tsyn = read_syn(gcfg)
gs = readRDS('~/projects/grn/data/09.gs.rds')
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
cols100v = viridis_pal(direction=-1,option='magma')(100)

read_tf_info <- function(fi='~/projects/grn/data/tf_info.xlsx')
    read_xlsx(fi)
read_tfbs_regulations <- function(fi='~/projects/grn/data/04_tfbs/15.regulations.tsv')
    read_tsv(fi)
read_ko <- function(fd = file.path(dird, '07_mutants', 'degs.rds'))
    readRDS(fd)

read_tf_ko_bs <- function(fi = file.path(dird, '13_eval', '05.tf.ko.bs.rds')) {
    #{{{
    ev = readRDS(fi)
    ev1 = ev$tf %>% mutate(eopt='tf') %>%
        gather(key, score, -gopt,-eopt,-nid,-ctag,-gene_alias,-tissue)
    ctags1 = sort(unique(ev$tf$ctag))
    ev2 = ev$ko %>% mutate(eopt='ko') %>% select(-kid) %>%
        gather(key, score, -gopt,-eopt,-nid,-ctag,-gene_alias,-tissue)
    ctags2 = sort(unique(ev$ko$ctag))
    ev3 = ev$bs %>% mutate(eopt='bs') %>% mutate(gene_alias=NA,tissue=NA) %>%
        gather(key, score, -gopt,-eopt,-nid,-ctag,-gene_alias,-tissue)
    ctags3 = as.character(sort(unique(ev$bs$ctag)))
    ctags = c(ctags1, ctags2, ctags3)
    #
    ev = rbind(ev1, ev2, ev3) %>%
        mutate(gopt = str_to_upper(gopt)) %>%
        mutate(score = ifelse(key=='auroc', score*1000, score)) %>%
        mutate(ctag = factor(ctag, levels=ctags)) %>%
        mutate(lab = number(score, accuracy=1)) %>%
        mutate(lab = str_remove(lab, '^0+')) %>%
        mutate(lab = ifelse(key=='spc.pval' & score < -log10(.05), '', lab))
    ev
    #}}}
}
plot_tile <- function(tp, t_cfg, lgd.opt=1, col.opt=1, faceting=F, ytext=T) {
    #{{{
    tp = tp %>%
        inner_join(t_cfg, by = 'nid') %>%
        mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
    tps = t_cfg %>% distinct(lgd, col)
    swit = max(tp$score) / 2
    #
    if(lgd.opt == 1)
        lgd = bquote('AUC'~x10^3~'(FPR<=0.1)')
    else if(lgd.opt == 2)
        lgd = bquote(-log[10]~'(Wilcox P-value)')
    #
    if(col.opt == 1)
        cols = cols100v
    else if(col.opt == 2)
        cols = cols100
    #
    p1 = ggplot(tp, aes(x=ctag, y=lgd, fill=score)) +
        geom_tile() +
        geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
        scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
        scale_y_discrete(drop=F, expand=c(0,0)) +
        scale_fill_gradientn(name=lgd, colors=cols) +
        scale_color_manual(values=c('black','white')) +
        otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
               ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
        guides(color=F)
    if(ytext)
        p1 = p1 + theme(axis.text.y = element_text(color=rev(tps$col), size=7))
    else
        p1 = p1 + theme(axis.text.y = element_blank())
    if(faceting)
        p1 = p1 + facet_grid(.~gopt)
    p1
    #}}}
}

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
run_deseq2 <- function(gene_alias, Tissue, tm, th) {
    #{{{
    require(DESeq2)
    require(edgeR)
    cat(sprintf('--> working on %s - %s\n', gene_alias, Tissue))
    th1 = th %>% filter(gene_alias == !!gene_alias, Tissue == !!Tissue)
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    #{{{ prepare data
    vh = th1 %>% mutate(Genotype = factor(Genotype)) %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~Genotype)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res1 = results(dds, contrast=c("Genotype","WT",'mutant'), pAdjustMethod="fdr")
    stopifnot(rownames(res1) == gids)
    #
    t_ds = tibble(gid = gids, disp = disp,
                padj = res1$padj, log2fc = res1$log2FoldChange,
                ) %>%
        replace_na(list(padj = 1))
    t_ds
    #}}}
}

call_deg_spe <- function(th, tm, tc, tm_m) { # th should have 'group' column
    #{{{
    require(DESeq2)
    require(edgeR)
    th1 = th; tm1 = tm
    #{{{ prepare data
    vh = th1 %>% mutate(Genotype = factor(Genotype)) %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount>=10)) %>%
        filter(n.sam > .2 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    x = readcount_norm(vm)
    mean.lib.size = mean(x$tl$libSize)
    vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~group)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    get_results <- function(group1, group2, dds) {
        res = results(dds, contrast=c("group",group1,group2), pAdjustMethod="fdr")
        as_tibble(res) %>% mutate(gid=rownames(res)) %>%
            select(gid, everything())
    }
    res = tc %>%
        mutate(g1 = str_c(cond, group1, sep="_")) %>%
        mutate(g2 = str_c(cond, group2, sep="_")) %>%
        mutate(data=map2(g1, g2, get_results, dds=dds)) %>%
        select(-g1, -g2) %>% unnest() %>%
        select(-baseMean,l2fc=log2FoldChange,-stat,-lfcSE,-pvalue)
    t_ds = res %>% replace_na(list(padj = 1)) %>%
        inner_join(tm_m, by=c('cond'='cond','group1'='group','gid'='gid')) %>%
        rename(cpm1=CPM) %>%
        inner_join(tm_m, by=c('cond'='cond','group2'='group','gid'='gid')) %>%
        rename(cpm2=CPM)
    DEtags = c("non_DE",'DE1-2','DE2-4','DE4+','SPE')
    dirtags = c('up','dn')
    t_ds %>%
        mutate(DE=ifelse(padj<.01, 'DE1-2', 'non_DE')) %>%
        mutate(DE=ifelse(padj<.01 & abs(l2fc)>=1, 'DE2-4', DE)) %>%
        mutate(DE=ifelse(padj<.01 & abs(l2fc)>=2, 'DE4+', DE)) %>%
        mutate(DE=ifelse(padj<.01 & ((cpm1<=.1 & cpm2 >=1) | (cpm1>=1 & cpm2<=.1)), 'SPE', DE)) %>%
        mutate(DEdir=ifelse(l2fc < 0, 'up', 'dn')) %>%
        mutate(DE = factor(DE, levels=DEtags)) %>%
        mutate(DEdir = factor(DEdir, levels=dirtags))
    #}}}
}

sum_degs <- function(t_ds) {
    #{{{
    t_ds %>% group_by(yid, cond, group1, group2) %>%
        summarise(ntot = n(),
                  deg_1p = sum(padj<.01 & lfc>0),
                  deg_1n = sum(padj<.01 & lfc<0),
                  deg_2p = sum(padj<.05 & lfc>1),
                  deg_2n = sum(padj<.05 & lfc< -1)) %>%
        ungroup() %>% print(n=50)
    #}}}
}



