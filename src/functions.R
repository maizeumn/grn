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
symb = read_symbol() %>% mutate(symbol=ifelse(gid=='Zm00001d026147','R1',symbol))
gs = readRDS('~/projects/grn/data/09.gs.rds')
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
cols100v = viridis_pal(direction=-1,option='magma')(100)

trim_gid <- function(gid, opt=1) str_replace(str_replace(gid,'^Zm00001[de]',''), '^GRMZM2G','')
read_tf_info <- function(fi='~/projects/grn/data/tf_info.xlsx')
    read_xlsx(fi)
read_tfbs_regulations <- function(fi='~/projects/grn/data/04_tfbs/15.regulations.tsv')
    read_tsv(fi)
read_ko <- function(fd = file.path(dird, '07_known_tf', 'degs.rds'))
    readRDS(fd)

read_eval_ko <- function(gopts=c("et","rf",'xgb')) {
    #{{{
    nrow_positive <- function(ti) sum(ti$response == 1)
    ko = read_ko() %>% filter(!tf %in% c('RA3')) %>%
        mutate(kid=1:n()) %>%
        select(kid,tf,tissue,reg.gid,ds) %>%
        unnest(ds) %>%
        mutate(response=ifelse(padj < .01, 1, 0)) %>%
        rename(tgt.gid = gid) %>%
        group_by(kid,tf,tissue) %>% nest() %>% rename(res = data) %>%
        mutate(n_tot = map_dbl(res, nrow)) %>%
        mutate(n_de = map_dbl(res, nrow_positive)) %>%
        mutate(prop_de=n_de/n_tot) %>%
        mutate(ctag=sprintf("%s_%s [%s] [%s]", tf, tissue, number(n_de), percent(prop_de, accuracy=.1))) %>%
        select(kid, tf, tissue, ctag, n_de)
    ctags = ko$ctag
    #
    ev = tibble(gopt = gopts) %>%
        mutate(fi = sprintf('%s/raw/%s.ko.rds', dird, gopt)) %>%
        mutate(data = map(fi, readRDS)) %>% select(-fi) %>% unnest(data) %>%
        select(gopt,nid,net_size,tf,tissue,auroc,pval,auroc0,auprc,tp,fp) %>%
        inner_join(ko, by=c("tf", "tissue")) %>%
        mutate(gopt = str_to_upper(gopt)) %>%
        mutate(ctag = factor(ctag, levels=ctags)) %>%
        group_by(gopt, net_size, ctag, n_de) %>%
        mutate(padj = p.adjust(pval, method='BH')) %>%
        ungroup() %>%
        mutate(score1 = auroc * 1000) %>%
        mutate(lab1 = str_remove(number(score1, accuracy=1), '^0+')) %>%
        mutate(score2 = -log10(padj)) %>%
        mutate(lab2 = str_remove(number(score2, accuracy=1), '^0+')) %>%
        mutate(lab2 = ifelse(padj < .05, lab2, '')) %>%
        mutate(score3 = auroc0) %>%
        mutate(lab3 = str_remove(number(score3, accuracy=.01), '^0+')) %>%
        mutate(score4 = auprc) %>%
        mutate(lab4 = str_remove(number(score4, accuracy=.01), '^0+')) %>%
        mutate(score5 = tp / n_de) %>%
        mutate(lab5 = number(tp, big.mark=','))
    ev
    #}}}
}
read_eval_bs <- function(gopts=c("et",'rf','xgb')) {
    #{{{
    tf = read_tf_info() %>% distinct(tf, gid) %>% rename(reg.gid=gid, name=tf)
    bs = gs$bs %>% count(ctag, reg.gid) %>%
        left_join(tf, by='reg.gid') %>%
        rename(tf = reg.gid) %>%
        mutate(name = ifelse(str_detect(ctag, '^REF'), str_replace(ctag, 'REF\\|', ''), name)) %>%
        mutate(name = ifelse(str_detect(ctag, '^P1'), ctag, name)) %>%
        mutate(xlab = sprintf("%s (%s)", name, number(n))) %>%
        select(ctag, tf, n, xlab)
    #
    ev = tibble(gopt = gopts) %>%
        mutate(fi = sprintf('%s/raw/%s.bs.rds', dird, gopt)) %>%
        mutate(data = map(fi, readRDS)) %>% select(-fi) %>% unnest(data) %>%
        select(gopt, nid, net_size, ctag, tf, auroc, pval, auroc0, auprc,tp,fp) %>%
        inner_join(bs, by=c("ctag",'tf')) %>%
        mutate(gopt = str_to_upper(gopt)) %>%
        group_by(gopt, net_size, ctag, tf, n) %>%
        mutate(padj = p.adjust(pval, method='BH')) %>%
        ungroup() %>%
        mutate(score1 = auroc * 1000) %>%
        mutate(lab1 = str_remove(number(score1, accuracy=1), '^0+')) %>%
        mutate(score2 = -log10(padj)) %>%
        mutate(lab2 = str_remove(number(score2, accuracy=1), '^0+')) %>%
        mutate(lab2 = ifelse(padj < .05, lab2, '')) %>%
        mutate(score3 = auroc0) %>%
        mutate(lab3 = str_remove(number(score3, accuracy=.01), '^0+')) %>%
        mutate(score4 = auprc) %>%
        mutate(lab4 = str_remove(number(score4, accuracy=.01), '^0+')) %>%
        mutate(score5 = tp / n) %>%
        mutate(lab5 = number(tp, big.mark=','))
    ev
    #}}}
}
read_eval_go <- function(gopts=c('et','rf','xgb')) {
    #{{{
    ev = tibble(gopt = gopts) %>%
        mutate(fi = sprintf("%s/%s.go.rds", dirr, gopt)) %>%
        mutate(data = map(fi, readRDS)) %>% select(-fi) %>% unnest(data)
    ev
    #}}}
}
plot_tile <- function(tp, t_cfg, lgd.opt=1, col.opt=1, faceting=F, ytext=T) {
    #{{{
    tp = tp %>% filter(!is.na(score)) %>%
        inner_join(t_cfg, by = 'nid') %>%
        mutate(lgd = factor(lgd, levels=rev(t_cfg$lgd)))
    tps = t_cfg %>% distinct(lgd, col)
    swit = (min(tp$score) + max(tp$score)) / 2
    #
    if(lgd.opt == 1)
        lgd = bquote('AUC'~x10^3~'(FPR<=0.1)')
    else if(lgd.opt == 2)
        lgd = bquote(-log[10]~'(Wilcox P-value adjusted)')
    else if(lgd.opt == 3)
        lgd = 'AUROC'
    else if(lgd.opt == 4)
        lgd = 'AUPRC'
    else if(lgd.opt == 5)
        lgd = 'Prop. True Targets'
    #
    if(col.opt == 1)
        cols = cols100v
    else if(col.opt == 2)
        cols = cols100
    #
    p1 = ggplot(tp, aes(x=xlab, y=lgd, fill=score)) +
        geom_tile() +
        geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
        scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
        scale_y_discrete(drop=F, expand=c(0,0)) +
        scale_fill_gradientn(name=lgd, colors=cols) +
        scale_color_manual(values=c('black','white')) +
        otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
               ygrid=T, xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T) +
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=7)) +
        theme(legend.title=element_text(size=8)) +
        theme(legend.text=element_text(size=7)) +
        #theme(legend.key.size = unit(.5, 'lines')) +
        guides(color=F, fill = guide_colourbar(barheight=.6))
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
get_ds <- function(condR, cond, dds, gids) {
    #{{{
    res1 = results(dds, contrast = c("cond",condR,cond), pAdjustMethod='fdr')
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}
run_deseq2 <- function(th, tm) {
    #{{{
    #cat(sprintf('--> working on %s - %s\n', gene_alias, Tissue))
    th1 = th %>% mutate(Genotype=str_replace(Genotype, 'WT', 'wt')) %>%
        mutate(cond = str_c(Tissue, Genotype, sep="."))
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    ct = th1 %>% distinct(Tissue,Genotype,cond) %>% filter(Genotype != 'wt') %>%
        mutate(condR = str_c(Tissue, 'wt', sep="."))
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
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = ct %>% mutate(ds = map2(condR, cond, get_ds, dds = dds, gids = gids))
    res
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



