source("functions.R")
ncfg = th
dirw = file.path(dird, '14_eval_sum')
diri = '~/projects/rnaseq'
gcfg = read_genome_conf()

#{{{ take top50k edges from each GRN and save
read_grn <- function(nid='n13a',net_size=5e4,diri='~/projects/grn/data/12_output') {
    #{{{
    cat(nid,'\n')
    f_net = sprintf("%s/%s.rda", diri, nid)
    x = load(f_net)
    tn %>% filter(row_number() <= net_size)
    #}}}
}
t_net = ncfg %>%
    select(nid) %>%
    mutate(data=map(nid, read_grn)) %>% unnest()
fo = file.path(dirw, '01.top50k.rds')
saveRDS(t_net, file=fo)
#}}}

fi = file.path(dirw, '01.top50k.rds')
t_net = readRDS(fi)

fi = file.path(dirw, '01.tf.rds')
ev_tf = readRDS(fi)

#{{{ general stats
tp0 = ev_tf %>%
    select(nid, nstat) %>% unnest() %>%
    filter(net_size == 5e4) %>% select(-net_size)
#{{{ topology stats
tp1 = tp0 %>%
    mutate(deg.reg.q50 = map_dbl(deg.reg, .f <- function(x) median(x$n))) %>%
    select(nid, n.reg, n.tgt, deg.reg.q50) %>%
    inner_join(th, by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt))) %>%
    arrange(txt)
#
ymax = 150
tp2 = tp0 %>%
    mutate(res = map(deg.reg, .f <- function(x) desc_stat(x$n))) %>%
    select(nid, res) %>%
    unnest() %>%
    spread(statK, statV) %>%
    mutate(q95 = ifelse(q95 > ymax, ymax, q95)) %>%
    inner_join(th, by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt))) %>%
    arrange(txt)

p1 = ggplot(tp1) +
    geom_point(aes(x=n.reg, y=n.tgt, size=deg.reg.q50, color=net_type)) +
    geom_text_repel(aes(x=n.reg, y=n.tgt, label=txt, color=net_type), size=2.5) +
    scale_x_continuous(name = '# TFs in network') +
    scale_y_continuous(name = '# Targets in network') +
    scale_color_aaas() +
    otheme(legend.pos = 'top.right', legend.dir = 'v',
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T)
p2 = ggplot(tp2) +
    geom_errorbar(aes(x=txt, y=mean, ymin=q5, ymax=q95, color=net_type), linetype='solid', size=.2, width=.4) +
    geom_crossbar(aes(x=txt, ymin=q25, y=q50, ymax=q75, color=net_type), fill='white', alpha=1, width=.8) +
    scale_x_discrete(expand = c(.02,0)) +
    scale_y_continuous(name = '# Targets per TF', limits=c(0,ymax), expand=expand_scale(mult=c(.01,.03))) +
    coord_flip() +
    scale_color_aaas() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
        xtick=T, ytick=T,
        xtitle=T, ytitle=F, xtext=T, ytext=T) +
    theme(axis.text.y = element_text(color=tp2$col))
fo = file.path(dirw, '03.stat.pdf')
ggpubr::ggarrange(p1, p2,
    nrow=1, ncol=2, widths = c(2.5,2), heights = c(1),
    labels=LETTERS[1:2]) %>%
    ggpubr::ggexport(filename = fo, width = 12, height = 8)
#}}}

#{{{ TF degree plot
tz = tp0 %>% select(nid, deg.reg) %>%
    unnest() %>% group_by(reg.gid) %>%
    summarise(n.net = sum(n>0), cv = sd(n)/mean(n)) %>%
    ungroup() %>%
    filter(n.net >= 20, cv > 1)
#}}}
#}}}

#{{{ Y1H eval
dirw = file.path(dird, '08_y1h')
fi = file.path(dirw, '01.rds')
y1h = readRDS(fi)
tn = y1h$t_y1h %>% transmute(reg.gid=reg, tgt.gid=tgt, Y1H='Yes')

tp = ev_tf %>% select(nid, ystat) %>% unnest() %>%
    filter(net_size == 50000) %>% select(-net_size) %>%
    unnest() %>%
    left_join(tn, by = c('reg.gid','tgt.gid'))
tp %>% count(reg.gid)
tp %>% count(tgt.gid)
tp %>% count(reg.gid,tgt.gid)
tp %>% filter(!is.na(Y1H)) %>% count(reg.gid,tgt.gid) %>% count(n)
tp %>% count(reg.gid,tgt.gid) %>% count(n)
tp %>% distinct(reg.gid,tgt.gid,Y1H) %>% count(Y1H)

to = tp %>%
    inner_join(th, by='nid') %>%
    select(reg.gid, tgt.gid, Y1H, txt) %>%
    mutate(txt = factor(txt, levels=th$txt)) %>%
    mutate(s = T) %>%
    spread(txt, s)
fo = file.path(dirw, '11.support.edges.tsv')
write_tsv(to, fo, na='')
#}}}

#{{{ known TF AUROC
dirw = file.path(dird, '14_eval_sum')
tp0 = ev_tf %>% select(nid, tfstat) %>% unnest()

#{{{ selected roc/prc plot
nids = c("n16b","n16c","n99b_1","nc03")
cols.dev = c(pal_npg()(4)[2:4], brewer.pal(6,"Paired")[6])
cols.dev = c(pal_npg()(8))
ss = th %>% distinct(nid,study) %>% filter(nid %in% nids) %>%
    arrange(nid) %>% pull(study)

tp1 = tp0 %>% inner_join(th, by='nid') %>%
    filter(nid %in% nids) %>%
    select(nid,study,note,ctag,roc) %>% unnest() %>%
    mutate(ctag = sprintf("%s AUROC", ctag)) %>%
    mutate(study=factor(study,levels=rev(ss)))
p1 = ggplot(tp1) +
    geom_line(mapping = aes(x = TPR, y = FPR, color = study)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'FPR: FP/(FP+TN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'TPR: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    otheme(strip.size = 9, legend.pos = 'none', margin = c(.5,.5,.1,.5),
           xtitle=T, ytitle=T, xtext=T, ytext=T, 
           xgrid=T, ygrid=T, xtick=T, ytick=T)
#
tp2 = tp0 %>% inner_join(th, by='nid') %>%
    filter(nid %in% nids) %>%
    select(nid,study,note,ctag,prc) %>% unnest() %>%
    mutate(ctag = sprintf("%s AUPR", ctag)) %>%
    mutate(study=factor(study,levels=rev(ss)))
p2 = ggplot(tp2) +
    geom_line(mapping = aes(x = recall, y = precision, color = study)) +
    #geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'Recall: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'Precision: TP/(TP+FP)', breaks=c(.25,.5,.75), limits=c(0,1), expand=c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    otheme(strip.size = 9, margin = c(.1,.5,.5,.5),
           xtitle=T, ytitle=T, xtext=T, ytext=T, 
           xgrid=T, ygrid=T, xtick=T, ytick=T) +
    theme(legend.position = c(.85,.25), legend.justification = c(0,0)) +
    guides(direction = 'vertical', color = guide_legend(ncol = 1, byrow = F)) 
#
fo = file.path(dirw, "05.roc_prc.dev.pdf")
ggpubr::ggarrange(p1, p2, nrow = 2, ncol = 1, heights = c(1,1))  %>%
    ggpubr::ggexport(filename = fo, width = 10, height = 5)
#}}}

#{{{ aupr/auroc bar-plot
fp = file.path(dirw, "05.auc.pdf")
wd = 8; ht = 12
tp = tp0 %>% select(nid,ctag,auroc,auprc) %>%
    #filter(ctag %in% ctags, nid %in% nids) %>%
    rename(AUPR = auprc, AUROC = auroc) %>%
    gather(type, auc, -nid,  -ctag) %>%
    mutate(auc = as.numeric(auc)) %>% filter(!is.na(auc)) %>%
    mutate(type = factor(type, levels = c("AUROC", "AUPR"))) %>%
    mutate(lab = str_remove(sprintf("%.03f", auc), '^0+')) %>%
    mutate(nid = factor(nid, levels = rev(th$nid)))
tpl = tp %>% distinct(type, ctag) %>% filter(type == 'AUROC') %>%
    mutate(itc.y = .5)
p1 = ggplot(tp, aes(x=nid, y=auc, fill=type)) +
    geom_bar(stat='identity', width=.75, alpha=.7) +
    geom_text(aes(label=lab), hjust=1, size=2) +
    geom_hline(data=tpl, aes(yintercept=itc.y), size=.3, alpha=.5) +
    scale_x_discrete(breaks=th$nid, labels=th$txt, expand=c(0,0)) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    scale_fill_npg() +
    coord_flip() +
    facet_wrap(type~ctag, scale = 'free_x', nrow = 2) +
    otheme(strip.size=7, legend.pos='none', margin=c(.2,.2,.2,.2),
           ygrid=T, xtitle=F, ytext=T) +
    theme(axis.text.y=element_text(color=th$col))
ggpubr::ggarrange(p1, nrow = 1, ncol = 1, labels = '', heights = c(2,2)) %>%
    ggpubr::ggexport(filename = fp, width = wd, height = ht)
#}}}
#}}}

#{{{ evaluate using GO/CornCyc
gs = read_gs()
fi = file.path(dirw, '01.go.rds')
ev_go = readRDS(fi)
net_sizes = c(1e4,5e4,1e5,5e5)
net_sizes = c(5e4,5e5)

#{{{ enrichment
tp1 = ev_go %>% select(nid, enrich) %>% unnest()
tp2 = ev_go %>% select(nid, enrich_term) %>% unnest() %>%
    group_by(nid, net_size, ctag) %>%
    summarise(n_grp=length(grp), n_grp_sig=sum(pval<.05)) %>%
    ungroup() %>%
    mutate(sigtxt = str_c(n_grp_sig,n_grp,sep='/')) %>%
    select(nid,net_size,ctag,sigtxt)
ctags = c("GO_HC", "GO_arabidopsis","GO_Interproscan5","CornCyc")
tp = tp1 %>% left_join(tp2, by=c('nid','net_size','ctag')) %>%
    mutate(sig = ifelse(pval<.05, 1, 0)) %>%
    filter(ctag %in% ctags) %>%
    mutate(ctag = factor(ctag, levels = ctags)) %>%
    inner_join(th, by='nid') %>%
    filter(net_size %in% net_sizes) %>%
    mutate(net_size=factor(net_size,levels=net_sizes)) %>%
    mutate(txt = factor(txt, levels=rev(th$txt)))
#
p1 = ggplot(tp, aes(txt, fc)) +
    geom_point(aes(color=net_size, shape=net_size, alpha=sig), size=2.5) +
    geom_text_repel(aes(label=sigtxt), size=2) +
    #geom_hline(yintercept = 1, alpha= .5, linetype='dotted') +
    scale_x_discrete(expand = expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name = 'Fold Enrichment', expand = c(.05,0)) +
    coord_flip() +
    facet_grid(.~ctag, scale='free') +
    scale_fill_npg() +
    scale_color_npg() +
    scale_alpha(range = c(.4,1), breaks=c(0,1), labels=c("pval > 0.05", "pval < 0.05")) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           strip.size = 8,
           xtitle=T, xtext=T, ytext=T, ygrid=T, xtick=T) +
    theme(axis.text.y = element_text(color = rev(th$col))) +
    theme(legend.box = "horizontal") +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.3)) +
    guides(color = guide_legend("network size:", nrow=1, order=1),
           shape = guide_legend("network size:", nrow=1, order=1),
           alpha = guide_legend("signifance:", nrow=1, order=2))
    #theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
fo = file.path(dirw, "14.go.pdf")
ggsave(p1, filename = fo, width = 10, height = 10)
#}}}

#{{{ GO heatmap
gotag = 'GO_uniprot.plants'
goname = gs$fun_ann %>% filter(ctag==gotag) %>% distinct(grp,note) %>%
    mutate(note = str_sub(note, 1, 65))
tx = ev_go %>% select(nid,enrich_term) %>%
    unnest() %>%
    filter(net_size==50000,ctag==gotag) %>%
    select(-net_size, -ctag) %>%
    filter(n>=50, pval<.05, fc>=2) %>%
    mutate(fc = log2(fc))

tx %>% inner_join(th, by = 'nid') %>%
    count(grp,net_type) %>% spread(net_type,nn) %>%
    select(-ril,-liftover) %>% inner_join(goname, by='grp') %>%
    print(n=40)

#{{{ #targets per TF histogram
tps = tp %>% count(grp)
p = ggplot(tps, aes(x=nn)) +
    geom_histogram() +
    scale_x_continuous(name = '# Enriched GRNs', expand=c(0,0)) +
    scale_y_continuous(name = '# GO Terms', expand=expand_scale(mult=c(0,.05))) +
    scale_color_aaas() +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '14.2.dist.pdf')
ggsave(p, file=fo, width=6, height=4)
#}}}

require(cluster)
require(ape)
grps = tx %>% count(grp) %>% filter(nn>=5) %>% pull(grp)
tp = tx %>% filter(grp %in% grps)

#{{{ GRN order
e = tp %>% select(nid,grp,fc) %>% spread(nid, fc) %>% select(-grp)
e[is.na(e)] = NA
dim(e)
#
cor_opt = "pearson"
hc_opt = "ward.D"
#edist <- as.dist(1-cor(e, method = cor_opt))
edist = daisy(t(e), metric = 'gower')
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
net_names = ehc$labels[ehc$order]
#}}}
#{{{ GO term order
e = tp %>% select(nid,grp,fc) %>% spread(grp, fc) %>% select(-nid)
e[is.na(e)] = NA
dim(e)
#
cor_opt = "pearson"
hc_opt = "ward.D"
#edist <- as.dist(1-cor(e, method = cor_opt))
edist = daisy(t(e), metric = 'gower')
edist[is.na(edist)] = 0
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
go_names = ehc$labels[ehc$order]
#}}}

tpg = tp %>% distinct(grp) %>% inner_join(goname,by='grp') %>%
    arrange(grp)
tpn = th %>% mutate(nid=factor(nid,levels=lnames)) %>% arrange(nid) %>%
    mutate(txt = factor(txt,levels=txt))
tpf = tp %>% mutate(nid = factor(nid,levels=lnames)) %>%
    mutate(grp = factor(grp, levels=go_names))
p = ggplot(tpf) +
    geom_tile(aes(x = nid, y = grp, fill = fc)) +
    #geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    #geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    #geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_discrete(breaks=tpn$nid, labels=tpn$txt, position='top', expand=c(0,0)) +
    scale_y_discrete(breaks=tpg$grp, labels=tpg$note, expand = c(0,0)) +
    #scale_fill_gradientn(name='log2FC', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    scale_fill_viridis(name='log2FC', direction=-1) +
    otheme(legend.pos='right', legend.dir='v',
           xtick=T,ytick=T,ytext = T) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x=element_text(color=tpn$col,angle=45,size=7,hjust=0,vjust=0)) +
    theme(legend.title = element_text())
fo = sprintf("%s/14.go.heat.pdf", dirw)
ggsave(p, file = fo, width = 10, height = 12)

tp = th %>% mutate(taxa = nid, lab = txt) %>%
    select(taxa, everything())
cols1 = c('gray80','black','red','seagreen3', pal_d3()(5))
p1 = ggtree(tree, layout = 'rectangular') +
    #geom_tiplab(size = labsize, color = 'black') +
    scale_x_continuous(expand = expand_scale(.9,.03)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp + geom_tiplab(aes(label = lab, color=net_type), family='mono') +
    scale_color_aaas()
fo = file.path(dirw, "gotree.pdf")
ggsave(p1, filename = fo, width=7, height=8)
#}}}
#}}}

t_link_sup = t_net %>% unnest() %>%
    right_join(t_link, by = c("reg.gid",'tgt.gid'))

t_link = gs$tf %>% filter(reg.gid == gid) %>% select(reg.gid,tgt.gid,binding)
#{{{ case study
gs = read_gs()
gid = gs$tfs %>% filter(ctag=='O2') %>% pull(reg.gid)
tgts = gs$tf %>% filter(reg.gid==gid) %>% pull(tgt.gid)
gcfg$gene.desc %>% filter(id %in% tgts) %>% print(n=40)
tq = tibble(gid = c(gid, tgts), tf = c(T, rep(F,length(tgts)))) %>%
    inner_join(gcfg$gene.desc, by = c('gid'='id')) %>%
    mutate(gname = ifelse(is.na(note1),
                          ifelse(is.na(note2), gid, note2), note1)) %>%
    select(gid, tf, gname)

study='mec03'
study='me16a'
th = rnaseq_sample_meta(study)
tm = rnaseq_sample_cpm(study)
ths = th %>%
    distinct(Tissue,Genotype,Treatment) %>%
    arrange(Tissue,Genotype,Treatment) %>%
    #mutate(txt = sprintf("%s--%s",Tissue,Treatment))
    mutate(txt = Genotype)

tp = tm %>% filter(gid %in% tq$gid) %>%
    select(gid,SampleID,CPM) %>%
    inner_join(th, by='SampleID') %>%
    group_by(gid,Genotype,Tissue,Treatment) %>%
    summarise(CPM = asinh(mean(CPM))) %>% ungroup() %>%
    #mutate(txt = sprintf("%s--%s",Tissue,Treatment)) %>%
    mutate(txt = Genotype) %>%
    inner_join(tq, by='gid') %>%
    mutate(gid = factor(gid, levels=tq$gid)) %>%
    mutate(txt = factor(txt, levels=ths$txt))

p = ggplot(tp) +
    geom_tile(aes(x = gid, y = txt, fill = CPM)) +
    #geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    #geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    #geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_discrete(expand = c(0,0), breaks=tq$gid, labels=tq$gname) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    #scale_fill_viridis() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xtext=T, ytext = T) +
    theme(axis.text.x=element_text(angle=30,size=8,hjust=1,vjust=1)) +
    theme(axis.text.y=element_text(size=7)) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 7))
fo = sprintf("%s/15.cpm.2.pdf", dirw)
ggsave(p, file = fo, width = 8, height = 15)
#
#}}}

#{{{ evaluate briggs data
#{{{ filter GRN using briggs DE info
fi = file.path(dirw, '01.br.rds')
ev_br = readRDS(fi)
ev_br_filt = ev_br %>%
    mutate(tgt.DE = ifelse(tgt.DE == 'non_DE', 'non_DE', 'DE')) %>%
    filter(!reg.DE %in% c('non_DE','DE1-2','DE2-4'), tgt.DE == 'DE') %>%
    mutate(drc = ifelse(reg.DEdir==tgt.DEdir, 1, -1)) %>%
    group_by(nid, reg.gid, tgt.gid) %>%
    summarise(n.tissue = n(), m.drc = sum(drc)/n.tissue)  %>%
    ungroup()
ev_br_filt %>% count(nid, n.tissue)

ev_br_sum = ev_br %>%
    group_by(nid, tissue, reg.DE) %>%
    summarise(nl = n(), n.reg = length(unique(reg.gid)),
              nl.tgt.de = sum(tgt.DE != 'non_DE'),
              pl.tgt.de = nl.tgt.de / nl) %>%
    ungroup()

fo = file.path(dirw, '01.br.filt.rds')
saveRDS(ev_br_filt, ev_br_sum, file = fo)
#}}}

br = read_briggs()
tissues = c('auricle_v12','ear_v14','embryo_27DAP','kernel_14DAP','root_0DAP',
            'seedlingmeristem_11DAS')
tissues = br$tissues
#
fi = file.path(dirw, '01.br.filt.rds')
ev_br_filt = readRDS(fi)
net_sizes = c(1e4, 5e4, 1e5)
net_size_map = c('1e4'=1e4, '5e4'=5e4, '1e5'=1e5)
net_size = 5e4
des = c("non_DE","DE1-2","DE2-4","DE4+","SPE")

#{{{ show support in Briggs dataset
th2 = th %>% transmute(nid=nid, lgd=txt)
tags = c("DE",'SPE')
tags = c("non_DE", "DE1-2",'DE2-4','DE4+','SPE')
tp = ev_br_filt$sum %>%
    mutate(txt = sprintf("%d", n.reg)) %>%
    mutate(lab = str_remove(sprintf("%.02f", pl.tgt.de), '^0+')) %>%
    rename(tag = reg.DE) %>%
    mutate(tag = factor(tag, levels = tags)) %>%
    #mutate(net_size = factor(net_size_map[as.character(net_size)], levels = net_sizes)) %>%
    inner_join(br$des, by = c("tissue"="Tissue")) %>%
    mutate(fc = pl.tgt.de/propDE) %>%
    mutate(tissue = factor(tissue, levels = tissues)) %>%
    inner_join(th2, by='nid') %>%
    mutate(nid = factor(nid, levels=th2$nid)) %>%
    mutate(lgd = factor(lgd, levels=th2$lgd))

#{{{ all
tp0 = tp %>%
    mutate(lgd = factor(lgd, levels=rev(th2$lgd)))
p1 = ggplot(tp0, aes(x=lgd)) +
    geom_point(aes(y=fc, color=tag, shape=tag), size=1) +
    scale_x_discrete(expand = c(.01,0)) +
    scale_y_continuous(name = 'Enrichment (fold change) in Target DE',
                       breaks = c(1,2,4),
                       expand = expand_scale(mult=c(0,.1))) +
    scale_fill_d3() +
    scale_shape_manual(values=pal_shapes()) +
    scale_color_d3() +
    coord_flip() +
    facet_wrap(~tissue, nrow=2) +
    otheme(legend.pos = 'top.center.out',
           margin = c(.2,.2,.2,1),
           strip.size = 9,
           xtitle=T, xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(axis.text.y=element_text(color=rev(th$col))) +
    theme(strip.text.y = element_text(angle=0,hjust=.5)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', color=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12.br.1.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}
#{{{ all 2
tp0 = tp %>%
    mutate(tissue = factor(tissue, levels=rev(tissues)))
p1 = ggplot(tp0, aes(x=tissue)) +
    geom_point(aes(y=fc, color=tag, shape=tag), size=1.5) +
    #scale_x_discrete(breaks=th$nid, labels=th$txt, expand = c(.01,0)) +
    scale_y_continuous(name = 'Enrichment (fold change) in Target DE',
                       breaks = c(1,2,4),
                       expand = expand_scale(mult=c(.1,.1))) +
    scale_fill_d3() +
    scale_shape_manual(values=pal_shapes()) +
    scale_color_d3() +
    coord_flip() +
    facet_wrap(~lgd, ncol=10) +
    otheme(legend.pos = 'top.center.out',
           strip.size = 7,
           xtitle=T, xtext=T, xtick=T, ytext=T, xgrid=T, ygrid=T) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', color=guide_legend(nrow=1),
        shape=guide_legend(nrow=1))
fp = sprintf("%s/12.br.2.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 14)
#}}}

#{{{ tissue network
nid = 'n99a_1' # kaeppler endosperm
nid = 'n15a' # leiboff SAM
nid = 'n14a' # hirsch seedling
nid = 'n16a' # jin kernel
txt = th %>% filter(nid==!!nid) %>% pull(txt)
tp0 = tp %>% filter(nid == !!nid)
p1 = ggplot(tp0, aes(x=tissue)) +
    geom_point(aes(y=pl.tgt.de, color=tag, shape=tag), size=2) +
    #geom_text(aes(y=pl.tgt.de+.01, group=tag, label=lab), hjust=0, position=position_dodge(width=1), size=2) +
    #geom_text(aes(y=.01, group=tag, label=txt), hjust=0, position=position_dodge(width=1), size=2) +
    #geom_hline(data=tpl, aes(yintercept=prop.de), size=.3, alpha=.5) +
    #scale_x_discrete(breaks=th$nid, labels=th$txt, expand = c(.01,0)) +
    scale_y_continuous(name = 'Proportion DE Targets', expand = expand_scale(mult=c(0,.1))) +
    scale_shape_manual(name = txt, values=pal_shapes()) +
    scale_color_d3(name = txt) +
    otheme(legend.pos = 'top.center.out', legend.title = T,
           margin = c(.5,.2,.2,1),
           ytitle=T, xtext=T, ytick=T, ytext=T, xgrid=T, ygrid=F) +
    theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1)) +
    #theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction='horizontal',
           color=guide_legend(title.position='left',nrow=1),
           shape=guide_legend(title.position='left',nrow=1))
fp = sprintf("%s/13.br.%s.pdf", dirw, nid)
ggsave(p1, filename = fp, width = 9, height = 4)

#}}}
#}}}


#{{{ how often is a link observed in multiple tissues? is it consistent?
tp = ev_br_filt %>% filter(n.tissue >= 5) %>%
    inner_join(th[,c('nid','txt','col')], by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt))) #%>%
    #mutate(reg.DE = factor(reg.DE, levels = rev(des)))
tps = tp %>% count(nid) %>% inner_join(th, by='nid') %>%
    mutate(txt = factor(txt, levels=th$txt))
p = ggplot(tp, aes(x=txt,y=m.drc)) +
    geom_violin() +
    geom_text(data=tps, aes(x=txt, y=0, label=n), hjust=.5, size=2.5) +
    scale_y_continuous(limits=c(-1,1), expand=expand_scale(mult=c(.01,.01))) +
    scale_color_aaas() +
    coord_flip() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '12.br.2.dir.tissue.pdf')
ggsave(p, file=fo, width=6, height=8)
#}}}

#{{{ if a TF has multiple targets, does it show consistent +/- effect?
tp = ev_br_filt %>% group_by(nid, reg.gid) %>%
    summarise(n.tissues=sum(n.tissue), ntgt = n(),
              m.drc = sum(n.tissue*m.drc)/sum(n.tissue)) %>%
    ungroup() %>%
    filter(n.tissues >= 10) %>%
    inner_join(th[,c('nid','txt','col')], by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt)))
tps = tp %>% count(nid) %>% inner_join(th, by='nid') %>%
    mutate(txt = factor(txt, levels=th$txt))
p = ggplot(tp, aes(x=txt,y=m.drc)) +
    geom_violin() +
    geom_text(data=tps, aes(x=txt, y=0, label=n), hjust=.5, size=2.5) +
    scale_y_continuous(limits=c(-1,1), expand=expand_scale(mult=c(.01,.01))) +
    scale_color_aaas() +
    coord_flip() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '12.br.2b.dir.tissue.pdf')
ggsave(p, file=fo, width=6, height=8)

gid0 = 'Zm00001d004230'
br$de %>% filter(gid==gid0) %>% print(n=23)

ev_br %>% filter(reg.gid==gid0, tgt.DE != 'non_DE', nid == 'n17a') %>%
    group_by(nid,tissue,tgt.DEdir) %>%
    summarise(nt = n())

tx %>% filter(reg.gid==gid0) %>%
    mutate(drc = ifelse(nt1<nt2, 'B<M', 'B>M')) %>%
    group_by(nid) %>%
    summarise(p.tgt = sum(drc=='B<M')/n()) %>% ungroup()
#}}}

#{{{ if a link is observed in multiple networks, is it consistent?
tp = ev_br_filt %>% group_by(reg.gid,tgt.gid) %>%
    summarise(nnet = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    filter(nnet >= 3)
p = ggplot(tp, aes(x=0, y=m.drc)) +
    geom_violin() +
    #geom_text(data=tps, aes(x=reg.DE, y=0, label=n), hjust=.5, size=3) +
    geom_text(x=0, y=0, label=nrow(tp), hjust=.5, size=3) +
    scale_y_continuous(limits=c(-1,1), expand=expand_scale(mult=c(.01,.01))) +
    scale_color_aaas() +
    coord_flip() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '12.br.3.dir.net.pdf')
ggsave(p, file=fo, width=6, height=5)
#}}}

#{{{ if a TF is observed in multiple networks, is it consistent?
tg = ev_br_filt %>% group_by(reg.gid, tgt.gid) %>%
    summarise(nnet = n(), m.drc = mean(m.drc),
              nn1 = sum(nid %in% nids_geno),
              nn2 = sum(nid %in% nids_dev)) %>%
    ungroup() %>%
    filter(nn1 >= 1, nn2 >= 1, nnet >= 2) %>%
    filter(abs(m.drc) > .5)
tg %>% count(nnet)

tp = tg %>% group_by(reg.gid) %>%
    summarise(n.tgt = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    filter(n.tgt >= 3)
tp %>% arrange(desc(n.tgt)) %>% print(n=40)
p = ggplot(tp, aes(x=m.drc)) +
    geom_histogram() +
    scale_x_continuous(name = 'effect on target (positive or negative)', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name = '# TFs', expand=expand_scale(mult=c(0,.05))) +
    scale_color_aaas() +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '12.br.2b.dir.tissue.pdf')
ggsave(p, file=fo, width=6, height=4)
#}}}

#{{{ obtain high quality edges
tg = ev_br_filt %>% group_by(reg.gid, tgt.gid) %>%
    summarise(nnet = n(), m.drc = mean(m.drc)) %>% ungroup() %>%
    filter(nnet >= 3) %>%
    filter(abs(m.drc) > .5)
tg %>% count(nnet)
#}}}

#}}}

#{{{ eval using biomap data - correlation
fi = file.path(dirw, '01.bm.rds')
ev_bm = readRDS(fi)

ev_bm_nr = ev_bm %>% distinct(reg.gid,tgt.gid,Tissue,pcc)

tp = tg %>% inner_join(ev_bm_nr, by=c('reg.gid','tgt.gid')) 
p = ggplot(tp, aes(x = pcc)) +
    geom_histogram() +
    facet_wrap(~Tissue, ncol=1)+
    scale_fill_futurama() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v') +
    theme(strip.text.y = element_text(angle=0))
fo = file.path(dirw, '15.biomap.1.pcc.pdf')
ggsave(p, file=fo, width=6, height=8)

tp = ev_br_filt %>% inner_join(ev_bm, by = c('nid','reg.gid','tgt.gid')) %>%
    inner_join(th[,c('nid','txt','col')], by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt)))
tps = tp %>% count(nid,Tissue) %>%
    inner_join(th[,c('nid','txt','col')], by = 'nid') %>%
    mutate(txt = factor(txt, levels=th$txt))

p = ggplot(tp, aes(x = txt, y = pcc)) +
    geom_violin() +
    geom_text(data=tps, aes(x=txt, y=0, label=n), hjust=.5, size=3) +
    coord_flip() +
    facet_wrap(~Tissue, ncol=5) +
    scale_fill_futurama() +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v') +
    theme(strip.text.y = element_text(angle=0))
fo = file.path(dirw, '15.biomap.1.pcc.pdf')
ggsave(p, file=fo, width=10, height=8)
#}}}


