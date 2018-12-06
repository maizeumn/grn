source("functions.R")
dirw = file.path(dird, '14_eval_sum')
diri = '~/projects/rnaseq'

#{{{ # manually merge tissue GRNs
if(F) {
study = "lin2017"
study = "kremling2018"
ths = th %>% filter(study == !!study, nid != 'n50', nid != 'n60')

tno = tibble()
for (i in 1:nrow(ths)) {
    nid = ths$nid[i]
    fgn = ths$fgn[i]
    x = load(fgn)
    x
    if( i == 1 ) {
        tno = tn %>%
            mutate(weight.a = weight, weight.m = weight) %>%
            select(-weight)
    } else {
        tno = tno %>% full_join(tn, by = c('reg.gid', 'tgt.gid')) %>%
            replace_na(list(weight = 0)) %>%
            mutate(weight.a = weight.a + weight,
                   weight.m = weight.m * weight) %>%
            select(-weight)
    }
    cat(nid,"\n")
}

fo1 = file.path(dirw, "12_output/n58.rda")
fo2 = file.path(dirw, "12_output/n59.rda")
fo1 = file.path(dirw, "12_output/n68.rda")
fo2 = file.path(dirw, "12_output/n69.rda")
rids = tno %>% distinct(reg.gid) %>% pull(reg.gid)
tids = tno %>% distinct(tgt.gid) %>% pull(tgt.gid)
tn1 = tno %>% rename(weight = weight.a)
tn2 = tno %>% rename(weight = weight.m)
save(rids, tids, tn1, file = fo1)
save(rids, tids, tn2, file = fo2)
}
#}}}

fi = file.path(dirw, '01.tf.rds')
ev_tf = readRDS(fi)

#{{{ general network stats
tp = ev_tf %>% select(nid, nstat) %>% unnest() %>%
    filter(net_size==50000) %>%
    inner_join(th, by = 'nid') %>%
    mutate(txt = factor(txt, levels=rev(th$txt))) %>%
    arrange(txt)
ymax = 150
tp$q95[tp$q95 > ymax] = ymax

p1 = ggplot(tp) +
    geom_point(aes(x=n.reg, y=n.tgt, size=q50, color=net_type)) +
    geom_text_repel(aes(x=n.reg, y=n.tgt, label=txt, color=net_type), size=2.5) +
    scale_x_continuous(name = '# TFs in network') +
    scale_y_continuous(name = '# Targets in network') +
    scale_color_aaas() +
    otheme(legend.pos = 'top.right', legend.dir = 'v',
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T)
p2 = ggplot(tp) +
    geom_errorbar(aes(x=txt, y=mean, ymin=q5, ymax=q95, color=net_type), linetype='solid', size=.2, width=.4) +
    geom_crossbar(aes(x=txt, ymin=q25, y=q50, ymax=q75, color=net_type), fill='white', alpha=1, width=.8) +
    scale_x_discrete(expand = c(.02,0)) +
    scale_y_continuous(name = '# Targets per TF', limits=c(0,ymax), expand=expand_scale(mult=c(.01,.03))) +
    coord_flip() +
    scale_color_aaas() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
        xtick=T, ytick=T,
        xtitle=T, ytitle=F, xtext=T, ytext=T) +
    theme(axis.text.y = element_text(color=tp$col))
fo = file.path(dirw, '03.stat.pdf')
ggpubr::ggarrange(p1, p2,
    nrow=1, ncol=2, widths = c(2.5,2), heights = c(1),
    labels=LETTERS[1:2]) %>%
    ggpubr::ggexport(filename = fo, width = 12, height = 8)
#}}}

#{{{ Y1H eval
dirw = file.path(dird, '08_y1h')
fi = file.path(dirw, '01.rds')
y1h = readRDS(fi)
tn = y1h$tn %>% transmute(reg.gid=reg, tgt.gid=tgt, Y1H='Yes')

tp = ev_tf %>% select(nid, ystat) %>% unnest() %>%
    filter(net_size == 50000) %>% select(-net_size) %>%
    left_join(tn, by = c('reg.gid','tgt.gid'))
tp %>% count(reg.gid)
tp %>% count(tgt.gid)
tp %>% count(reg.gid,tgt.gid)
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
fi = file.path(dirw, '01.tf.rds')
ev_tf = readRDS(fi)

#{{{ selected roc/pr plot
nids = c("n16b","n16c","n99b_1","nc03")
cols.dev = c(pal_npg()(4)[2:4], brewer.pal(6,"Paired")[6])
cols.dev = c(pal_npg()(8))
ss = th %>% distinct(nid,study) %>% filter(nid %in% nids) %>%
    arrange(nid) %>% pull(study)

tp1 = ev_tf %>% inner_join(th, by='nid') %>%
    filter(nid %in% nids) %>%
    select(nid,study,note,roc) %>% unnest() %>%
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
tp2 = ev_tf %>% inner_join(th, by='nid') %>%
    filter(nid %in% nids) %>%
    select(nid,study,note,pr) %>% unnest() %>%
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
fo = file.path(dirw, "10.dev.atlas.pdf")
ggarrange(p1, p2, nrow = 2, ncol = 1, heights = c(1,1))  %>%
    ggexport(filename = fo, width = 10, height = 5)
#}}}

#{{{ aupr/auroc bar-plot
fp = file.path(dirw, "05.auc.pdf")
wd = 8; ht = 12
tp = ev_tf %>% select(nid,auroc,aupr) %>% unnest() %>% select(-ctag1) %>%
    #filter(ctag %in% ctags, nid %in% nids) %>%
    rename(AUPR = aupr, AUROC = auroc) %>%
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
ggarrange(p1, nrow = 1, ncol = 1, labels = '', heights = c(2,2)) %>%
    ggexport(filename = fp, width = wd, height = ht)
#}}}
#}}}

#{{{ evaluate using GO/CornCyc
fi = file.path(dirw, '01.go.rds')
ev_go = readRDS(fi)
net_sizes = c(1e4,5e4,1e5,5e5)
net_sizes = c(5e4,1e5,5e5)

#{{{ enrichment
tps = ev_go %>% distinct(nid,txt,col)
tp = ev_go %>% select(-study,-note,-sample_size,-fi) %>% unnest() %>%
    filter(n.tgt >= 10, net_size %in% net_sizes) %>%
    group_by(nid,net_size,grp_tag) %>%
    summarise(fc = sum(pairs.coreg.obs)/sum(pairs.coreg.exp)) %>%
    ungroup() %>%
    mutate(net_size=factor(net_size,levels=net_sizes),
           nid = factor(nid, levels = rev(tps$nid)))
#
p1 = ggplot(tp) +
    geom_point(aes(nid,fc,color=net_size,shape=net_size), size=2.5) +
    geom_hline(yintercept = 1, alpha= .5, linetype='dotted') +
    scale_x_discrete(breaks = tps$nid, labels = tps$txt, expand = expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name = 'Fold Enrichment', expand = c(.05,0)) +
    coord_flip() +
    facet_grid(.~grp_tag) +
    scale_fill_npg() +
    scale_color_npg() +
    otheme(legend.pos='bottom.right',
           strip.size = 8, margin = c(.5,.5,.5,.5),
           xtitle=T, xtext=T, ytext=T, ygrid=T, xtick=T) +
    theme(axis.text.y = element_text(color = tps$col)) +
    guides(color = guide_legend(ncol = 1, byrow = F))
    #theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
fo = file.path(dirw, "14.go.pdf")
ggsave(p1, filename = fo, width = 8, height = 10)
#}}}

#{{{
tp = ev_go %>% select(net_type,nid,col,txt,res) %>%
    unnest() %>%
    filter(net_size==50000,grp_tag=='GO',n.tgt>=50,fc>=3, fc!=Inf) %>%
    select(net_type,nid,col,txt,grp,pairs.coreg.obs,pairs.coreg.exp,fc)

goname = gs$grp %>% distinct(grp,note)
tp %>% count(grp,net_type) %>% spread(net_type,n) %>%
    select(-ril,-liftover) %>% inner_join(goname, by='grp') %>%
    pull(note)

tps = tp %>% count(grp)
p = ggplot(tps, aes(x=n)) +
    geom_histogram() +
    scale_x_continuous(name = '# Enriched GRNs', expand=c(0,0)) +
    scale_y_continuous(name = '# GO Terms', expand=expand_scale(mult=c(0,.05))) +
    scale_color_aaas() +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
           legend.pos = 'top.right', legend.dir = 'v')
fo = file.path(dirw, '14.2.dist.pdf')
ggsave(p, file=fo, width=6, height=4)


to = tp %>% mutate(txt = factor(txt,levels=th$txt)) %>%
    select(txt,grp,fc) %>%
    spread(txt,fc) %>%
    inner_join(goname, by = 'grp') %>%
    select(grp,note, everything())
fo = file.path(dirw, '41.enriched.go.tsv')
write_tsv(to, fo, na = '')
#}}}

#{{{
require(cluster)
goname = gs$grp %>% distinct(grp,note)
gos = ev_go %>% select(net_type,nid,col,txt,res) %>% unnest() %>%
    filter(net_size==50000,grp_tag=='GO',n.tgt>=50,fc>=3, fc!=Inf) %>%
    count(grp) %>% filter(n>=5, n<10) %>%
    pull(grp)
tp = ev_go %>% select(nid,res) %>% unnest() %>%
    filter(net_size==50000,grp_tag=='GO',grp %in% gos) %>%
    select(nid, grp, fc) %>%
    mutate(fc = log2(fc)) %>%
    mutate(fc = ifelse(is.finite(fc), fc, NA))
e = tp %>% spread(nid, fc) %>% select(-grp)
e[is.na(e)] = NA
dim(e)

cor_opt = "pearson"
hc_opt = "ward.D"
hc_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
#edist <- as.dist(1-cor(e, method = cor_opt))
edist = daisy(t(e), metric = 'gower')
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]

tpg = tp %>% distinct(grp) %>% inner_join(goname,by='grp')
tpn = th %>% mutate(nid=factor(nid,levels=lnames)) %>% arrange(nid) %>%
    mutate(txt = factor(txt,levels=txt))
tp = tp %>% mutate(nid = factor(nid,levels=lnames))
p = ggplot(tp) +
    geom_tile(aes(x = nid, y = grp, fill = fc)) +
    #geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    #geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    #geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_discrete(breaks=tpn$nid, labels=tpn$txt, position='top', expand=c(0,0)) +
    scale_y_discrete(breaks=tpg$grp, labels=tpg$note, expand = c(0,0)) +
    scale_fill_gradientn(name='log2FC', colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    #scale_fill_viridis() +
    otheme(legend.pos='right', legend.dir='v',
           xtick=T,ytick=T,ytext = T) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x=element_text(color=tpn$col,angle=45,size=7,hjust=0,vjust=0)) +
    theme(legend.title = element_text())
fo = sprintf("%s/14.go.heat.pdf", dirw)
ggsave(p, file = fo, width = 12, height = 12)


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
fo = file.path(dirw, '01.br.filt.rds')
saveRDS(ev_br_filt, file = fo)
#}}}

fi = file.path(dirw, '01.br.filt.rds')
ev_br_filt = readRDS(fi)
tissues6 = c('auricle_v12','ear_v14','embryo_27DAP','kernel_14DAP','root_0DAP',
            'seedlingmeristem_11DAS')
net_sizes = c(1e4, 5e4, 1e5)
net_size_map = c('1e4'=1e4, '1e5'=1e5, '1e6'=1e6)
net_size = 1e4
des = c("non_DE","DE1-2","DE2-4","DE4+","SPE")

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

#{{{
tp0 = ev_br %>% #filter(net_size == !!net_size) %>%
    filter(tissue %in% tissues6) %>%
    group_by(nid, tissue, net_size) %>%
    summarise(n.links = sum(n),
              n.reg.de = sum(n[reg.DE]),
              p.reg.de = sum(n[reg.DE & tgt.DE])/n.reg.de
    ungroup()

tp1 = tp0 %>% select(nid,tissue,net_size,n.links,starts_with("n.reg.")) %>%
    gather(tag, n.reg, -nid,-tissue,-net_size,-n.links) %>%
    mutate(tag = str_replace(tag, "^n\\.reg\\.", ""))
tp2 = tp0 %>% select(nid,tissue,net_size,n.links,starts_with("p.reg.")) %>%
    gather(tag, p.reg, -nid,-tissue,-net_size,-n.links) %>%
    mutate(tag = str_replace(tag, "^p\\.reg\\.", ""))
tags = c("DE",'DE2','DE4','DE8','SPE')
tags = c("DE",'SPE')
tp = tp1 %>%
    inner_join(tp2, by = c('nid','tissue','net_size','n.links','tag')) %>%
    filter(net_size == !!net_size) %>%
    mutate(txt = sprintf("%d", n.reg)) %>%
    mutate(lab = str_remove(sprintf("%.02f", p.reg), '^0+')) %>%
    mutate(tag = toupper(tag)) %>%
    filter(tag %in% tags) %>%
    mutate(tag = factor(tag, levels = tags)) %>%
    #mutate(net_size = factor(net_size_map[as.character(net_size)], levels = net_sizes)) %>%
    mutate(tissue = factor(tissue, levels = tissues6)) %>%
    filter(nid %in% th$nid) %>%
    mutate(nid = factor(nid, levels=rev(th$nid)))
tpl = br$des %>% transmute(tissue = Tissue, prop.de = propDE) %>%
    filter(tissue %in% tissues6)
    #complete(study, tissue, net_size, reg.DE, nesting(tgt.DE),
    #         fill = list(n=0, nreg=0, ntgt=0))

p1 = ggplot(tp, aes(x=nid)) +
    geom_bar(aes(y=p.reg, fill=tag), stat='identity', position='dodge', width=.9, alpha=.7) +
    geom_text(aes(y=p.reg+.01, group=tag, label=lab), hjust=0, position=position_dodge(width=1), size=2) +
    geom_text(aes(y=.01, group=tag, label=txt), hjust=0, position=position_dodge(width=1), size=2) +
    geom_hline(data=tpl, aes(yintercept=prop.de), size=.3, alpha=.5) +
    scale_x_discrete(breaks=th$nid, labels=th$txt, expand = c(.01,0)) +
    scale_y_continuous(name = 'Prop. DE Targets', expand = expand_scale(mult=c(0,.1))) +
    scale_fill_d3() +
    coord_flip() +
    facet_wrap(~tissue, nrow = 1, scale = 'free_x') +
    otheme(strip.size = 9, margin = c(1.5,.2,.2,.2),
           xtitle = T, ytext = T, xgrid = F, ygrid = T) +
    theme(axis.text.y=element_text(color=th$col)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F))
fp = sprintf("%s/12.briggs.%d.pdf", dirw, net_size)
ggsave(p1, filename = fp, width = 8, height = 8)
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



#{{{ # use PCC to predict TF targets
if(F) {
dirw = file.path(dird, '13_auroc')
#{{{ old codes
    pcc.matrix = cor(asinh(expr), method = 'pearson')
    pcc = pcc.matrix[lower.tri(pcc.matrix)]

    ii = which(lower.tri(pcc.matrix))
    colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
    rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

    pcc[pcc == 1] = 0.999999
    pcc[pcc == -1] = -0.999999
    #pcc2 = log((1+pcc) / (1-pcc)) / 2
    pcc2 = atanh(pcc)
    coexv = (pcc2 - mean(pcc2)) / sd(pcc2)

    ord = order(-coexv)
    idxs = ord[1:(1*1000000)]
    dw = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], coex = coexv[idxs])
#}}}

#{{{ run
th2 = th %>% mutate(fi = sprintf("%s/11_input/%s.rda", dird, nid)) %>%
    filter(file.exists(fi))
#th2 = th2 %>% filter(nid %in% c("n41", "n42", "n71", "n81", "n82")) 

t_auc = tibble()
for (i in 1:nrow(th2)) {
    nid = th2$nid[i]; study = th2$study[i]; tag = th2$tag[i]; fi = th2$fi[i]
    x = load(fi)
    x
    #
    t2_tf = t_exp %>% filter(gid %in% t_gss$reg.gid) %>%
        group_by(gid) %>%
        summarise(ncond = sum(CPM>=1), pcond = ncond/n(), expressed = pcond >= .1,
                  CPM.mean = mean(CPM), CPM.var = var(CPM),
                  CPM.median = median(CPM), CPM.mad = mad(CPM)) %>%
        ungroup() %>% inner_join(t_gss, by = c('gid'='reg.gid')) %>%
        select(-ncond,-pcond)
    #
    n_cond = length(unique(t_exp$condition))
    gids = t_exp %>% group_by(gid) %>%
        summarise(n_cond_exp = sum(CPM >= 1)) %>%
        filter(n_cond_exp >= n_cond * .1) %>%
        pull(gid)
    t_flt = t_exp %>% filter(gid %in% gids)
    #
    use_cpm = T
    if(use_cpm) {
        t_flt = t_flt %>% mutate(exp.val = asinh(CPM))
    } else {
        t_flt = t_flt %>% mutate(exp.val = asinh(FPKM))
    }
    #
    et_b = t_flt %>% 
        select(condition, gid, exp.val) %>%
        spread(condition, exp.val)
    em_b = as.matrix(et_b[,-1])
    rownames(em_b) = et_b$gid
    #
    gids = rownames(em_b)
    pcc.matrix = cor(t(em_b), method = 'pearson')
    #
    rids = t_gss$reg.gid
    rids = rids[rids %in% gids]
    t_pcc = pcc.matrix[,rids] %>% as_tibble() 
    if(length(rids) == 1) colnames(t_pcc)[1] = rids[1]
    t_pcc = t_pcc %>% mutate(tgt.gid = gids) %>%
        gather(reg.gid, pcc, -tgt.gid) %>%
        filter(reg.gid != tgt.gid) %>%
        inner_join(t_gss, by = 'reg.gid') %>%
        left_join(t_gs, by = c('ctag','reg.gid','tgt.gid')) %>%
        replace_na(list(binding = 0)) %>%
        mutate(categ = factor(binding))
    #
    ctagss = t_gss %>% filter(reg.gid %in% rids) %>% pull(ctag)
    for (ctag in ctagss) {
        ty = t_pcc %>% filter(ctag == !!ctag)
        categ = as.integer(as.character(ty$categ))
        score = as.numeric(ty$pcc)
        score = as.numeric(abs(ty$pcc))
        if(max(score) == 0) score[1] = 0.1
        pr = pr.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        roc = roc.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        aupr = pr$auc.integral
        auroc = roc$auc
        t_auc1 = tibble(nid = nid, study = study, tag = tag, 
                        ctag = ctag, aupr = aupr, auroc = auroc)
        t_auc = rbind(t_auc, t_auc1)
        cat(sprintf("%s %10s aupr [%.04f] auroc [%.04f]\n", nid, ctag, aupr, auroc))
    } 
}
fo = file.path(dirw, '21.pcc.rda')
save(t_auc, file = fo)
#}}}

#{{{ plot
fi = file.path(dirw, '21.pcc.rda')
x = load(fi)
tpa = t_auc %>% 
    mutate(tag = sprintf("%s %s", study, tag)) %>%
    select(-study) 

for (opt in names(lst_eval)) {
    nids = lst_eval[[opt]]
    tp = tpa %>% filter(nid %in% nids) %>% 
        mutate(nid = factor(nid, levels = rev(nids)))
    fp = sprintf("%s/06.pcc.%s.pdf", dirw, opt)
    wd = 7; ht = 5
    if (opt == 'all') ht = 7
    auc_barplot(tp, fp, wd, ht)
}
#}}}
}
#}}}

