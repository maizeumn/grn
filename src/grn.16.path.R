source("functions.R")
require(ggpmisc)
require(ggraph)
require(tidygraph)
require(gridExtra)
diri = '~/projects/rnaseq'
dirw = file.path(dird, '16_pathways')
th = t_cfg %>% mutate(lgd=str_replace(lgd, " \\[.*\\]", '')) %>% select(nid,lgd,col)
tf = read_tf_info()
fp = file.path(dirw, 'pathways.xlsx')
#
fv = sprintf("%s/rf.100k.rds", dirr)
tx = readRDS(fv) %>% select(nid, tn) %>% unnest(tn)

#{{{ plot functions
str_split_s <- function(s, sep='') str_split(s, sep)[[1]]
make_cc_step <- function(reacs, prods, gids, spon, ta=symb, extra_filt=c()) {
    #{{{
    known_filt = c('ATP','CTP','AMP','CMP','ADP','UDP',
                   'NADPH','NADP+','NAD(P)+','NAD(P)H',
                   'H2O','H+','oxygen','CO2',
                   'Mg2+','phosphate','diphosphate',
                   'pyruvate','fumarate','succinate')
    filt = unique(c(extra_filt, known_filt))
    reacs = reacs[!reacs %in% filt]
    prods = prods[!prods %in% filt]
    sym_str = ''
    if(spon | (length(gids) == 1 & gids[1] == ''))
        sym_str = ''
    else {
        symbs = tibble(gid = gids) %>% distinct(gid) %>%
            left_join(ta, by='gid') %>%
            mutate(symbol = ifelse(is.na(symbol), trim_gid(gid), symbol)) %>%
            distinct(symbol) %>%
            mutate(s2 = ifelse(str_detect(symbol,'^0'), str_c('z',symbol,sep=''), symbol)) %>%
            arrange(s2) %>% pull(symbol)
        sym_str = str_c(symbs, sep=' ')
    }
    sym_str = str_to_upper(sym_str)
    if(length(reacs) == 0) reactants = 'None'
    if(length(prods) == 0) products = 'None'
    crossing(snode=reacs,enode=prods,lab=sym_str) %>% as_tibble() %>%
        distinct(snode,enode,lab) %>%
        mutate(lab2 = ifelse(str_detect(lab,'^0'), str_c('z',lab,sep=''), lab)) %>%
        arrange(snode,enode,lab2) %>%
        group_by(snode,enode) %>% summarise(lab=str_c(lab, collapse=' ')) %>%
        ungroup()
    #}}}
}
plotNet_cc <- function(path, cc, ta=symb, extra_filt=c(), lay='kk', angle='along') {
    #{{{
    ccm = cc %>% filter(pathway==!!path) %>%
        mutate(net = pmap(list(reactants,products,gids,spontaneous), make_cc_step, symb, extra_filt=!!extra_filt)) %>%
        select(net) %>% unnest(net)
    # ad hoc adjustments
    ccm = ccm %>% filter(snode != 'isopentenyl diphosphate' | enode != 'dimethylallyl diphosphate')
    #
    e0 = ccm %>% select(from=snode, to=enode, lab)
    g = as_tbl_graph(e0)
    circular = ifelse(lay == 'linear', T, F)
    pg = ggraph(g, layout=lay, circular=!!circular) +
        geom_edge_link(aes(label = lab), angle_calc=angle,
                       #label_dodge = unit(.2, 'cm'),
                       arrow=arrow(length=unit(.1,'cm'), angle=20, type='closed'),
                       lwd=.1, label_size=2.5, label_colour='darkblue',
                       start_cap=circle(.2,'cm'), end_cap=circle(.15,'cm')) +
        geom_node_text(aes(label=name), vjust=.5, size=2.5, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        otheme(panel.border=F, margin=c(0,0,0,0)) +
        guides(color=F, shape=T)
    pg
    #}}}
}
plotTable <- function(tbl, nr=1, size=8) {
    #{{{
    #est = tibble(x=.5,y=.5, tb=list(es))
    #pt = ggplot() +
        #geom_table_npc(data=est, aes(npcx=x,npcy=y,label=tb), size=2.5, color='black', vjust=.5, hjust=.5) +
        #otheme(margin=c(0,0,0,0), panel.border=F)
    tt1 = ttheme_default(base_size=size,
        core=list(fg_params=list()),
        colhead=list(fg_params=list(col="navyblue", fontface=4L))
    )
    r = nrow(tbl)
    if(nr == 1) {
        grob1 = tableGrob(tbl, rows=NULL, theme=tt1)
        ggarrange(grob1)
    } else if(nr == 2) {
        pr = ceiling(r/2)
        grob1 = tableGrob(tbl[1:pr,], rows=NULL, theme=tt1)
        grob2 = tableGrob(tbl[(pr+1):r,], rows=NULL, theme=tt1)
        ggarrange(grob1, grob2, nrow=1, ncol=2)
    } else if(nr == 3) {
        pr = ceiling(r/3)
        grob1 = tableGrob(tbl[1:pr,], rows=NULL, theme=tt1)
        grob2 = tableGrob(tbl[(pr+1):2*pr,], rows=NULL, theme=tt1)
        grob3 = tableGrob(tbl[(2*pr+1):r,], rows=NULL, theme=tt1)
        ggarrange(grob1, grob2, grob3, nrow=1, ncol=3)
    }
    #}}}
}
plotNets1 <- function(ti, th, ta=symb, lays='kk', angle='along') {
    #{{{
    x = ti %>%
        left_join(ta, by=c('reg.gid'='gid')) %>% rename(reg=symbol) %>%
        mutate(reg = ifelse(is.na(reg), trim_gid(reg.gid), reg)) %>%
        left_join(ta, by=c('tgt.gid'='gid')) %>% rename(tgt=symbol) %>%
        mutate(tgt = ifelse(is.na(tgt), trim_gid(tgt.gid), tgt)) %>%
        mutate(reg = str_to_sentence(reg), tgt = str_to_sentence(tgt)) %>%
        left_join(th, by='nid') %>%
        mutate(lgd = factor(lgd, levels=th$lgd))
    e0 = x %>% select(path, from=reg, to=tgt, lgd, score)
    chs = c(letters, LETTERS)
    es = e0 %>% distinct(lgd) %>% filter(!is.na(lgd)) %>%
        arrange(lgd) %>%
        mutate(key=chs[1:n()]) %>% select(key, network=lgd)
    #
    paths = unique(ti$path)
    if(length(lays)==1) lays=rep(lays, length(paths))
    pl = list()
    for(i in 1:length(paths)) {
        path = paths[i]; lay = lays[i]
    e1 = e0 %>% filter(path == !!path) %>% select(-path)
    e2 = e1 %>%
        left_join(es, by=c('lgd'='network')) %>%
        replace_na(list(key = '')) %>%
        group_by(from, to) %>%
        summarise(keys=str_c(key, collapse=' ')) %>% ungroup() %>%
        mutate(keys = str_replace(keys, '^ +%', ''))
    tfs = unique(e1$from)
    g = as_tbl_graph(e2)
    p1 = ggraph(g, layout=lay) +
        geom_edge_link(aes(label=keys), angle_calc=angle,
                       label_dodge = unit(.2, 'cm'),
                       arrow=arrow(length=unit(.2,'cm'), angle=20, type='closed'),
                       lwd=.1, label_size=2.5,
                       start_cap=circle(.4,'cm'), end_cap=circle(.4,'cm')) +
        geom_node_text(aes(label=name, color=name %in% tfs), vjust=.4, size=3, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        otheme(panel.border = F) +
        guides(color=F, shape=T)
    pg = ggplot_gtable(ggplot_build(p1))
    pl[[path]] = pg
    }
    list(pl=pl, tbl=es)
    #}}}
}
plotNets2 <- function(ti, th, ta=symb, nc=3, lay='kk', angle='along') {
    #{{{
    x = ti %>% filter(nid != '') %>%
        left_join(ta, by=c('reg.gid'='gid')) %>% rename(reg=symbol) %>%
        mutate(reg = ifelse(is.na(reg), trim_gid(reg.gid), reg)) %>%
        left_join(ta, by=c('tgt.gid'='gid')) %>% rename(tgt=symbol) %>%
        mutate(tgt = ifelse(is.na(tgt), trim_gid(tgt.gid), tgt)) %>%
        mutate(reg = str_to_sentence(reg), tgt = str_to_sentence(tgt)) %>%
        inner_join(th, by='nid') %>%
        mutate(lgd = factor(lgd, levels=th$lgd))
    e0 = x %>% select(path, from=reg, to=tgt, lgd, score)
    #
    paths = unique(ti$path)
    pl = list()
    for (path in paths) {
    e1 = e0 %>% filter(path == !!path) %>% select(-path)
    e2 = e1 %>% group_by(from, to) %>%
        summarise(lgd = 'All_combined', score=max(score)) %>% ungroup()
    tfs = unique(e1$from)
    #
    g = as_tbl_graph(rbind(e1, e2))
    p1 = ggraph(g, layout=lay) +
        #geom_node_point(aes(color=name %in% at.tfs), alpha=.2) +
        geom_edge_link(arrow=arrow(length=unit(.1,'cm'), angle=20, type='closed'), lwd=.1, start_cap=circle(.25,'cm'), end_cap=circle(.25,'cm')) +
        geom_node_text(aes(label=name, color=name %in% tfs), vjust=.4, size=2.5, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        facet_edges(~lgd, ncol=nc) +
        otheme(margin=c(0,0,0,0)) +
        guides(color=F)
    pg = ggplot_gtable(ggplot_build(p1))
    i = 1
    r = floor((i-1)/nc) + 1; c = (i-1) %% nc + 1
    strip = sprintf("strip-t-%d-%d", c, r)
    idx = grep(pattern=strip, pg$layout$name)
    pg$grobs[[idx]]$grobs[[1]]$children[[1]] = rectGrob(gp=gpar(fill='yellow',lty=0))
    pg$grobs[[idx]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col='black'
    pl[[path]] = pg
    }
    pl
    #}}}
}
#}}}

#{{{ At interactions
#{{{ read, process and write
fi = file.path(dirw, 'Regulations_in_ATRM.xlsx')
ti = read_xlsx(fi) %>% select(-note)
fm = '~/projects/genome/data2/ortholog/10.maize.arath.tsv'
tm = read_tsv(fm, col_types='ccccc') %>% select(-fam) %>%
    mutate(tag = str_c(code, tag, sep='_')) %>% select(-code)

tp = ti %>% select(at.rid=tf,at.tid=target,at.tf=tf_name,at.tgt=target_name,mode) %>%
    inner_join(tm,by=c('at.rid'='gid2')) %>% rename(zm.rid=gid1,tag1=tag) %>%
    inner_join(tm,by=c('at.tid'='gid2')) %>% rename(zm.tid=gid1,tag2=tag) %>%
    arrange(at.rid, at.tid)

to = tp %>% inner_join(tx, by=c('zm.rid'='reg.gid','zm.tid'='tgt.gid')) %>%
    group_by(at.tf,at.tgt,at.rid,at.tid,mode,zm.rid,zm.tid,tag1,tag2) %>%
    nest() %>%
    mutate(nnet = map_int(data, nrow)) %>% ungroup()
to %>% filter(nnet>0) %>% select(at.tf,at.tgt,mode,zm.rid,zm.tid) %>% print(n=80)

fo = file.path(dirw, '08.at.mapped.to.maize.tsv')
write_tsv(tp, fo)

fo = file.path(dirw, '10.at.maize.rds')
saveRDS(to, fo)
#}}}

#{{{ permutation
gids = gcfg$gene %>% filter(ttype=='mRNA') %>% pull(gid)
rids = tp$zm.rid
rids = tx %>% distinct(reg.gid) %>% pull(reg.gid)
tids = tx %>% distinct(tgt.gid) %>% pull(tgt.gid)
psample <- function(i, rids, tids, size=100)
    tibble(reg.gid=sample(rids,size), tgt.gid=sample(tids,size))

tm = tibble(perm=1:1000) %>%
    mutate(data=map(perm, psample, rids=rids, tids=tids, size=nrow(tp))) %>%
    unnest(data) %>%
    inner_join(tx, by=c('reg.gid','tgt.gid')) %>%
    distinct(perm, reg.gid, tgt.gid) %>%
    count(perm)
summary(tm$n)

fo = file.path(dirw, '19.perm.tsv')
write_tsv(tm, fo)
#}}}

#{{{ plot permutation
fi = file.path(dirw, '19.perm.tsv')
tm = read_tsv(fi)

tms = tibble(x=75,y=0,label='observed')
pm = ggplot(tm) +
    geom_histogram(aes(x=n), binwidth=2, color='black', fill='white') +
    geom_point(data=tibble(x=median(tm$n),y=0-3), aes(x=x, y=y), color='blue', shape=17, size=2) +
    geom_point(data=tms, aes(x=x, y=y), color='red') +
    scale_x_continuous(name = '# TF-target interactions supported by GRN') +
    scale_y_continuous(name = 'number of permutation runs') +
    geom_vline(xintercept = tms$x, linetype='dashed', color='red') +
    geom_label_repel(data=tms, aes(x=x,y=y,label=label),
                     direction='both', nudge_x=-15,nudge_y=20) +
    otheme(legend.pos = 'none', legend.dir = 'v', legend.title = T,
        xtick=T, ytick=T,
        xtitle=T, ytitle=T, xtext=T, ytext=T)
fo = file.path(dirw, '19.perm.pdf')
#ggsave(pm, file=fo, width=4, height=4)
#}}}

#{{{ ABI & HY5
#{{{ read in, build
fi = file.path(dirw, '08.at.mapped.to.maize.tsv')
ti1 = read_tsv(fi)
fi = file.path(dirw, '10.at.maize.rds')
ti2 = readRDS(fi)
build_atrm_net <- function(path, at.tfs, ti1, ti2) {
    #{{{
    x0 = ti2 %>% filter(at.tf %in% tfs, at.tf != at.tgt) %>%
        arrange(at.tf) %>% select(at.tf, at.tgt, data) %>%
        unnest(data) %>% select(reg.gid=at.tf, tgt.gid=at.tgt, nid, score)
    x = ti1 %>% filter(at.tf %in% at.tfs) %>% mutate(path=!!path) %>%
        select(path, reg.gid=at.tf, tgt.gid=at.tgt) %>%
        left_join(x0, by=c('reg.gid','tgt.gid')) %>%
        replace_na(list(nid = '')) %>% arrange(reg.gid, tgt.gid)
    #}}}
}
#}}}

path='aba'; tfs=c('ABI3','ABI4','ABI5'); x1=build_atrm_net(path, tfs, ti1, ti2)
path='hy5'; tfs=c('HY5'); x2=build_atrm_net(path, tfs, ti1, ti2)
x = rbind(x1, x2)

res = plotNets1(x, th, symb, 'circle')
pt = plotTable(res$tbl, nr=2, size=8)
fo = sprintf("%s/22.pdf", dirw)
ggarrange(pm, res$pl$hy5, res$pl$aba, pt,
        widths=c(1.2,1), heights=c(1,1), labels=c('A','C','B'),
    nrow=2, ncol=2) %>%
    ggexport(filename = fo, width=8, height=7.3)
#
p2 = plotNets2(x, th, symb, 5)
fo = sprintf("%s/22b.pdf", dirw)
ggarrange(p2$aba, p2$hy5,
        widths=c(1,1), heights=c(3,4), labels=c('A','B'),
    nrow=2, ncol=1) %>%
    ggexport(filename = fo, width=8, height=10)

#{{{ #[obsolete] heatmap
xs = x %>% distinct(lgd, col) %>% arrange(lgd)
swit = (max(x$score)+min(x$score)) / 2
p1 = ggplot(x, aes(x=lgd, y=pair, fill=score)) +
    geom_tile() +
    geom_text(aes(label=lab, color=score>swit), hjust=.5, size=2) +
    scale_x_discrete(expand=expand_scale(mult=c(0,0))) +
    scale_y_discrete(limits = rev(levels(as.factor(x$pair))), expand=c(0,0)) +
    scale_fill_viridis(name='feature importance score', direction=-1) +
    scale_color_manual(values=c('black','white')) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           ygrid=F, xtick=T, ytick=F, xtitle=F, xtext=T, ytext=T) +
    theme(axis.text.x = element_text(angle=30, color=rev(xs$col), hjust=1, vjust=1, size=7)) +
    guides(color=F)
fo = sprintf("%s/21.%s.pdf", dirw, tit)
ggsave(p1, file=fo, width=wd, height=ht)
#}}}
#}}}
#}}}

#{{{ CornCyc examples (pathway plot)
#{{{ read in
fv = sprintf("%s/rf.go.rds", dirr)
ev_go = readRDS(fv)
ev = ev_go %>% select(nid, enrich_grp) %>% unnest(enrich_grp) %>%
    filter(ctag=='CornCyc', pval < .01, score==10)
#
fi = '~/projects/genome/data/Zmays_B73/61_functional/07.corncyc.rds'
cc0 = readRDS(fi)
cc = gs$fun_ann %>% filter(ctag=='CornCyc') %>% select(-ctag)
#}}}
#{{{ explore
tz = ev %>% count(grp, max.reg.gid) %>% filter(n>2) %>% select(grp,max.reg.gid)
ev %>% inner_join(tz, by=c('grp','max.reg.gid')) %>% select(grp,nid,fc,pval,max.reg.gid,max.reg.size) %>%
    arrange(grp,max.reg.gid) %>% print(n=100)
#}}}

cca = read_xlsx(fp, 'cc2') %>% fill(pathway) %>%
    mutate(reactants=map(reactants, str_split_s, sep='~')) %>%
    mutate(products=map(products, str_split_s, sep='~')) %>%
    mutate(gids=map(gids, str_split_s, sep=' '))
cc9 = rbind(cc0, cca)

tp2 = read_xlsx(fp, 'cc')
idxs = 1:5
pn = list()
for (i in idxs) {
    path = tp2$pathway[i]
    extra_filt = str_split(tp2$filt[i], ',') %>% unlist()
    pn1 = plotNet_cc(path, cc9, symb, extra_filt, tp2$lay[i], tp2$angle[i])
    pn[[tp2$path[i]]] = pn1
fo = sprintf("%s/x.pdf", dirw)
#ggsave(pn, file=fo, width=5, height=5)
}

fo = sprintf("%s/24.cc.pdf", dirw)
ggarrange(pn[[1]], pn[[2]], pn[[3]], pn[[4]], pn[[5]],
    nrow=1, ncol=5, widths=c(1,1,1,1,1)) %>%
    ggexport(filename = fo, width=8, height=6)
#}}}

#{{{ known & cc pathways
tp1 = read_xlsx(fp, 'path') %>% mutate(tfs = map(tf, str_split_s, ' ')) %>%
    mutate(tgts = map(tgt, str_split_s, ' '))
makeNets <- function(idxs, tp, tx, cc) {
    #{{{
    x = tibble()
    for (i in idxs) {
        opt = tp$opt[i]
        path=tp$path[i]; tfs=tp$tfs[[i]]; tgts=tp$tgts[[i]]
        if(opt == 'cc')
            tgts = cc %>% filter(grp == tp$grp[i]) %>% pull(gid)
        x1 = tx %>% filter(reg.gid %in% tfs, tgt.gid %in% tgts) %>%
            mutate(path=!!path) %>% select(path, everything())
        x = rbind(x, x1)
    }
    x
    #}}}
}

#{{{ fig 4
idxs = c(1,2,8); lays = tp1$lay1[idxs]
x = makeNets(idxs, tp1, tx, cc)
#
res = plotNets1(x, th, symb, lays)
pt = plotTable(res$tbl, nr=1, size=8)
fo = sprintf("%s/25.pdf", dirw)
ggarrange(
    ggarrange(pn$antho, res$pl$antho,
        pn$dimboa, res$pl$dimboa,
        pn$hb26, res$pl$hb26,
        widths=c(1,1.3), heights=c(1,1,1), labels=LETTERS[1:6],
        nrow=3, ncol=2),
    pt, widths=c(4,1), nrow=1, ncol=2) %>%
    ggexport(filename = fo, width=8.5, height=10)

p2 = plotNets2(x, th, symb, 5)
for(i in idxs) {
    path=tp1$path[i]; wd2=tp1$wd2[i]; ht2=tp1$ht2[i]
    fo = sprintf("%s/25.%s.pdf", dirw, path)
    ggsave(p2[[path]], file=fo, width=wd2, height=ht2)
}

fo = sprintf("%s/25b.pdf", dirw)
ggarrange(p2$antho, p2$dimboa, p2$hb26,
        widths=c(1,1), heights=c(3,3,4), labels=LETTERS,
    nrow=3, ncol=1) %>%
    ggexport(filename = fo, width=9, height=11)
#}}}

#{{{ fig 4 - sup
idxs = c(3,4,6,7); lays = tp1$lay1[idxs]
x = makeNets(idxs, tp1, tx, cc)
#
res = plotNets1(x, th, symb, lays)
pt = plotTable(res$tbl, nr=1, size=8)
fo = sprintf("%s/26.pdf", dirw)
ggarrange(
    ggarrange(res$pl$tb1, res$pl$glu,
        pn$wrky34, res$pl$wrky34,
        pn$col11, res$pl$col11,
        widths=c(1,1.2), heights=c(1,1,1), labels=LETTERS[1:6],
        nrow=3, ncol=2),
    pt, widths=c(4,1), nrow=1, ncol=2, labels=c('','G')) %>%
    ggexport(filename = fo, width=8, height=10)
#
p2 = plotNets2(x, th, symb, 6)
fo = sprintf("%s/26b.pdf", dirw)
ggarrange(p2$tb1, p2$glu, p2$wrky34, p2$col11,
        widths=c(1,1), heights=c(1,2,1,1), labels=LETTERS,
    nrow=4, ncol=1) %>%
    ggexport(filename = fo, width=8, height=10)
#}}}
#}}}




