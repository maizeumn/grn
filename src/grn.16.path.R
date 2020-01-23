source("functions.R")
require(ggpmisc)
require(ggraph)
require(tidygraph)
diri = '~/projects/rnaseq'
dirw = file.path(dird, '16_pathways')
th = t_cfg %>% mutate(lgd=str_replace(lgd, " \\[.*\\]", '')) %>% select(nid,lgd,col)
tf = read_tf_info()
fp = file.path(dirw, 'pathways.xlsx')
#
fv = sprintf("%s/rf.100k.rds", dirr)
tx = readRDS(fv) %>% select(nid, tn) %>% unnest(tn)

#{{{ plot functions
make_cc_step <- function(reacs, prods, gids, spon, ta=symb, extra_filt=c()) {
    #{{{
    known_filt = c('ATP','CTP','AMP','CMP','ADP','UDP','NADPH','NADP+',
                   'H2O','H+','oxygen','CO2',
                   'Mg2+','phosphate','diphosphate','pyruvate','fumarate')
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
plot_cc_net <- function(path, cc0, ta=symb, extra_filt=c(), lay='kk', angle='along') {
    #{{{
    cc2 = cc0 %>% filter(pathway==!!path) %>%
        mutate(net = pmap(list(reactants,products,gids,spontaneous), make_cc_step, symb, extra_filt=!!extra_filt)) %>%
        select(net) %>% unnest(net)
    # ad hoc adjustments
    cc2 = cc2 %>% filter(snode != 'isopentenyl diphosphate' | enode != 'dimethylallyl diphosphate')
    #
    e0 = cc2 %>% select(from=snode, to=enode, lab)
    g = as_tbl_graph(e0)
    circular = ifelse(lay == 'linear', T, F)
    pg = ggraph(g, layout=lay, circular=!!circular) +
        geom_edge_link2(aes(label = lab), angle_calc=angle,
                       #label_dodge = unit(.2, 'cm'),
                       arrow=arrow(length=unit(.2,'cm'), angle=20, type='closed'),
                       lwd=.1, label_size=2.5, label_colour='darkblue',
                       start_cap=circle(.2,'cm'), end_cap=circle(.15,'cm')) +
        geom_node_text(aes(label=name), vjust=.4, size=2.5, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        otheme(panel.border=F, margin=c(0,0,0,0)) +
        guides(color=F, shape=T)
    pg
    #}}}
}
plot_network_1 <- function(ti, th, tfs, ta=symb, wr=.7, lay='kk') {
    #{{{
    x = ti %>%
        left_join(ta, by=c('reg.gid'='gid')) %>% rename(reg=symbol) %>%
        mutate(reg = ifelse(is.na(reg), trim_gid(reg.gid), reg)) %>%
        left_join(ta, by=c('tgt.gid'='gid')) %>% rename(tgt=symbol) %>%
        mutate(tgt = ifelse(is.na(tgt), trim_gid(tgt.gid), tgt)) %>%
        left_join(th, by='nid') %>%
        mutate(lgd = factor(lgd, levels=th$lgd))
    tfs = unique(c(tfs, x$reg))
    e0 = x %>% select(from=reg, to=tgt, lgd, score)
    es = e0 %>% distinct(lgd) %>% filter(!is.na(lgd)) %>%
        arrange(lgd) %>%
        mutate(key=letters[1:n()]) %>% select(key, network=lgd)
    est = tibble(x=.5,y=.5, tb=list(es))
    e2 = e0 %>% left_join(es, by=c('lgd'='network')) %>%
        replace_na(list(key = '')) %>%
        group_by(from, to) %>%
        summarise(keys=str_c(key, collapse=' ')) %>% ungroup() %>%
        mutate(keys = str_replace(keys, '^ +%', ''))
    g = as_tbl_graph(e2)
    pg = ggraph(g, layout=lay) +
        geom_edge_link(aes(label=keys), angle_calc='along',
                       label_dodge = unit(.2, 'cm'),
                       arrow=arrow(length=unit(.2,'cm'), angle=20, type='closed'),
                       lwd=.1, label_size=2.5,
                       start_cap=circle(.4,'cm'), end_cap=circle(.4,'cm')) +
        geom_node_text(aes(label=name, color=name %in% tfs), vjust=.4, size=3, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        otheme(panel.border = F) +
        guides(color=F, shape=T)
    pt = ggplot() +
        geom_table_npc(data=est, aes(npcx=x,npcy=y,label=tb), size=2.5, color='black', vjust=.5, hjust=.5) +
        otheme(margin=c(0,0,0,0), panel.border=F)
    list(pg=pg, pt=pt)
    #}}}
}
plot_network_2 <- function(ti, th, tfs, ta=symb, nc=3, lay='kk') {
    #{{{
    x = ti %>% filter(nid != '') %>%
        left_join(ta, by=c('reg.gid'='gid')) %>% rename(reg=symbol) %>%
        mutate(reg = ifelse(is.na(reg), trim_gid(reg.gid), reg)) %>%
        left_join(ta, by=c('tgt.gid'='gid')) %>% rename(tgt=symbol) %>%
        mutate(tgt = ifelse(is.na(tgt), trim_gid(tgt.gid), tgt)) %>%
        inner_join(th, by='nid') %>%
        mutate(lgd = factor(lgd, levels=th$lgd))
    tfs = unique(c(tfs, x$reg))
    e0 = x %>% select(from=reg, to=tgt, lgd, score)
    e1 = e0 %>% group_by(from, to) %>% summarise(lgd = 'All_combined', score=max(score)) %>% ungroup()
    #
    g = as_tbl_graph(rbind(e0, e1))
    pg = ggraph(g, layout=lay) +
        #geom_node_point(aes(color=name %in% at.tfs), alpha=.2) +
        geom_edge_link(arrow=arrow(length=unit(.2,'cm'), angle=20, type='closed'), lwd=.1, start_cap=circle(.25,'cm'), end_cap=circle(.25,'cm')) +
        geom_node_text(aes(label=name, color=name %in% tfs), vjust=.4, size=2.5, repel=T, point.padding=NA, box.padding=0, force=.1)+
        scale_color_aaas() +
        facet_edges(~lgd) +
        otheme() +
        guides(color=F)
    pg
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
fi = file.path(dirw, '08.at.mapped.to.maize.tsv')
ti1 = read_tsv(fi)
fi = file.path(dirw, '10.at.maize.rds')
ti2 = readRDS(fi)

tit = 'hy5'; at.tfs = c("HY5"); wr=.7; nc=4; wd=8; ht=8
tit = 'aba'; at.tfs = c("ABI2","ABI3","ABI4","ABI5"); wr=.7; nc=4; wd=8; ht=8
x0 = ti2 %>% filter(at.tf %in% at.tfs, at.tf != at.tgt) %>%
    arrange(at.tf) %>% select(at.tf, at.tgt, data) %>%
    unnest(data) %>% select(reg.gid=at.tf, tgt.gid=at.tgt, nid, score)
x = ti1 %>% filter(at.tf %in% at.tfs) %>%
    select(reg.gid=at.tf, tgt.gid=at.tgt) %>%
    left_join(x0, by=c('reg.gid','tgt.gid')) %>%
    replace_na(list(nid = '')) %>% arrange(reg.gid, tgt.gid)
#
p1 = plot_network_1(x, th, at.tfs, symb, wr)
p2 = plot_network_2(x, th, at.tfs, symb, nc, 'kk')
fo = sprintf("%s/21.%s_1.pdf", dirw, tit)
#ggarrange(p1$pg, p1$pt, nrow=1, ncol=2, widths=c(wr,1-wr)) %>%
#    ggexport(filename = fo, width=6, height=6)
fo = sprintf("%s/21.%s_2.pdf", dirw, tit)
ggsave(p2, file=fo, width=wd, height=ht)
#{{{ save for multi-panel plot
if(tit == 'aba') {
    p11 = p1$pg; p12 = p1$pt
} else if(tit == 'hy5') {
    p21 = p1$pg; p22 = p1$pt
}
#}}}

fo = sprintf("%s/22.pdf", dirw)
ggarrange(pm,
    ggarrange(
        ggarrange(p11, p12, ncol=2, widths=c(1,.5)),
        ggarrange(p21, p22, ncol=2, widths=c(1,.5)),
        widths=c(1.2,1), labels=c("B",'C')),
    nrow=2, ncol=1, heights=c(1,1.3), labels=LETTERS[1]) %>%
    ggexport(filename = fo, width=8, height=7.3)

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

#{{{ CornCyc examples
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

tp = read_xlsx(fp, 'cc')
i=3
reg.gid = tp$reg.gid[i]; grp = tp$grp[i]; tf = tp$tf[i]; path = tp$pathway[i]
extra_filt = str_split(tp$filt[i], ',') %>% unlist()
pn = plot_cc_net(path, cc0, symb, extra_filt, tp$lay[i], tp$angle[i])
fo = sprintf("%s/x.pdf", dirw)
ggsave(pn, file=fo, width=5, height=5)
nc=tp$nc[i]; wr=tp$wr[i]
wd1=tp$wd1[i]; ht1=tp$ht1[i]; wd2=tp$wd2[i]; ht2=tp$ht2[i]
tgt.gids = cc %>% filter(grp == !!grp) %>% pull(gid)
#
x = tx %>% filter(reg.gid==!!reg.gid, tgt.gid %in% tgt.gids)
p1 = plot_network_1(x, th, tf, symb, wr)
p2 = plot_network_2(x, th, tf, symb, nc)
fo = sprintf("%s/25.%s_1.pdf", dirw, tf)
ggarrange(p1$pg, p1$pt, nrow=1, ncol=2, widths=c(wr,1-wr)) %>%
    ggexport(filename = fo, width=wd1, height=ht1)
fo = sprintf("%s/25.%s_2.pdf", dirw, tf)
ggsave(p2, file=fo, width=wd2, height=ht2)
#{{{ save for multi-panel plot
if(i == 1) {
    p11 = pn; p12 = p1$pg; p13 = p1$pt
} else if(i == 2) {
    p21 = pn; p22 = p1$pg; p23 = p1$pt
} else if(i == 3) {
    p31 = pn; p32 = p1$pg; p33 = p1$pt
}
#}}}

fo = sprintf("%s/26.pdf", dirw)
ggarrange(
    ggarrange(p11, p12, p13, ncol=3, labels=c('A','B'), widths=c(1,1,.5)),
    ggarrange(p21, p22, p23, ncol=3, labels=c('C','D'), widths=c(1,1,.5)),
    ggarrange(p31, p32, p33, ncol=3, labels=c('E','F'), widths=c(1,1,.5)),
    nrow=3, ncol=1, heights=c(1,1.4,3)) %>%
    ggexport(filename = fo, width=8, height=10.5)
#}}}

#{{{ known pathways
str_split_s <- function(s, sep='') str_split(s, sep)[[1]]
tp = read_xlsx(fp, 'path') %>% mutate(tfs = map(tf, str_split_s, ' ')) %>%
    mutate(tgts = map(tgt, str_split_s, ' '))
i=1
path=tp$path[i]; tfs=tp$tfs[[i]]; tgts=tp$tgts[[i]]
nc=tp$nc[i]; wr=tp$wr[i]; lay1=tp$lay1[i]; lay2=tp$lay2[i]
wd1=tp$wd1[i]; ht1=tp$ht1[i]; wd2=tp$wd2[i]; ht2=tp$ht2[i]
x = tx %>% filter(reg.gid %in% tfs, tgt.gid %in% tgts)
#
p1 = plot_network_1(x, th, tfs, symb, wr, lay1)
p2 = plot_network_2(x, th, tfs, symb, nc, lay2)
fo = sprintf("%s/23.%s_1.pdf", dirw, path)
ggarrange(p1$pg, p1$pt, nrow=1, ncol=2, widths=c(wr,1-wr)) %>%
    ggexport(filename = fo, width=wd1, height=ht1)
fo = sprintf("%s/23.%s_2.pdf", dirw, path)
ggsave(p2, file=fo, width=wd2, height=ht2)
#{{{ save for multi-panel plot
if(i == 1) {
    p11 = p1$pg; p12 = p1$pt
} else if(i == 2) {
    p21 = p1$pg; p22 = p1$pt
} else if(i == 3) {
    p31 = p1$pg; p32 = p1$pt
} else if(i == 4) {
    p41 = p1$pg; p42 = p1$pt
}
#}}}

fo = sprintf("%s/24.pdf", dirw)
ggarrange(
    ggarrange(
        ggarrange(p11, p12, ncol=2, widths=c(1,.5)),
        ggarrange(p41, p42, ncol=2, widths=c(1,.5)),
        labels=c('A','B'), widths=c(1.4,1)),
    ggarrange(
        ggarrange(p21, p22, ncol=2, widths=c(1,.5)),
        ggarrange(p31, p32, ncol=2, widths=c(1,.5)),
        labels=c('C','D'), widths=c(1,1)),
    nrow=2, ncol=1, heights=c(1,1.7)) %>%
    ggexport(filename = fo, width=8, height=10)
#}}}



