source("functions.R")
diri = '~/projects/rnaseq'
dirw = file.path(dird, '16_pathways')
th = t_cfg %>% mutate(lgd=str_replace(lgd, " \\[.*\\]", '')) %>% select(nid,lgd,col)

fv = sprintf("%s/rf.100k.rds", dirr)
ev = readRDS(fv)
tx = ev %>% select(nid, tn) %>% unnest()

fg = file.path(dirw, 'pathways.xlsx')
tg = read_xlsx(fg) %>% left_join(gcfg$gene[,c('gid','note2')]) %>%
    mutate(name = ifelse(is.na(name), gid, name)) %>%
    replace_na(list(role='tgt')) %>%
    print(n=50)
tf = read_tf_info()

#{{{ At interactions
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
    mutate(nnet = map_int(data, nrow))
to %>% filter(nnet>0) %>% select(at.tf,at.tgt,mode,zm.rid,zm.tid) %>% print(n=80)

#{{{ write
fo = file.path(dirw, '08.at.mapped.to.maize.tsv')
write_tsv(tp, fo)

to1 = to %>% select(-data)
fo = file.path(dirw, '10.at.maize.tsv')
write_tsv(to1, fo)
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
    unnest() %>%
    inner_join(tx, by=c('reg.gid','tgt.gid')) %>%
    distinct(perm, reg.gid, tgt.gid) %>%
    count(perm)
summary(tm$n)

tms = tibble(x=75,y=0,label='observed')
p1 = ggplot(tm) +
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
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ ABI & HY5
tit = 'hy'; at.tfs = c("HY5"); wd=6; ht=3
tit = 'aba'; at.tfs = c("ABI3","ABI2","ABI4","ABI5"); wd=5; ht=4
to %>% filter(at.tf %in% at.tfs) %>%
    arrange(at.tgt) %>% print(n=30)

x = to %>% filter(at.tf %in% at.tfs, at.tf != at.tgt) %>%
    arrange(at.tf) %>% select(at.tf, at.tgt, data) %>%
    unnest() %>% mutate(lab = sprintf("%.02f\n%.0e", pcc, pval.p)) %>%
    inner_join(th, by='nid') %>%
    mutate(lgd = factor(lgd, levels=th$lgd)) %>%
    mutate(pair = str_c(at.tf, at.tgt, sep=' -> ')) %>%
    arrange(at.tf, at.tgt)
#    select(at.tf, at.tgt, lgd, lab) %>%
#    spread(lgd, lab)

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

fo = file.path(dirw, 'x.tsv')
write_tsv(x, fo, na='')
#}}}

#{{{ known pathways
path = 'antho'; wd=6; ht=4
path = 'dimboa'; wd=8; ht=8
path = 'tb1'; wd=4; ht=3
tg0 = tg %>% filter(pathway == path)
gidsR = tg0 %>% filter(role=='tf') %>% pull(gid)
gidsT = tg0 %>% filter(role=='tgt') %>% pull(gid)
#
x = tx %>% filter(reg.gid %in% gidsR, tgt.gid %in% gidsT) %>%
    inner_join(tg0[,c('gid','name')], by=c('reg.gid'='gid')) %>% rename(tf=name) %>%
    inner_join(tg0[,c('gid','name')], by=c('tgt.gid'='gid')) %>% rename(tgt=name) %>%
    unnest() %>% mutate(lab = sprintf("%.02f\n%.0e", pcc, pval.p)) %>%
    inner_join(th, by='nid') %>%
    mutate(lgd = factor(lgd, levels=th$lgd)) %>%
    mutate(pair = str_c(tf, tgt, sep=' -> ')) %>%
    arrange(tf,tgt)

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
fo = sprintf("%s/22.%s.pdf", dirw, path)
ggsave(p1, file=fo, width=wd, height=ht)
#}}}


#{{{ get DIMBOA genes - obsolete
tg = gcfg$gene %>% select(gid, chrom, start, end, note1, note2) %>%
    filter(str_detect(note1, 'benzoxazinone synthesis') |
        str_detect(note2, 'benzoxazinone synthesis') |
        str_detect(note2, 'benzoxazinless')) %>%
    mutate(note2 = str_replace(note1, 'benzoxazinone synthesis', '')) %>%
    mutate(note2 = str_replace(note2, 'benzoxazinless', '')) %>%
    mutate(note2 = as.integer(note2)) %>% arrange(note2) %>%
    mutate(note1 = sprintf("bx%d", note2))
#}}}




