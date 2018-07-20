#{{{ head
source("grn.fun.r")
#source("enrich.R")
fi = file.path(dird, '05.previous.grns/10.RData')
x = load(fi)
fi = file.path(dird, '07.known.tf/10.RData')
x = load(fi)
# read RNA-Seq data
dirw = file.path("~/projects/briggs/data", "49.coop")
fi = file.path(dirw, "01.master.RData")
x = load(fi)
fd = file.path(dirw, "03.sharing.RData")
x = load(fd)
dirw = file.path(dira, "06_rna_de")
#}}}

#{{{ TF
ff = '/home/springer/zhoux379/data/genome/Zmays_v4/TF/11.tsv'
tf = read_tsv(ff)

to = tt %>% 
    left_join(tf, by = 'gid') %>%
    group_by(ctag, tag, tsTag) %>%
    summarise(ngene.total = n(),
              n.tf = sum(!is.na(fam)),
              prop.tf = n.tf / ngene.total) %>%
    print(n=50)

tt %>% 
    inner_join(tf, by = 'gid') %>%
    filter(ctag == 'HC') %>%
    group_by(tag, fam) %>%
    summarise(ng = n()) %>%
    arrange(tag, -ng) %>%
    print(n=50)
#}}}

#{{{ remake tm w.o. using DE 2-FC cutoff
fi = file.path(dirp, "41.qc/10.rep.merged.RData")
x = load(fi)
nrow(t_rep_merged)/69
fi = file.path(dirp, "42.de/11.de.RData")
x = load(fi)
fi = file.path(dirp, "44.ase/10.ase.2.RData")
x = load(fi)
#fi = file.path(dirp, "45.doa/10.dom.RData")
#x = load(fi)
t_de = t_de %>% mutate(Tissue = factor(Tissue, levels = tissues23))
t_ase = t_ase %>% mutate(Tissue = factor(Tissue, levels = tissues23))
#t_dom = t_dom %>% mutate(Tissue = factor(Tissue, levels = tissues23)) %>% select(Tissue, gid, Dom, DoA)
taglst = list(
    pDE = c("DE_B", "DE_M", "non_DE"),
    hDE = levels(t_de$hDE),
    Dom = levels(t_de$Dom),
    Reg1 = levels(t_ase$Reg1),
    Reg2 = levels(t_ase$Reg2)
)

tm = t_rep_merged %>% filter(Genotype != 'MxB') %>%
    select(Tissue, gid, Genotype, CPM) %>%
    spread(Genotype, CPM)
nrow(tm)/23
tm2 = tm %>%
    left_join(t_de, by = c("Tissue", 'gid')) %>%
    #mutate(pDE = pDE2) %>%  #!!! use 2-FC DE cutoff
    mutate(silent = ifelse(is.na(pDE), 1, 0)) %>%
    select(Tissue, gid, B73, Mo17, BxM, silent, log2MB, pDE, log2FM, hDE, Dom, DoA)
tm2 %>% filter(is.nan(DoA) | is.infinite(DoA))
tm2 %>% count(pDE, Dom) %>% spread(Dom, n)
tm2 %>% count(pDE, hDE) %>% spread(hDE, n)

tm3 = tm2 %>% left_join(t_ase, by = c("Tissue", 'gid')) 
nrow(tm3)/23
tm3 %>% count(pDE, Reg1) %>% spread(Reg1, n)
tm3 %>% count(hDE, Reg2) %>% spread(Reg2, n)
tm3 = tm3 %>% 
    mutate(Reg1 = as.character(Reg1),
           Reg2 = as.character(Reg2)) %>%
    mutate(Reg1 = ifelse(!is.na(pDE) & pDE != 'non_DE', Reg1, NA),
           Reg2 = ifelse(!is.na(pDE) & pDE == 'non_DE', Reg2, NA)) %>%
    mutate(Reg1 = factor(Reg1, levels = taglst$Reg1),
           Reg2 = factor(Reg2, levels = taglst$Reg2))
tm3 %>% count(pDE, Reg1) %>% spread(pDE, n)
tm3 %>% count(pDE, Reg2) %>% spread(pDE, n)

tm4 = tm3 %>% 
    mutate(MP = (B73 + Mo17) / 2,
           SPE = ifelse(silent == 1, NA,
                 ifelse(B73>=1 & Mo17<0.1 & pDE=='DE_B', 'SPE_B',
                 ifelse(B73<0.1 & Mo17>=1 & pDE=='DE_M', 'SPE_M', 'non_SPE'))),
           HC = ifelse(is.na(SPE), NA,
                ifelse(SPE=='SPE_B' & (BxM>=1 | BxM >= MP), 'HC_B',
                ifelse(SPE=='SPE_M' & (BxM>=1 | BxM >= MP), 'HC_M', 'non_HC')))) %>%
    select(-MP)
tm4 %>% count(Tissue, SPE) %>% spread(SPE, n) %>% print(n=23)
tm = tm4
tm %>% count(Tissue, pDE) %>% spread(pDE, n) %>% print(n=23)
#}}}

#{{{ GRN
require(network)
ctags = t_grn %>% distinct(ctag) %>% pull(ctag)
ctags

for (grn in ctags) {
    #{{{
tn = t_grn %>% filter(ctag == grn) %>% select(-ctag) %>% 
    top_n(10000, score) %>%
    transmute(r.gid = regulator, t.gid = target, r.fam = fam) 
#
tm.reg = tm %>% filter(!is.na(pDE)) %>%
    transmute(Tissue = Tissue, r.gid = gid, r.DE = pDE, r.log2MB = log2MB, r.SPE = SPE)
tm.tgt = tm %>% filter(!is.na(pDE)) %>%
    transmute(Tissue = Tissue, t.gid = gid, t.DE = pDE, t.log2MB = log2MB, r.reg = Reg1)
#
tags = c("non_DE", "DE_B|2-4", "DE_M|2-4", "DE_B|4+", "DE_M|4+", "SPE_B", "SPE_M")
tn1 = tn %>% distinct(r.gid, r.fam) %>% inner_join(tm.reg, by = 'r.gid')
tn1a = tn1 %>% filter(r.DE == 'non_DE') %>% mutate(tag = r.DE) %>%
    select(Tissue, tag, r.gid)
tn1b = tn1 %>% filter(r.DE != 'non_DE') %>% 
    mutate(tag = ifelse(abs(r.log2MB)<2, '2-4', '4+')) %>%
    mutate(tag = sprintf("%s|%s", r.DE, tag)) %>%
    select(Tissue, tag, r.gid)
tn1c = tn1 %>% filter(r.SPE != 'non_SPE') %>% mutate(tag = r.SPE) %>%
    select(Tissue, tag, r.gid)
tn = rbind(tn1a, tn1b, tn1c) %>%
    mutate(tag = factor(tag, levels = rev(tags))) %>%
    inner_join(tn, by = 'r.gid') %>%
    inner_join(tm.tgt, by = c('Tissue', 't.gid')) 
#
tps = tn %>% group_by(Tissue, tag) %>%
    summarise(n.reg = length(unique(r.gid)),
              n.tgt = length(unique(t.gid))) %>%
    mutate(lab = sprintf("%d regulating %d", n.reg, n.tgt))
tp = tn %>% distinct(Tissue, tag, t.gid, t.DE, t.log2MB, t.reg) %>%
    group_by(Tissue, tag, t.DE) %>%
    summarise(n = n()) %>% mutate(p = n/sum(n))
tp %>%
    select(-n) %>% spread(t.DE, p) %>%
    print(n=69)
#
pa = ggplot(tp) +
    geom_bar(aes(x = tag, y = p, fill = t.DE), position = position_stack(reverse = T), stat = 'identity', width = .8) +
    geom_text(data = tps, aes(x = tag, y = 1-.03, label = lab), color = 'white', size = 3, hjust = 1) +
    scale_x_discrete(name = 'Regulator', expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Prop. Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_jco(name = 'Target') +
    coord_flip() +
    facet_wrap(~Tissue, nrow = 6, strip.position = 'top') +
    theme_bw() +
    theme(legend.position=c(.5,1), legend.justification=c(.5,-.7), legend.direction = 'horizontal', legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
#
fo = sprintf("%s/01.%s.pdf", dirw, grn)
ggarrange(pa, nrow = 1, ncol = 1)  %>% 
    ggexport(filename = fo, width = 10, height = 10)
    #}}}
}

grn = 'huang_sam'
tissue = 'seedlingmeristem_11DAS'
#grn = 'huang_root'
#tissue = 'seedlingroot_11DAS'
tn0 = t_grn %>% filter(ctag == grn) %>% select(-ctag) %>% 
    transmute(r.gid = regulator, t.gid = target, r.fam = fam, score = score) %>%
    top_n(20000, score)
#
tm.reg = tm %>% filter(!is.na(pDE)) %>%
    transmute(Tissue = Tissue, r.gid = gid, r.DE = pDE, r.log2MB = log2MB, r.SPE = SPE)
tm.tgt = tm %>% filter(!is.na(pDE)) %>%
    transmute(Tissue = Tissue, t.gid = gid, t.DE = pDE, t.log2MB = log2MB, r.reg = Reg1)
#
ctags = c("non_DE", "DE_B", "DE_M", "non_SPE", "SPE_B", "SPE_M")
tn1 = tn0 %>% distinct(r.gid, r.fam) %>% inner_join(tm.reg, by = 'r.gid') 
tn2 = tn1 %>% transmute(Tissue = Tissue, ctag = r.DE, r.gid = r.gid)
tn3 = tn1 %>% transmute(Tissue = Tissue, ctag = r.SPE, r.gid = r.gid)
tn = rbind(tn2, tn3) %>%
    mutate(ctag = factor(ctag, levels = rev(ctags))) %>%
    inner_join(tn0, by = 'r.gid') %>%
    inner_join(tm.tgt, by = c('Tissue', 't.gid')) 
#
tn.spe = tn %>% filter(ctag %in% c("DE_B", "DE_M")) %>% 
    group_by(Tissue, r.gid) %>%
    summarise(n.tgt = n(), 
              n.de.b = sum(t.DE == 'DE_B'),
              n.de.m = sum(t.DE == 'DE_M'),
              p.de.b = sum(t.DE == 'DE_B')/n.tgt,
              p.de.m = sum(t.DE == 'DE_M')/n.tgt) %>%
    filter(n.tgt >= 3, p.de.b > .7 | p.de.m > .7) %>%
    filter(Tissue == tissue) %>%
    print(n = 30)

i = 1
tnz = tn %>% filter(ctag %in% c("DE_B", "DE_M"), 
                    Tissue == tn.spe$Tissue[i], 
                    r.gid == tn.spe$r.gid[i]) %>% select(-Tissue) %>% print(n=60)
#tq %>% filter(gid %in% tnz$t.gid) %>% select(gid, ADD, qchr, qpos, type) %>%
#    print(n=50)
rid = tn.spe$r.gid[i]; tids = tnz %>% pull(t.gid)
vids = 1:(1+length(tids))
names(vids) = c(rid, tids)
tno = tn0 %>% filter(r.gid %in% tm.reg$r.gid, t.gid %in% tm.tgt$t.gid,
                     (r.gid == rid & t.gid %in% tids) | (r.gid %in% tids & t.gid %in% tids)) %>%
    transmute(tails = vids[r.gid], heads = vids[t.gid], 
              r.gid = r.gid, t.gid = t.gid, score = score) %>%
    print(n=50)

gn0 = as.data.frame(tno)
gn = network.initialize(length(vids))
gn = network.edgelist(gn0, gn, directed = T)
network.vertex.names(gn) = names(vids)
gn %v% "color" = ifelse(names(vids) == rid, 'tomato', 'steelblue')
gn %e% "score" = tno$score * 10
pn = ggnet2(gn, size = 10, alpha = 0, label = T, label.size = 3, color = 'color',
            arrow.size = 12, arrow.gap = 0.025,
            edge.size = 'score') +
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
fo = file.path(dirw, 'test3.pdf')
ggarrange(pn, NULL, 
    nrow = 1, ncol = 2)  %>% 
    ggexport(filename = fo, width = 8, height = 4)
#}}}

#{{{ eQTL & GRN
fg = '~/data/genome/B73/v32/t5.gtb'
tg = read_tsv(fg) %>% transmute(gid = par, chrom = chr, start = beg, end = end, srd = srd)

fq = '~/data/misc1/li2013/10.eQTL.v4.tsv'
tq = read_tsv(fq)
tql = tq %>%
    transmute(t.gid = gid, r.chr = qchr, r.pos = qpos, reg = type, addi = ADD) %>%
    filter(r.chr %in% 1:10) %>%
    mutate(r.chr = as.integer(r.chr)) %>%
    arrange(r.chr, r.pos)
tql %>% count(reg)
tqlo = tql %>% transmute(chrom = r.chr, 
                         start = r.pos - 1, end = r.pos, 
                         t.gid = t.gid, reg = reg)

#
grn = 'huang_sam'
tn = t_grn %>% filter(ctag == grn) %>% select(-ctag) %>% 
    top_n(100000, score) %>%
    transmute(r.gid = regulator, t.gid = target, r.fam = fam)
tno = tn %>% distinct(r.gid) %>% inner_join(tg, by = c('r.gid'='gid')) %>%
    arrange(chrom, start) %>%
    transmute(chrom = chrom, 
              start = pmax(0, start - 10000001), 
              end = end + 10000000, gid = r.gid)

fo = file.path(dirw, 't1.tf.bed')
write_tsv(tno, fo, col_names = F)
fo = file.path(dirw, 't2.eqtl.bed')
write_tsv(tqlo, fo, col_names = F)
system("intersectBed -a t1.tf.bed -b t2.eqtl.bed -wo > t3.bed")

fx = file.path(dirw, 't3.bed')
tx = read_tsv(fx, col_names = F)[,c(4,5,7,8,9)]
colnames(tx) = c('r.gid', 'r.chr', 'r.pos', 't.gid', 'reg')
tx %>% distinct(r.gid)
tx %>% distinct(t.gid)
tn2 = tn %>% left_join(tx, by = c('r.gid','t.gid'))
tn2 %>% count(reg)
#}}}


