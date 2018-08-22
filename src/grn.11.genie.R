#{{{ head
source("grn.fun.r")
require(PRROC)
#source("enrich.R")
diri = '~/projects/maize.expression'
dirw = '~/projects/maize.grn/data'
f_cfg = file.path(dirw, '10.genie3.tsv')
t_cfg = read_tsv(f_cfg)
studies = t_cfg %>% distinct(study) %>% pull(study)

fi = file.path(dird, '05.previous.grns/10.RData')
x = load(fi)
# read RNA-Seq data
dirw = file.path("~/projects/briggs/data", "49.coop")
fi = file.path(dirw, "01.master.RData")
x = load(fi)
fd = file.path(dirw, "03.sharing.RData")
x = load(fd)

# TFs
ff = '~/data/genome/Zmays_v4/TF/11.tsv'
tf = read_tsv(ff)
tf_ids = tf$gid
#}}}

#{{{ read known TF targets
fi = file.path(dirw, '07_known_tf/10.rda')
x = load(fi)
x
ctags = unique(t_gs$ctag)
reg.gids = unique(t_gs$reg.gid)
t_gss = t_gs %>% distinct(ctag, reg.gid)
t_gs = t_gs %>% mutate(binding = 1)
#}}}

#{{{ prepare genie3 input
#study = 'li2013'
#study = 'hirsch2014'
#study = 'leiboff2015'
#study = 'jin2016'
#study = 'stelpflug2016'
#study = 'walley2016'
#study = 'lin2017'
#study = 'kremling2018'
#study = 'briggs'
study = 'dev41'
#study = 'dev64'
study
diri1 = file.path(diri, study, 'data')
if(study %in% c("dev41","dev64")) diri1 = file.path(diri, 'data', study)
fi = file.path(diri1, '20.rc.norm.rda')
x = load(fi)
fh1 = file.path(diri1, '01.reads.tsv')
fh2 = file.path(diri1, '02.reads.corrected.tsv')
fh = ifelse(file.exists(fh2), fh2, fh1) 
th = read_tsv(fh)
diro = file.path(dird, '11_input')
ngene = 46117
#
if(study == 'li2013') {
    #{{{ li2013
    th = th %>% select(SampleID, Genotype, Replicate) %>%
        mutate(Replicate = sprintf("SAM%d", Replicate))
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = Replicate,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'hirsch2014') {
    #{{{ hirsch2014
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'seedling_503')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'leiboff2015') {
    #{{{ leiboff2015
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'SAM_380')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Genotype, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'jin2016') {
    #{{{
    th = th %>% select(SampleID, Genotype) %>%
        mutate(tag = 'kernel_368')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'stelpflug2016') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73_18')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'walley2016') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73_23')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'lin2017') {
    #{{{
    th = th %>% select(SampleID, Genotype, Tissue) %>%
        mutate(tag = Tissue)
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Genotype, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    t_exp2 = t_exp %>% 
        mutate(condition = sprintf("%s_%s", tag, condition),
               tag = '5_tissues')
    t_exp = t_exp %>% bind_rows(t_exp2)
    #}}}
} else if(study == 'kremling2018') {
    #{{{
    th = th %>% select(SampleID, Genotype, Tissue) %>%
        mutate(tag = Tissue)
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Genotype, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Genotype,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    t_exp2 = t_exp %>% 
        mutate(condition = sprintf("%s_%s", tag, condition),
               tag = '7_tissues')
    t_exp = t_exp %>% bind_rows(t_exp2)
    #}}}
} else if(study == 'briggs') {
    #{{{
    th = th %>% select(SampleID, Genotype, Tissue) %>%
        mutate(tag = Genotype)
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else if(study == 'dev41' | study == 'dev64') {
    #{{{
    th = th %>% select(SampleID, Tissue) %>%
        mutate(tag = 'B73')
    t_exp = tm %>% select(gid, SampleID, ReadCount, CPM, FPKM) %>% 
        inner_join(th, by = 'SampleID') %>%
        group_by(tag, Tissue, gid) %>%
        summarise(ReadCount = sum(ReadCount), CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup() %>%
        transmute(tag = tag,
                  condition = Tissue,
                  gid = gid, RC = ReadCount, CPM = CPM, FPKM = FPKM)
    #}}}
} else {
    cat('unknown study', study, "\n")
}
nrow(t_exp)/ngene
tt = t_exp
#
tags = unique(tt$tag)
for (tag in tags) {
    tc1 = t_cfg %>% filter(study == !!study, tag == !!tag)
    if(nrow(tc1) != 1)
        stop(sprintf("not 1 row: %s %s\n", study, tag))
    nid = tc1$nid[1]
    t_exp = tt %>% filter(tag == !!tag) %>% select(-tag)
    fo = sprintf("%s/%s.rda", diro, nid)
    save(t_exp, file = fo)
}
#}}}

th = t_cfg %>%
    mutate(fgn = sprintf("%s/12_output/%s.rda", dirw, nid),
    #mutate(fgn = sprintf("%s/12_output_cpm/%s.rda", dirw, nid),
           done = file.exists(fgn)) %>%
    filter(done) %>%
    select(-done)
#th = th %>% filter(nid %in% c('n41', 'n42', 'n71', 'n81', 'n82'))

#{{{ compute AUROC using known TF data
t_pr = tibble(); t_roc = tibble()
t_aupr = tibble(); t_auroc = tibble()
for (i in 1:nrow(th)) {
    nid = th$nid[i]
    fgn = th$fgn[i]
    #nid = 'n71'
    #fgn = th$fgn[th$nid==nid]
    x = load(fgn)
    x
    # 
    sum(rids %in% t_gs$reg.gid)
    tz = expand.grid(ctag = ctags, tgt.gid = tids, stringsAsFactors = F) %>%
        as_tibble() %>%
        left_join(t_gss, by = 'ctag') %>%
        left_join(t_gs, by = c('ctag','reg.gid','tgt.gid')) %>%
        left_join(tn, by = c('reg.gid','tgt.gid')) %>%
        replace_na(list(binding = 0, weight = 0)) %>%
        mutate(categ = factor(binding))
    tz %>% count(ctag, categ)
    # 
    for (ctag in ctags) {
        ty = tz %>% filter(ctag == !!ctag)
        categ = as.integer(as.character(ty$categ))
        score = as.numeric(ty$weight)
        if(max(score) == 0) score[1] = 0.1
        pr = pr.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        roc = roc.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        t_pr1 = pr$curve %>% as_tibble() %>%
            transmute(nid = nid,
                      ctag = ctag,
                      recall = V1,
                      precision = V2,
                      score = V3)
        t_roc1 = roc$curve %>% as_tibble() %>%
            transmute(nid = nid,
                      ctag = ctag,
                      TPR = V1,
                      FPR = V2,
                      score = V3)
        t_pr = rbind(t_pr, t_pr1)
        t_roc = rbind(t_roc, t_roc1)
        aupr = pr$auc.integral
        auroc = roc$auc
        t_aupr1 = tibble(nid = nid, ctag = ctag, aupr = aupr)
        t_auroc1 = tibble(nid = nid, ctag = ctag, auroc = auroc)
        t_aupr = rbind(t_aupr, t_aupr1)
        t_auroc = rbind(t_auroc, t_auroc1)
        cat(sprintf("%s %10s aupr [%.04f] auroc [%.04f]\n", nid, ctag, aupr, auroc))
    }
}

fo = file.path(dirw, '13_auroc/03.rda')
save(t_pr, t_roc, t_aupr, t_auroc, file = fo)
#}}}

#{{{ plot AUROC
fi = file.path(dirw, '13_auroc/03.rda')
x = load(fi)


#{{{ selected roc/pr plot
nids = c("n41","n42","n71","n81","n82")
cols.dev = c(pal_npg()(4)[2:4], brewer.pal(6,"Paired")[5:6])
tps = t_cfg %>% mutate(lab = sprintf("%s %s", study, tag)) %>%
    filter(nid %in% nids)
labs = tps %>% mutate(nid = factor(nid, levels = nids)) %>%
    arrange(nid) %>% pull(lab)

tp1 = t_roc %>% filter(nid %in% nids) %>%
    inner_join(tps, by = 'nid') %>%
    mutate(ctag = sprintf("%s AUROC", ctag),
           lab = factor(lab, levels = labs))
p1 = ggplot(tp1) +
    geom_line(mapping = aes(x = TPR, y = FPR, color = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'FPR: FP/(FP+TN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'TPR: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = 'none') +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))

tp2 = t_pr %>% filter(nid %in% nids) %>%
    inner_join(tps, by = 'nid') %>%
    mutate(lab = factor(lab, levels = labs),
           ctag = sprintf("%s AUPR", ctag))
p2 = ggplot(tp2) +
    geom_line(mapping = aes(x = recall, y = precission, color = lab)) +
    #geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'Recall: TP/(TP+FN)', breaks=c(.25,.5,.75), limits=c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = 'Precision: TP/(TP+FP)', breaks=c(.25,.5,.75), limits=c(0,1), expand=c(0,0)) +
    scale_color_manual(values = cols.dev) +
    facet_wrap(~ctag, nrow = 1, strip.position = 'top') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = c(.85,.25), legend.justification = c(0,0), legend.background = element_blank()) +
    guides(direction = 'vertical', color = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))

fo = sprintf("%s/13_auroc/10.dev.atlas.pdf", dirw, ctag)
ggarrange(p1, p2, nrow = 2, ncol = 1, heights = c(1,1))  %>% 
    ggexport(filename = fo, width = 10, height = 5)
#}}}

#{{{ aupr/auroc bar-plot
#nids = t_auroc %>% filter(ctag == 'KN1_ear') %>% arrange(auroc) %>% pull(nid)
nids = t_auroc %>% arrange(nid) %>% distinct(nid) %>% pull(nid)
tp = t_aupr %>% left_join(t_auroc, by = c('nid','ctag')) %>%
    mutate(AUPR = aupr, AUROC = auroc) %>% select(-aupr, -auroc) %>%
    gather(type, auc, -nid, -ctag) %>%
    mutate(ctag = sprintf("%s %s", ctag, type)) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lab = str_remove(sprintf("%.03f", auc), '^0+')) %>%
    mutate(nid = factor(nid, levels = rev(nids)))
tps = t_cfg %>% mutate(lab = sprintf("%s_%s", study, tag)) 

p1 = ggplot(tp) +
    #geom_hline(yintercept = .5, size = .3, alpha = .5) +
    geom_bar(mapping = aes(x = nid, y = auc, fill = type), stat = 'identity', width = .75, alpha = .7) +
    geom_text(mapping = aes(x = nid, y = auc , label = lab), hjust = 1, size = 2) +
    scale_x_discrete(breaks = tps$nid, labels = tps$lab) +
    scale_y_continuous(expand = expand_scale(mult=c(0,.05))) +
    scale_fill_npg() +
    coord_flip() +
    facet_grid(.~ctag, scale = 'free') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = 'none') +
    theme(axis.ticks = element_blank()) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
fp = sprintf("%s/13_auroc/05.auc.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 6)
#}}}
#}}}


#{{{ evluate using Briggs B73 & Mo17 data
#{{{ read Briggs DE data
diri = file.path("~/projects/briggs/data", "49.coop")
fi = file.path(diri, "01.master.RData")
x = load(fi)
fd = file.path(diri, "03.sharing.RData")
x = load(fd)
tissues = unique(tm$Tissue)
t_de = tm %>% filter(silent == 0) %>% 
    transmute(Tissue = Tissue, gid = gid, DE = pDE) 
t_de %>% count(Tissue, DE)
#}}}

net_sizes = c(10000, 100000, 1000000)
max_net_size = max(net_sizes)
t_bv = tibble()
for (i in 1:nrow(th)) {
#for (i in 1:3) {
    nid = th$nid[i]
    fgn = th$fgn[i]
    cat(nid, "\n")
    x = load(fgn)
    x
    for (tissue in tissues) {
        t_de1 = t_de %>% filter(Tissue == tissue, gid %in% tids) %>% select(-Tissue)
        gids = t_de1$gid
        rids1 = rids[rids %in% gids]
        tids1 = tids[tids %in% gids]
        tn1 = tn %>% filter(reg.gid %in% rids1, tgt.gid %in% tids1) %>%
            filter(row_number() <= max_net_size) %>%
            inner_join(t_de1, by = c('reg.gid' = 'gid')) %>%
            rename(reg.DE = DE) %>%
            inner_join(t_de1, by = c('tgt.gid' = 'gid')) %>%
            rename(tgt.DE = DE)
        for (net_size in net_sizes) {
            tn2 = tn1 %>% filter(row_number() <= net_size) %>%
                count(reg.DE, tgt.DE) %>%
                mutate(nid = nid, tissue = tissue, net_size = net_size) %>%
                select(nid, tissue, net_size, everything())
            t_bv = rbind(t_bv, tn2)
        }
    }
}

net_size = 100000
for (net_size in net_sizes) {
tp = t_bv %>% filter(net_size == !!net_size) %>%
    mutate(reg.DE = ifelse(reg.DE == 'non_DE', 'non_DE', 'DE'),
           tgt.DE = ifelse(tgt.DE == 'non_DE', 'non_DE', 'DE')) %>%
    group_by(nid, tissue, net_size, reg.DE, tgt.DE) %>%
    summarise(n = sum(n)) %>% ungroup() %>%
    spread(tgt.DE, n) %>%
    mutate(prop.de = DE/(DE+non_DE)) %>% select(-DE,-non_DE) %>%
    mutate(nid = factor(nid, levels = rev(t_cfg$nid)),
           lab = sprintf("%.02f", prop.de))
    #spread(reg.DE, prop.de) %>%
    #mutate(fc = DE/non_DE)
tps = t_cfg %>% mutate(lab = sprintf("%s_%s", study, tag)) 
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = nid, y = prop.de, fill = reg.DE), stat = 'identity', position = 'dodge', width = .9, alpha = .7) +
    geom_text(mapping = aes(x = nid, y = prop.de + .01, group = reg.DE, label = lab), hjust = 0, position = position_dodge(width = 1), size = 2) +
    scale_x_discrete(breaks = tps$nid, labels = tps$lab, expand = c(.01,0)) +
    scale_y_continuous(name = 'Prop. DE Targets', limits = c(0, 1), expand = c(0,0)) +
    scale_fill_npg() +
    coord_flip() +
    facet_wrap(~tissue, nrow = 3) +
    theme_bw() +
    #theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.6), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks = element_blank()) +
    theme(plot.margin = unit(c(1.5,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
fp = sprintf("%s/15_bv_%dk.pdf", dirw, net_size/1000)
ggsave(p1, filename = fp, width = 10, height = 15)
}
#}}}

#{{{ compare GRNs
t_grn = list()
#for (i in 1:nrow(th)) {
for (i in c(2,3,14)) {
    nid = th$nid[i]; fgn = th$fgn[i]
    x = load(fgn)
    if(! nid %in% t_grn) 
        t_grn[[nid]] = list(
            rids = rids,
            tids = tids,
            tn = tn
        )
}
names(t_grn)

nid1 = 'n53'
nid2 = 'n55'
regs1 = t_grn[[nid1]][['rids']] 
regs2 = t_grn[[nid2]][['rids']] 
tgts1 = t_grn[[nid1]][['tids']] 
tgts2 = t_grn[[nid2]][['tids']] 
regs = regs1[regs1 %in% regs2]
tgts = tgts1[tgts1 %in% tgts2]
nreg = length(regs)
ntgt = length(tgts)
cat(sprintf("%d regs x %d tgts: %d total edges\n", nreg, ntgt, nreg*ntgt))

nsize = 1000000
tn1 = t_grn[[nid1]][['tn']] %>% filter(reg.gid %in% regs, tgt.gid %in% tgts) %>%
    filter(row_number() <= nsize)
tn2 = t_grn[[nid2]][['tn']] %>% filter(reg.gid %in% regs, tgt.gid %in% tgts) %>%
    filter(row_number() <= nsize)
tno = tn1 %>% inner_join(tn2, by = c('reg.gid','tgt.gid'))
novlp = nrow(tno)
novlp

#}}}


