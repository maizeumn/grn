#{{{ head
source("grn.fun.r")
require(pROC)
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
t_roc = tibble()
t_auroc = tibble()
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
        roc_obj <- roc(ty$categ, ty$weight)
        sensi = roc_obj$sensitivities
        speci = roc_obj$specificities
        #sensi = c(sensi, rep(0, nrow(ty)-length(sensi)))
        #speci = c(speci, rep(1, nrow(ty)-length(speci)))
        t_roc1 = tibble(nid = nid,
                        ctag = ctag,
                        TPR = rev(sensi), 
                        FPR = rev(1 - speci))
            #lab = roc_obj$response, 
            #score = roc_obj$predictor)
        t_roc = rbind(t_roc, t_roc1)
        z = auc(roc_obj)
        t_auroc1 = tibble(nid = nid, ctag = ctag, auroc = z)
        t_auroc = rbind(t_auroc, t_auroc1)
        cat(sprintf("%s %s auroc: %f\n", nid, ctag, z))
    }
}

fo = file.path(dirw, '13_auroc/03.rda')
save(t_roc, t_auroc, file = fo)
#}}}

#{{{ plot AUROC
fi = file.path(dirw, '13_auroc/03.rda')
x = load(fi)

ctag = ctags[1]
#for (ctag in ctags) {
#tp = t_roc %>% filter(ctag == !!ctag) %>% select(-ctag)
tps = t_auroc %>% filter(ctag == !!ctag) %>% select(-ctag) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lab = sprintf("[AUC=%.03f] %s_%s", auroc, study, tag)) %>%
    #mutate(col = rep(pal_npg()(6), each = 4)[1:nrow(th)],
    #       lty = rep(c('solid','dotted','dashed','dotdash'), 6)[1:nrow(th)])
    mutate(lcol = pal_d3()(5))

tp = t_roc %>% filter(nid %in% c("n41","n42","n71","n81","n82")) %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lab = sprintf("%s %s", study, tag))
p1 = ggplot(tp) +
    geom_line(mapping = aes(x = FPR, y = TPR, color = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    scale_x_continuous(name = 'FPR: FP/(FP+TN)', breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_y_continuous(name = 'TPR: TP/(TP+FN)', breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_color_d3() +
    #scale_color_manual(name = toupper(ctag), values = tps$col, labels = tps$lab) +
    #scale_linetype_manual(name = toupper(ctag), values = tps$lty, labels = tps$lab) +
    facet_wrap(~ctag, ncol = 3, strip.position = 'top') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = c(.2,.1), legend.justification = c(0,0), legend.background = element_blank()) +
    guides(direction = 'vertical', color = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
fp = sprintf("%s/13_auroc/10.dev.atlas.pdf", dirw, ctag)
ggsave(p1, filename = fp, width = 8, height = 6)

nids = t_auroc %>% filter(ctag == 'KN1_ear') %>% arrange(auroc) %>% pull(nid)
tp = t_auroc %>%
    inner_join(t_cfg, by = 'nid') %>%
    mutate(lab = sprintf("%.03f", auroc)) %>%
    mutate(nid = factor(nid, levels = nids))
tps = t_cfg %>% mutate(lab = sprintf("%s_%s", study, tag)) 
p1 = ggplot(tp) +
    geom_hline(yintercept = .5, size = .3, alpha = .5) +
    geom_bar(mapping = aes(x = nid, y = auroc, fill = ctag), stat = 'identity', width = .75, alpha = .7) +
    geom_text(mapping = aes(x = nid, y = auroc - .01 , label = lab), hjust = 1, size = 2) +
    scale_x_discrete(breaks = tps$nid, labels = tps$lab) +
    scale_y_continuous(name = 'AUROC', breaks = c(.25,.5,.75), limits = c(0, 1), expand = c(0,0)) +
    scale_fill_npg() +
    coord_flip() +
    facet_grid(.~ctag) +
    theme_bw() +
    #theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = 'none') +
    theme(axis.ticks.y = element_blank()) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
fp = sprintf("%s/13_auroc/05.auroc.pdf", dirw, ctag)
ggsave(p1, filename = fp, width = 8, height = 6)
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
tp = t_bv %>% filter(net_size == net_size) %>%
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
ggsave(p1, filename = fp, width = 12, height = 8)
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


