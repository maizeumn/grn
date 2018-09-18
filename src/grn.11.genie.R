#{{{ head
source("grn.fun.r")
require(PRROC)
#source("enrich.R")
diri = '~/projects/maize.expression'
dirw = '~/projects/maize.grn/data'
f_cfg = file.path(dirw, '10.genie3.tsv')
t_cfg = read_tsv(f_cfg)
th = t_cfg %>% mutate(fgn = sprintf("%s/12_output/%s.rda", dirw, nid)) 
studies = t_cfg %>% distinct(study) %>% pull(study)
#
eval1 = 'lin2017'
nids1 = th %>% filter(study == !!eval1) %>% pull(nid)
eval2 = 'kremling2018'
nids2 = th %>% filter(study == !!eval2) %>% pull(nid)
eval3 = 'all'
nids3 = th %>% filter(!study %in% c("lin2017", "kremling2018")) %>% pull(nid)
lst_eval = list()
lst_eval[[eval1]] = nids1; lst_eval[[eval2]] = nids2; lst_eval[[eval3]] = nids3;

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

#{{{ read known TF targets / GGIs
fi = file.path(dirw, '09.gs.rda')
x = load(fi)
x
t_gs = t_gs %>% 
    filter(! ctag %in% c("KN1_any","KN1_ear","KN1_tassel","KN1_leaf")) %>%
    mutate(ctag = ifelse(ctag=='KN1_all', 'KN1', ctag)) %>% 
    mutate(binding = 1)
t_gs = t_gs %>% filter(! ctag %in% c("GO", "CornCyc", "PPIM"))
t_gs %>% count(ctag)
ctags = unique(t_gs$ctag)
reg.gids = unique(t_gs$reg.gid)
t_gss = t_gs %>% distinct(ctag, reg.gid)
#}}}

#{{{ prepare genie3 input
#study = 'li2013'
#study = 'hirsch2014'
#study = 'leiboff2015'
#study = 'jin2016'
#study = 'stelpflug2016'
study = 'walley2016'
#study = 'lin2017'
#study = 'kremling2018'
#study = 'briggs'
#study = 'dev41'
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

#th = th %>%
#    mutate(done = file.exists(fgn)) %>%
#    filter(done) %>%
#    select(-done)
#th = th %>% filter(nid %in% c('n41', 'n42', 'n71', 'n81', 'n82'))

#{{{ manually merge tissue GRNs
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
#}}}

#{{{ use PCC to predict TF targets
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
#}}}


#{{{ compute AUROC using known TF data
t_pr = tibble(); t_roc = tibble()
t_aupr = tibble(); t_auroc = tibble()
for (i in 1:nrow(th)) {
    nid = th$nid[i]; fgn = th$fgn[i]
    #nid = 'n71'
    #fgn = th$fgn[th$nid==nid]
    x = load(fgn)
    x
    if('score' %in% colnames(tn)) tn = tn %>% rename(weight = score)
    # 
    ctags_v = t_gss %>% filter(reg.gid %in% rids) %>% pull(ctag)
    tz = expand.grid(ctag = ctags_v, tgt.gid = tids, stringsAsFactors = F) %>%
        as_tibble() %>%
        left_join(t_gss, by = 'ctag') %>%
        left_join(t_gs, by = c('ctag','reg.gid','tgt.gid')) %>%
        left_join(tn, by = c('reg.gid','tgt.gid')) %>%
        replace_na(list(binding = 0, weight = 0)) %>%
        mutate(categ = factor(binding))
    tz %>% count(ctag, categ)
    # 
    for (ctag in ctags_v) {
        ty = tz %>% filter(ctag == !!ctag)
        categ = as.integer(as.character(ty$categ))
        score = as.numeric(ty$weight)
        if(max(score) == 0) score[1] = 0.1
        pr = pr.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        roc = roc.curve(scores.class0 = score, weights.class0 = categ, curve = T)
        t_pr1 = pr$curve %>% as_tibble() %>%
            transmute(nid = nid, ctag = ctag,
                      recall = V1, precision = V2, score = V3)
        t_roc1 = roc$curve %>% as_tibble() %>%
            transmute(nid = nid, ctag = ctag,
                      TPR = V1, FPR = V2, score = V3)
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
#
fo = file.path(dirw, '13_auroc/03.rda')
save(t_pr, t_roc, t_aupr, t_auroc, file = fo)
#}}}

#{{{ plot AUROC
dirw = file.path(dird, '13_auroc')
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
    geom_line(mapping = aes(x = recall, y = precision, color = lab)) +
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
tpa = t_aupr %>% filter(ctag %in% ctags) %>%
    left_join(t_auroc, by = c('nid','ctag')) %>%
    inner_join(t_cfg, by = 'nid') %>%
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
#nids = t_auroc %>% filter(ctag == 'KN1_ear') %>% arrange(auroc) %>% pull(nid)
#nids = t_auroc %>% arrange(nid) %>% distinct(nid) %>% pull(nid)
#}}}
#}}}


#{{{ evluate using Briggs B73 & Mo17 data
dirw = file.path(dird, '15_de_val')
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
t_des = t_de %>% group_by(Tissue) %>%
    summarise(propDE = sum(DE!='non_DE') / n())
#}}}

#{{{ run the pipeline
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
                group_by(reg.DE, tgt.DE) %>%
                summarise(n = n(), nreg = length(unique(reg.gid)),
                          ntgt = length(unique(tgt.gid))) %>%
                ungroup() %>%
                mutate(nid = nid, tissue = tissue, net_size = net_size) %>%
                select(nid, tissue, net_size, everything())
            t_bv = rbind(t_bv, tn2)
        }
    }
}
fo = file.path(dirw, "01.rda")
save(t_bv, file = fo)
#}}}

#{{{ plot
fi = file.path(dirw, "01.rda")
x = load(fi)
tissues6 = c('auricle_v12', 'ear_v14', 'embryo_27DAP', 'kernel_14DAP', 'root_0DAP',
            'seedlingmeristem_11DAS')

opt = 'kremling2018'
nids = lst_eval[[opt]]
net_sizes = c("10k", "100k", "1m")
net_size_map = c("10000"="10k", "1e+05"="100k", "1e+06"="1m")
tp = t_bv %>% #filter(net_size == !!net_size) %>%
    filter(tissue %in% tissues6, nid %in% nids) %>%
    mutate(reg.DE = ifelse(reg.DE == 'non_DE', 'non_DE', 'DE'),
           tgt.DE = ifelse(tgt.DE == 'non_DE', 'non_DE', 'DE')) %>%
    group_by(nid, tissue, net_size, reg.DE, tgt.DE) %>%
    summarise(n = sum(n), nreg = sum(nreg), ntgt = sum(ntgt)) %>% ungroup() %>%
    filter(reg.DE == 'DE') %>%
    group_by(nid, tissue, net_size) %>%
    summarise(prop.de = ntgt[tgt.DE == 'DE']/sum(ntgt), 
              nreg = sum(nreg), ntgt = sum(ntgt)) %>% ungroup() %>%
    mutate(txt = sprintf("%d", nreg)) %>% select(-nreg, -ntgt) %>%
    mutate(lab = str_remove(sprintf("%.02f", prop.de), '^0+')) %>%
    left_join(t_cfg, by = 'nid') %>%
    mutate(tag = sprintf("%s_%s", study, tag)) %>% select(-study) %>%
    mutate(nid = factor(nid, levels = rev(nids))) %>%
    mutate(net_size = factor(net_size_map[as.character(net_size)], levels = net_sizes)) %>%
    mutate(tissue = factor(tissue, levels = tissues6))
    #spread(reg.DE, prop.de) %>%
    #mutate(fc = DE/non_DE)
tps = tp %>% distinct(nid, tag)
tpl = t_des %>% transmute(tissue = Tissue, prop.de = propDE) %>%
    filter(tissue %in% tissues6)
#
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = nid, y = prop.de, fill = net_size), stat = 'identity', position = 'dodge', width = .9, alpha = .7) +
    geom_text(mapping = aes(x = nid, y = prop.de + .01, group = net_size, label = lab), hjust = 0, position = position_dodge(width = 1), size = 2) +
    geom_text(mapping = aes(x = nid, y = .005, group = net_size, label = txt), hjust = 0, position = position_dodge(width = 1), size = 2) +
    geom_hline(data = tpl, aes(yintercept = prop.de), size = .3, alpha = .5) +
    scale_x_discrete(breaks = tps$nid, labels = tps$tag, expand = c(.01,0)) +
    scale_y_continuous(name = 'Prop. DE Targets', expand = expand_scale(mult=c(0,.1))) +
    scale_fill_d3() +
    coord_flip() +
    facet_wrap(~tissue, nrow = 1, scale = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9, margin = margin(0,0,.2,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks = element_blank()) +
    theme(plot.margin = unit(c(1.5,.2,.2,.2), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
fp = sprintf("%s/03.%s.pdf", dirw, opt)
ggsave(p1, filename = fp, width = 8, height = 6)
#}}}
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


