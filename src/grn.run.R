source("functions.R")
require(future)
require(furrr)
dirw = file.path(dird, '14_eval_sum')
gopt = 'rf'

tgo = gs$fun_ann
#tgos = tgo %>% count(ctag, grp) %>% filter(n >= 3) %>% select(-n)
#tgo = tgo %>% inner_join(tgos, by = c("ctag", 'grp')) %>%
    #group_by(ctag) %>% nest() %>% rename(tgrp=data)

prep_net <- function(tn, bins=10, max_score=4) {
    #{{{
    tn %>%
        mutate(pcc_sign=ifelse(pcc < 0, '-', '+')) %>%
        mutate(score=as.integer(cut_interval(score,bins))) %>%
        mutate(score = ifelse(score >= max_score, max_score, score)) %>%
        select(group=reg.gid,gid=tgt.gid,score,pcc_sign)
    #}}}
}
filter_tn_by_score <- function(tn, score=1)
    tn %>% filter(score >= !!score) %>% select(-score)
fi = sprintf("%s/%s.50k.rds", dirr, gopt)
tn = readRDS(fi) %>% mutate(tn2 = map(tn, prep_net, bins=10, max_score=4)) %>%
    select(nid, tids, tn=tn2) %>%
    mutate(s1 = map(tn, filter_tn_by_score, score = 1)) %>%
    mutate(s2 = map(tn, filter_tn_by_score, score = 2)) %>%
    mutate(s3 = map(tn, filter_tn_by_score, score = 3)) %>%
    mutate(s4 = map(tn, filter_tn_by_score, score = 4)) %>%
    select(nid, tids, s1, s2, s3, s4) %>%
    gather(score, tn, -nid, -tids) %>%
    mutate(score = as.integer(str_replace(score, '^s', '')))

plan(multicore, workers=4)
options(future.globals.maxSize=4.29e9)

hyper_enrich_m <- function(tg, gids, tgrp) {
    #{{{ tg: group + gid;  tgrp: grp, gid, note
    popSize = length(gids)
    tg1 = tg %>% count(group) %>% filter(n>=3) %>% rename(sampleSize = n)
    n_group = nrow(tg1)
    tgrp = tgrp %>% filter(gid %in% gids)
    tgrp1 = tgrp %>% count(ctag, grp) %>% filter(n>3) %>% rename(hitInPop = n)
    tgrp2 = tgrp1 %>% count(ctag) %>% rename(n_grp = n)
    tg %>% inner_join(tgrp, by='gid') %>%
        count(group, ctag, grp) %>% rename(hitInSample=n) %>%
        inner_join(tg1, by='group') %>%
        inner_join(tgrp1, by=c('ctag','grp')) %>%
        mutate(n_group = n_group, popSize = popSize) %>%
        inner_join(tgrp2, by='ctag') %>%
        mutate(pval.raw = phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail=F)) %>%
        #mutate(ratioInSample = sprintf("%d/%d", hitInSample, sampleSize)) %>%
        #mutate(ratioInPop = sprintf("%d/%d", hitInPop, popSize)) %>%
        select(n_group,group,ctag,n_grp,grp,hitInSample,sampleSize,hitInPop,popSize,pval.raw)
    #}}}
}

tm = tn %>% #slice(1:5) %>%
    mutate(res = future_map2(tn, tids, hyper_enrich_m, tgrp=tgo, .progress=T)) %>%
    select(nid,score,res) %>% unnest() %>%
    group_by(nid,score,n_group,ctag,n_grp) %>%
    mutate(pval.adj = p.adjust(pval.raw, method='BH')) %>%
    ungroup()

tms = tm %>% filter(pval.adj < .05) %>%
    distinct(nid,score,n_group, group, ctag) %>%
    count(nid, score, n_group, ctag) %>%
    mutate(p_group_enrich=n/n_group)

fo = sprintf('%s/raw/%s.go2.rds', dird, gopt)
saveRDS(tm, file=fo)

if(FALSE) {
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

res=list(filt=ev_br_filt,sum=ev_br_sum)
fo = file.path(dirw, '01.br.filt.rds')
saveRDS(res, file = fo)
#}}}
}
