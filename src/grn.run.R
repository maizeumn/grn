source("functions.R")
dirw = file.path(dird, '14_eval_sum')
diri = '~/projects/rnaseq'


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
