source("functions.R")
dirw = file.path(dird, '92_te_gene')

#{{{ read, validate and transform
diri = '/home/springer/sna/te_expression/method_paper_combo'
fg1 = file.path(diri, 'B73_gene_RPM_Zhou.txt')
fg2 = file.path(diri, 'Mo17_gene_RPM_Zhou.txt')
ft1 = file.path(diri, 'B73_TEfamily_RPM_Zhou.txt')
ft2 = file.path(diri, 'Mo17_TEfamily_RPM_Zhou.txt')

tg1h = read_tsv(fg1, col_names = F, n_max = 1)
tg1 = read_tsv(fg1, col_names = c('gid', as.character(tg1h)), skip = 1)
tg2h = read_tsv(fg2, col_names = F, n_max = 1)
tg2 = read_tsv(fg2, col_names = c('gid', as.character(tg2h)), skip = 1)
tt1h = read_tsv(ft1, col_names = F, n_max = 1)
tt1 = read_tsv(ft1, col_names = c('gid', as.character(tt1h)), skip = 1)
tt2h = read_tsv(ft2, col_names = F, n_max = 1)
tt2 = read_tsv(ft2, col_names = c('gid', as.character(tt2h)), skip = 1)

identical(tg1$gid, tg2$gid)
identical(tt1$gid, tt2$gid)
gids = tg1$gid
tids = tt1$gid

colnames(tt1)[32] = colnames(tg1)[32]
identical(colnames(tg1), colnames(tt1))
identical(colnames(tg2), colnames(tt2))

tb = tg1 %>% bind_rows(tt1) %>% gather(tissue, rpm, -gid)
tm = tg2 %>% bind_rows(tt2) %>% gather(tissue, rpm, -gid)

gids_b = tb %>% group_by(gid) %>%
    summarise(psam = sum(rpm>1)/length(rpm), rpm_sd = sd(rpm)) %>%
    ungroup() %>% filter(psam > .1, rpm_sd > 0) %>% pull(gid)
gids_m = tm %>% group_by(gid) %>%
    summarise(psam = sum(rpm>1)/length(rpm), rpm_sd = sd(rpm)) %>%
    ungroup() %>% filter(psam > .1, rpm_sd > 0) %>% pull(gid)
tb = tb %>% filter(gid %in% gids_b)
tm = tb %>% filter(gid %in% gids_m)
#}}}

# write output for seidr
gt = 'B73'
gt = 'Mo17'
if(gt == 'B73') {
    to = tb
} else {
    to = tm
}
to = to %>% spread(gid, rpm)

aids = colnames(to)[-1]
fo = sprintf("%s/01.%s.tsv", dirw, gt)
write_tsv(to[-1], fo, col_names = F)
fo = sprintf("%s/01.%s.gid.tsv", dirw, gt)
write(aids, file=fo)
tids = tids[tids %in% aids]
fo = sprintf("%s/01.%s.tid.tsv", dirw, gt)
write(tids, file=fo)


