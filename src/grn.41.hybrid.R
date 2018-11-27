source(file.path(dirr, "grn.fun.r"))
require(PRROC)
#source("enrich.R")
diri = '~/projects/maize.expression'
dirw = '~/projects/maize.grn/data'
f_cfg = file.path(dirw, '10.dataset.tsv')
t_cfg = read_tsv(f_cfg)
t_cfg = t_cfg %>% mutate(note = ifelse(is.na(note), study, sprintf("%s_%s", study, note)))
studies = t_cfg %>% distinct(study) %>% pull(study)

#{{{ evaluate B73 and Mo17 network, hybrid, ovlp w. others
dirw = file.path(dird, '25_briggs')
fi = file.path(dird, "12_output/n99b_1.rda")
load(fi)
rids.b = rids; tids.b = tids; reg.mat.b = reg.mat
fi = file.path(dird, "12_output/n99b_2.rda")
load(fi)
rids.m = rids; tids.m = tids; reg.mat.m = reg.mat
fi = file.path(dird, "12_output/n99b_3.rda")
load(fi)
rids.h = rids; tids.h = tids; reg.mat.h = reg.mat

rids = rids.b[rids.b %in% rids.m]; rids = rids[rids %in% rids.h]
tids = tids.b[tids.b %in% tids.m]; tids = tids[tids %in% tids.h]
reg.mat.b = reg.mat.b[rids, tids]
reg.mat.m = reg.mat.m[rids, tids]
reg.mat.h = reg.mat.h[rids, tids]

max_net_size = 1000000
tnb = as_tibble(reg.mat.b) %>% mutate(rid = rids) %>%
    gather(tid, score, -rid) %>% filter(rid != tid) %>%
    arrange(desc(score)) %>% filter(row_number() <= max_net_size)
tnm = as_tibble(reg.mat.m) %>% mutate(rid = rids) %>%
    gather(tid, score, -rid) %>% filter(rid != tid) %>%
    arrange(desc(score)) %>% filter(row_number() <= max_net_size)
tnh = as_tibble(reg.mat.h) %>% mutate(rid = rids) %>%
    gather(tid, score, -rid) %>% filter(rid != tid) %>%
    arrange(desc(score)) %>% filter(row_number() <= max_net_size)

net_size = 1000000
fo = sprintf("%s/vennd.%d.pdf", dirw, net_size)
grn_venn(net_size, tnb, tnm, tnh, fo)

grn_venn <- function(net_size = 10000, tnb, tnm, tnh, fo) {
    require(eulerr)
    edges.b = tnb[1:net_size,] %>% mutate(ename = sprintf("%s_%s", rid, tid)) %>% pull(ename)
    edges.m = tnm[1:net_size,] %>% mutate(ename = sprintf("%s_%s", rid, tid)) %>% pull(ename)
    edges.h = tnh[1:net_size,] %>% mutate(ename = sprintf("%s_%s", rid, tid)) %>% pull(ename)
    teb = tibble(pair = edges.b, eb = 1) 
    tem = tibble(pair = edges.m, em = 1)
    teh = tibble(pair = edges.h, eh = 1)
    tp = teb %>% full_join(tem, by = 'pair') %>% full_join(teh, by = 'pair')
    input = list("B73" = edges.b, "Mo17" = edges.m, "BxM" = edges.h)
    mat = tp %>% select(-pair) %>% 
        mutate(eb = ifelse(is.na(eb), F, T)) %>%
        mutate(em = ifelse(is.na(em), F, T)) %>%
        mutate(eh = ifelse(is.na(eh), F, T)) %>%
        transmute(B73 = eb, Mo17 = em, BxM = eh)
    #
    #eulerr_options(pointsize = 8)
    #options(digits = 4)
    fit <- euler(mat)
    pdf(fo, width = 4, height = 4)
    print(plot(fit,
         quantities = T,
         fills =  T, #"transparent",
         lty = 1:3,
         #labels = c(list(font = 8))
         ))
    dev.off()
}

#{{{ upset visualization
require(UpSetR)
fo = file.path(dirw, 'venn.diagram.pdf') 
pdf(fo, width = 8, height = 8)
grid.newpage()
upset(fromList(input), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.off()
#}}}
#}}}



