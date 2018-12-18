source("functions.R")
dirw = file.path(dird, '08_y1h')
tm = v3_to_v4()

#{{{ save Top45(44) TFs and 123 targets to 01.rds
fi = file.path(dirw, "y1h.targets.tsv")
ti = read_tsv(fi)
colnames(ti) = c("reg.v3", "reg.v4", "tgt.v3", "tgt.v4")
ti = ti %>% fill(reg.v3, reg.v4, .direction = 'down')
ti %>% distinct(reg.v3)
ti %>% distinct(tgt.v3)

tch = ti %>% distinct(reg.v3, reg.v4) %>%
    left_join(tm, by = c('reg.v3' = 'ogid')) %>%
    print(n = 45)
tch %>% filter(is.na(gid) | reg.v4 != gid)

ti2 = ti %>% filter(reg.v4 != 'none', !is.na(tgt.v4))
t_y1h = ti2 %>%
    transmute(reg=reg.v4, tgt=tgt.v4, reg.v3=reg.v3, tgt.v3=tgt.v3)
t_reg = ti2 %>% distinct(reg.v3, reg.v4) %>%
    transmute(gid=reg.v4, gid_v3 = reg.v3)

fi = file.path(dirw, 'phenolic.xlsx')
ti = read_xlsx(fi,
    col_names = c('gid_v3','gname','yeast_prom','prom_mplant',
                  'gid','chrom','start','end','orig','note',
                  'ref','x1','coexp','syntelog','sub1','sub2_ref')) %>%
    filter(row_number() > 1) %>%
    filter(!is.na(gid) & gid != 'Na') %>%
    select(gid,gname,everything()) %>%
    filter(gid_v3 != 'GRMZM2G049424')
ti %>% count(gid) %>% count(n)
t_tgt = ti

fo = file.path(dirw, '01.rds')
res = list(reg.gids=t_reg$gid, tgt.gids=t_tgt$gid,
           t_reg=t_reg, t_tgt=t_tgt, t_y1h=t_y1h)
saveRDS(res, file=fo)
#}}}

#{{{ plot top45 TF GRN
fi = file.path(dirw, '10.tsv')
ti = read_tsv(fi)
ti2 = ti %>% mutate(ctag = 'y1h', score = NA)
ctags = c("y1h")
cols.e = c("gray70", pal_npg()(length(ctags)-1))
names(cols.e) = ctags
tp = rbind(ti2) %>% mutate(e.col = cols.e[ctag])

#{{{ plot GRN
cols.v = c("tomato", "gray50")
net = graph_from_data_frame(tp, directed = T)
#net <- simplify(net, remove.multiple = F, remove.loops = T)
deg <- degree(net, mode="all")
V(net)$size <- deg
vs = names(V(net))
V(net)$color = cols.v[2]
V(net)$color[vs %in% tp$reg] = cols.v[1]
vs.trimmed = str_sub(vs, 7, -1)
lab.locs = radian.rescale(x = 1:length(vs), direction = -1, start = 0)
la = layout_with_fr(net)
la = layout_in_circle(net)
fo = file.path(dirw, '91.pdf')
pdf(fo, width = 8, height = 8)
plot(net, edge.arrow.size=.5, 
     edge.color = E(net)$e.col,
     edge.curved = curve_multiple(net),
     vertex.label = NA, layout = la)
#
x = la[,1]*1.2
y = la[,2]*1.2
angle = ifelse(atan(-(la[,1]/la[,2]))*(180/pi) < 0,  90 + atan(-(la[,1]/la[,2]))*(180/pi), 270 + atan(-la[,1]/la[,2])*(180/pi))
for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=vs.trimmed[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
}
#
legend(x=-1.5, y=1.3, xjust = 0, c("Regulator/TF", "Target"), pch=21,
    col="#777777", pt.bg = cols.v, pt.cex=2, cex=.8, bty="n", ncol=1)
legend(x=-1.5, y=-1, xjust = 0, names(cols.e),
    col= cols.e, lty = 1, lwd = 2, cex=.8, bty="n", ncol=1)
dev.off()
#}}}
#}}}

#{{{ expression in dev139
fi = file.path(dirw, '01.rds')
res = readRDS(fi)
t_reg = res$t_reg
t_tgt = res$t_tgt %>% select(gid, gname, gid_v3)

study='mec03'
th = rnaseq_sample_meta(study)
tm = rnaseq_sample_cpm(study)

tor = tm %>% filter(gid %in% t_reg$gid) %>%
    select(gid,SampleID,CPM) %>%
    inner_join(th, by='SampleID') %>%
    mutate(Tissue=sprintf("%s|%s|%d", Tissue, Treatment,Replicate)) %>%
    select(gid,Tissue,CPM) %>%
    spread(Tissue,CPM)
#}}}


