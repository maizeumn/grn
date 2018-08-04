source("grn.fun.r")

#{{{ read in
fi = file.path(dird, '05.previous.grns/10.RData')
x = load(fi)
fi = file.path(dird, '07.known.tf/10.RData')
x = load(fi)

tta = t_grn %>% group_by(ctag) %>%
    summarise(n.edge = n(),
              n.reg = length(unique(regulator)),
              n.tgt = length(unique(target)),
              score.min = min(score),
              score.max = max(score),
              score.median = median(score))
tt = tta %>%
    transmute('GRN source' = ctag,
              N_Edge = n.edge, N_Regulator = n.reg, N_Target = n.tgt,
              min_score = sprintf("%.02e", score.min),
              max_score = sprintf("%.02e", score.max),
              median_score = sprintf("%.02e", score.median))
fo = file.path(dira, '02.net.stats.RData')
save(tt, file = fo)
#}}}

#{{{ check ovlp of y1h w. Huang et al and Walley et al GRN
fi = file.path(dird, '08.y1h.45/10.tsv')
ti = read_tsv(fi)
dirw = file.path(dira, "04_y1h")
ti2 = ti %>% mutate(ctag = 'y1h', score = NA)
tx = ti %>% inner_join(t_grn[,1:4], by = c("reg"="regulator",'tgt'='target'))
# color edges according to source GRN
ctags = c("y1h", unique(tx$ctag))
cols.e = c("gray70", pal_npg()(length(ctags)-1))
names(cols.e) = ctags
tp = rbind(ti2, tx) %>% mutate(e.col = cols.e[ctag])

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
fo = file.path(dirw, '01.pdf')
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

#{{{ output to table
ti %>% distinct(reg)
ti %>% distinct(tgt)
nrow(ti)
ty = ti %>% left_join(t_grn[,1:4], by = c("reg"="regulator",'tgt'='target'))
to = ty %>% spread(ctag, score) %>% select(-`<NA>`)
sum(!is.na(ty$ctag))
fo = file.path(dirw, '01.tsv')
write_tsv(to, fo)
#}}}
#}}}

