source("functions.R")
dirg = '/home/springer/zhoux379/data/genome/B73'
dirp = '/home/springer/zhoux379/projects/maize.grn'
dird = file.path(dirp, 'data')
dirw = file.path(dird, '08_y1h_45')
dirw

# read v3 to v4 mapping table
fm = file.path(dirg, "gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read_tsv(fm, col_names = F) %>%
    transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5) %>%
    select(ogid, gid, type)

#{{{ Top45 TFs and targets by Y1H
fi = file.path(dirw, "y1h.targets.tsv")
ti = read_tsv(fi)
colnames(ti) = c("reg.v3", "reg.v4", "tgt.v3", "tgt.v4")
ti = ti %>% fill(reg.v3, reg.v4, .direction = 'down')
ti %>% distinct(reg.v3)
ti %>% distinct(tgt.v3)

tch = ti %>% distinct(reg.v3, reg.v4) %>% 
    left_join(tmr, by = c('reg.v3' = 'ogid')) %>%
    print(n = 45)
tch %>% filter(is.na(gid) | reg.v4 != gid)

tr = ti %>% filter(reg.v4 != 'none', !is.na(tgt.v4)) %>%
    transmute(reg = reg.v4, tgt = tgt.v4)
fr = file.path(dirw, '10.tsv')
write_tsv(tr, fr)
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
