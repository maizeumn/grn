#require(plyranges)
source("functions.R")
dirw = file.path(dird, 'uniformmu')
gs = readRDS('~/projects/grn/data/09.gs.rds')

#{{{ fix UniformMu BED
fi = file.path(dirw, "00.uniformmu.gff3")
ti = read_tsv(fi, col_names = F) %>%
    mutate(X1 = str_replace(X1, "Chr([\\d+])", "\\1")) %>%
    mutate(chrom = X1, start = X4, end = X5, note = X9) %>%
    separate(note, c("namestr", "stockstr"), sep = ";") %>%
    separate(namestr, c("namet", "mid"), sep = "=") %>%
    separate(stockstr, c("stockt", "sids"), sep = "=") %>%
    mutate(start = start - 1) %>%
    select(chrom, start, end, mid, sids) %>%
    arrange(chrom, start, end)
tis = ti %>% count(mid)
mids1 = tis %>% filter(n==1) %>% pull(mid)
mids2 = tis %>% filter(n>1) %>% pull(mid)
ti1 = ti %>% filter(mid %in% mids1)
ti2 = ti %>% filter(mid %in% mids2)

tu = ti2 %>%
    mutate(seqnames=chrom) %>% select(-chrom) %>%
    as_granges() %>%
    group_by(mid) %>%
    reduce_ranges(n_mid = length(mid), sids = paste(sids, collapse=','))

str_dedup <- function(strs, sep=',') paste(unique(unlist(str_split(strs, sep))), collapse=sep)
tu2 = tu %>% as_tibble() %>% rename(chrom=seqnames) %>%
    select(-strand,-width,-n_mid) %>%
    mutate(sids = map_chr(sids, str_dedup)) %>%
    arrange(mid, chrom, start) %>%
    #group_by(mid) %>% mutate(idx = 1:n()) %>% ungroup() %>%
    #mutate(mid = ifelse(idx > 1, sprintf("%s_%d",mid,idx), mid)) %>%
    select(chrom,start,end,mid,sids)
nrow(ti2) - nrow(tu2)

to = rbind(ti1, tu2)
nrow(ti) - nrow(to)
fo = file.path(dirw, "01.v3.bed")
write_tsv(to, fo, col_names = F)
# liftOver -minMatch=0.1 01.v3.bed chain 02.v4.raw.bed unmap2
# liftOver 02.v4.raw.bed chain 04.v4.unsorted.bed unmap4
# sortBed -i 04.v4.unsorted.bed > 05.v4.bed
#}}}

#{{{ pick mutants
#{{{ create exon range BED
fi = '~/projects/genome/data/Zmays_B73/50_annotation/15.tsv'
ti = read_tsv(fi)

tx = ti %>% filter(etype=='exon') %>% select(-etype)
gids_2exon = tx %>% count(gid) %>% filter(n>1) %>% pull(gid)
t_intr = tx %>% filter(gid %in% gids_2exon) %>%
    arrange(chrom, start) %>%
    group_by(gid) %>% nest() %>%
    mutate(intron=map(data, infer_intron)) %>%
    select(-data) %>% unnest() %>% mutate(etype = 'intron')

etype_map = c("CDS"='cds','five_prime_UTR'='utr5','three_prime_UTR'='utr3')
to = ti %>% filter(ttype == 'mRNA', etype != 'exon') %>%
    mutate(etype=etype_map[etype]) %>%
    select(gid, etype, chrom, start, end, srd)
to = rbind(to, t_intr)
to1 = to %>% filter(srd == '+') %>%
    arrange(chrom, start) %>%
    group_by(gid, etype) %>%
    mutate(eidx = 1:n()) %>% ungroup()
to2 = to %>% filter(srd == '-') %>%
    arrange(chrom, desc(start)) %>%
    group_by(gid, etype) %>%
    mutate(eidx = 1:n()) %>% ungroup()

to = rbind(to1, to2) %>% arrange(chrom, start) %>%
    mutate(start = start - 1) %>%
    select(chrom,start,end,srd,gid,etype,eidx)
fo = file.path(dirw, "12.genic.bed")
write_tsv(to, fo, col_names = F)

# intersectBed -a 11 -b 12 -wao > 13

etypes = c('cds','utr5','utr3','intron')
fi = file.path(dirw, "13.ovlp.bed")
ti = read_tsv(fi, col_names = c('mchrom','mbeg','mend','mid','sids',
                                'echrom','ebeg','eend','esrd','gid',
                                'etype','eidx','bp')) %>%
    select(mid, sids, gid, etype, eidx, bp) %>%
    mutate(etype = factor(etype, levels=etypes)) %>%
    arrange(mid, sids, gid, etype, desc(bp), eidx) %>%
    group_by(mid, sids, gid) %>%
    summarise(etype=etype[1], eidx=eidx[1]) %>% ungroup()
ti %>% mutate(genic = ifelse(gid == '.', T, F)) %>% count(genic)
mu = ti %>% filter(gid != '.')

fo = file.path(dirw, '15.mu.genic.tsv')
write_tsv(mu, fo)
#}}}

#{{{ add meta
fi = file.path(dirw, "15.mu.genic.tsv")
mu = read_tsv(fi)
n_mus = c('1','2','>=3')
tg1 = mu %>%
    #filter(etype != 'intron') %>%
    filter(etype == 'cds') %>%
    group_by(gid) %>%
    summarise(n_mu = n(),
              mid = str_c(mid, sep = "|", collapse = "|"),
              sids = str_c(sids, sep = "|", collapse = "|")) %>% ungroup() %>%
    mutate(n_mu_p = ifelse(n_mu >= 3, ">=3", n_mu)) %>%
    mutate(n_mu_p = factor(n_mu_p, levels = n_mus))

# add TF info
ft = file.path(dirp, 'data/B73/61_functional/06.tf.tsv')
tt = read_tsv(ft)
tts = tt %>% count(fam) %>% rename(fam_size=n)
tt = tt %>% inner_join(tts, by='fam') %>%
    group_by(fam) %>% mutate(fam_idx = 1:n()) %>% ungroup() %>%
    mutate(fam_idx_size = sprintf("%d/%d", fam_idx, fam_size)) %>%
    mutate(tf = T)
tg2 = tg1 %>% left_join(tt, by = 'gid') %>%
    replace_na(list(tf = F))

# add W22 gene ID
fi = '~/projects/genome/data/B73/gene_mapping/B73_W22.tsv'
ti = read_tsv(fi,col_names=c('chr1','start1','end1','gid1','chr2','start2','end2','gid2','src','type')) %>%
    transmute(gid=gid1, gid_W22=gid2, map_type=type)
ti2 = ti %>% group_by(gid) %>%
    summarise(gid_W22 = paste(gid_W22,collapse=','),
              map_type = paste(sort(unique(map_type)), collapse=',')) %>% ungroup()
tg3 = tg2 %>% left_join(ti2, by='gid')

# add W22-B73 gene model change
fi = '~/projects/wgc/data/05_stats/10.B73_W22.tsv'
ti = read_tsv(fi)
impacts = c("no_change","low",'modifier','moderate','high','non-syntenic')
tg4 = tg3 %>% left_join(ti, by = 'gid') %>%
    replace_na(list(syn='non-syntenic')) %>%
    mutate(impact = ifelse(syn=='syntenic', impact, syn)) %>%
    mutate(impact = factor(impact, levels = impacts)) %>%
    select(-po, -syn,-tid)
tg4 %>% filter(tf) %>% count(map_type, impact) %>% print(n=30)

# output
to = tg4
fo = file.path(dirw, "16.gene.mu.tsv")
write_tsv(to, fo)

to %>% count(n_mu_p)
to %>% count(impact)
to %>% count(tf)
low_impacts = c("no_change","low",'modifier','moderate')
to %>% filter(impact %in% low_impacts)
to %>% filter(tf, n_mu != '1', impact %in% low_impacts)
to %>% filter(tf, n_mu != '1', impact %in% low_impacts)
#}}}

#{{{ process W22 expression
fi1 = file.path(dirw, 'Samples_B73_Ref_HTseq.txt')
fi2 = file.path(dirw, 'Samples_W22_Ref_HTseq.txt')
ti = read_tsv(fi1) %>% gather(cond, counts, -Genes) %>%
    rename(gid=Genes) %>%
    filter(!str_detect(gid, '^_'))
#
ti2 = ti %>% mutate(cond = str_replace(cond, "_Ref_counts.txt", "")) %>%
    separate(cond, c('cond','ref'), sep='_') %>%
    separate(cond, c('Genotype',"Tissue",'Rep'), sep='-') %>%
    group_by(gid,Genotype,Tissue) %>%
    summarize(counts = sum(counts)) %>% ungroup() %>%
    group_by(Genotype, Tissue) %>%
    mutate(rpm = counts/sum(counts) * 1000000) %>%
    ungroup()
#
tis_map = c('A'='Anther','En'='Endosperm','Em'='Embryo','I'='Internode',
'IE'='Ear', 'L'='Leaf','L10'='Leaf10', 'R'='Root', 'SC'='Shoot', 'T'='Tassel')
ti3 = ti2 %>% filter(Genotype %in% c("B","W")) %>%
    filter(Tissue %in% c("En","Em","L10","SC","R","I")) %>%
    mutate(Tissue=tis_map[Tissue]) %>%
    mutate(cond = sprintf("%s_%s", Genotype, Tissue)) %>%
    select(gid,cond,rpm) %>%
    mutate(rpm = sprintf("%.01f", rpm)) %>%
    spread(cond, rpm)
te = ti3
#}}}

#{{{ further characterize all TF mutants
fi = file.path(dirw, "16.gene.mu.tsv")
ta = read_tsv(fi)

# add TF45 info
fi = '~/projects/grn/data/08_y1h/10.tsv'
ti = read_tsv(fi) %>% distinct(reg) %>% transmute(gid=reg, TF45=T)
ta2 = ta %>% left_join(ti, by='gid') %>% replace_na(list(TF45=F))
ta2 = ta2 %>% filter(tf | TF45) %>% select(-tf)

# add biomap support
fi = '~/projects/grn/data/14_eval_sum/02.valid.bm.spc.rds'
ti = readRDS(fi)$tf
ta3 = ta2 %>% left_join(ti, by=c('gid'='reg.gid')) %>%
    replace_na(list(n.tgt=0))

# add eQTL support
fi = '~/projects/grn/data/14_eval_sum/02.hs.tsv'
ti = read_tsv(fi)
ta4 = ta3 %>% left_join(ti, by=c('gid'='reg.gid')) %>%
    rename(eQTL=qtags, eQTL_n=n_qtag, eQTL_fc=fc, eQTL_grp_size=max.grp.size) %>%
    replace_na(list(eQTL=''))

ta4 = ta2
# add W22 expression
ta5 = ta4 %>% left_join(te, by='gid') %>%
    mutate(max_W22_exp=pmax(W_Embryo,W_Endosperm,W_Internode,W_Leaf10,W_Root,W_Shoot))

to = ta5
to = to %>% rename(gid_B73=gid) %>%
    select(-n_mu_p, -fam_idx_size, -ttype)
fo = file.path(dirw, '20.tf.tsv')
write_tsv(to, fo)

fi = file.path(dirw, '20.tf.tsv')
ti = read_tsv(fi) %>%
    left_join(gcfg$gene[,c('gid','note1','note2')], by=c("gid_B73"='gid')) %>%
    select(gid_B73, fam, fam_size, fam_idx, note1, note2, n_mu, mid, sids, everything())
fo = file.path(dirw, '21.tf.note.tsv')
write_tsv(ti, fo, na='')
#}}}

#{{{ select TF mutants
fi = file.path(dirw, '20.tf.tsv')
ti = read_tsv(fi)
ti
ti2 = ti %>% filter(n_mu >= 1, map_type == 'One-to-One', max_W22_exp >= 2)

gids1 = ti2 %>% filter(eQTL != '') %>%
    #filter(str_detect(eQTL,',')) %>%
    pull(gid_B73)
gids2 = ti2 %>% filter(n.tgt >= 3) %>% arrange(desc(n.tgt)) %>%
    filter(row_number() <= 20) %>% pull(gid_B73)
gids3 = ti2 %>%
    filter(fam %in% c("HSF","LBD","SBP","TCP","WRKY","MYB"), n_mu >= 2) %>%
    pull(gid_B73)
to1 = tibble(select='eQTL', gid=gids1)
to2 = tibble(select='biomAP', gid=gids2)
to3 = tibble(select='multi-fam', gid=gids3)
to = rbind(to1,to2,to3) %>% group_by(gid) %>%
    summarise(select=paste(select,collapse=',')) %>% ungroup()
to = ti %>% inner_join(to, by=c('gid_B73'='gid')) %>%
    select(gid_B73, select, everything()) %>%
    arrange(select, fam, gid_B73)

fo = file.path(dirw, '30.tf.selected.tsv')
write_tsv(to, fo)
#}}}

#{{{ stock order
fi = file.path(dirw, '30.tf.selected.tsv')
ti = read_tsv(fi)
fi = file.path(dirw, "15.mu.genic.tsv")
mu = read_tsv(fi)
etypes = c("cds")
fi = file.path(dirw, "31.erika.selected.xlsx")
tk = read_xlsx(fi) %>% select(gid=gid_B73, mu) %>%
    replace_na(list(mu='')) %>%
    mutate(mu = map(mu, str_split, "[\\|]")) %>% mutate(mu = map(mu, 1)) %>%
    unnest()

tm = ti %>% select(gid=gid_B73, select_reason=select, n_mu=n_mu) %>%
    inner_join(mu, by='gid') %>% filter(etype %in% etypes)
tm %>% count(gid, n_mu) %>% mutate(nd = n-n_mu) %>% pull(nd)
# use Erika's list to filter
tm.1 = tm %>% filter(!gid %in% tk$gid)
tm.2 = tm %>% inner_join(tk, by=c('gid'='gid','mid'='mu'))
tm2 = rbind(tm.1, tm.2)
#
tm2 = tm2 %>% mutate(etype=fct_relevel(etype, etypes)) %>%
    arrange(gid, etype, eidx) %>%
    group_by(gid) %>%
    filter(row_number() <= 3) %>% ungroup()
tm3 = tm2 %>% select(-n_mu) %>%
    mutate(sid = str_split(sids, ',')) %>%
    mutate(sid = map_chr(sid, 1)) %>%
    arrange(select_reason, gid)

fo = file.path(dirw, '32.gene.stocks.tsv')
write_tsv(tm3, fo)

tmo = tm3 %>% group_by(sid) %>%
    summarise(mid=paste(mid, collapse=","),
              gid=paste(gid, collapse=','),
              select_reason = paste(select_reason, collapse=',')) %>%
    ungroup() %>%
    arrange(sid)

fo = file.path(dirw, '34.stocks.tsv')
write_tsv(tmo, fo)
#}}}

#{{{ gather info from existing/collaborator stocks
fi = file.path(dirw, '34.stocks.tsv')
ti = read_tsv(fi)

# previous stocks
f_od = file.path(dirw, 'Springer_UniformMu_orders.xlsx')
t_od = read_xlsx(f_od, col_names = c("row",'stock','season')) %>%
    mutate(note = sprintf("%s_%s",season, row)) %>%
    group_by(stock) %>%
    summarise(note = paste(note, collapse=",")) %>%
    ungroup()

# mgc reply
fi = file.path(dirw, '40.stocks.misc.xlsx')
t_mgc = read_xlsx(fi, col_names = c("sid",'note')) %>%
    mutate(note = factor(note, unique(note))) %>%
    arrange(sid, desc(note)) %>%
    group_by(sid) %>% summarise(note = note[1]) %>% ungroup()

to = ti %>%
    left_join(t_mgc, by='sid') %>% rename(note1=note) %>%
    left_join(t_od, by=c('sid'='stock')) %>% rename(note2=note) %>%
    arrange(note1,note2,sid)
to %>% count(is.na(note1), is.na(note2))

fo = file.path(dirw, '41.stocks.field.tsv')
write_tsv(to, fo, na='')
#}}}
#}}}


#{{{ characterize picked mutants
fr = file.path(dirw, '50_field_ufmu_2019.xlsx')
tr = read_xlsx(fr) %>%
    select(Row,sid=Genotype,mids=Event,primer_f_seq,primer_r_seq) %>%
    fill(sid) %>%
    mutate(mid = str_split(mids, ',')) %>% unnest() %>% select(-mids) %>%
    select(Row, sid, mid, primer_f_seq, primer_r_seq)
tr %>% count(sid, mid) %>% filter(n>1)
tr %>% count(sid, mid) %>% count(sid) %>% filter(n>1)
tr %>% count(sid, mid) %>% count(mid) %>% filter(n>1)
fo = file.path(dirw, "50.um.rows.tsv")
write_tsv(tr, fo, na='')

fi = file.path(dirw, '32.gene.stocks.tsv')
ti = read_tsv(fi)

fu = file.path(dirw, "10.uniformmu.bed")
tu = read_tsv(fu, col_names=c("chrom",'start','end','mid','sids')) %>%
    mutate(sid = str_split(sids, ',')) %>% unnest() %>% select(-sids)

fg = '~/projects/genome/data/Zmays_B73/50_annotation/15.tsv'
tg = read_tsv(fg) %>% distinct(gid, chrom) %>% rename(gchrom=chrom)

tp = ti %>% select(mid,sid,gid,select_reason) %>%
    filter(mid %in% tr$mid) %>%
    inner_join(tg, by='gid') %>%
    inner_join(tu, by=c('mid','sid')) %>% filter(gchrom==chrom) %>%
    select(mid,sid,gid,select_reason,chrom,start,end)

mids_r = tr %>% distinct(mid) %>% pull(mid)
mids_i = ti %>% distinct(mid) %>% pull(mid)
mids_p = tp %>% distinct(mid) %>% pull(mid)
mids_r[!mids_r %in% mids_i]
mids_p[!mids_p %in% mids_r]
mids_r[!mids_r %in% mids_p]

# lift to W22 coordinates
tc = tp %>% distinct(chrom,start,end,mid)
fc = file.path(dirw, '52.sites.B73.bed')
write_tsv(tc, fc, col_names=F)
# liftOver 52.sites.B73.bed ~/projects/wgc/data/raw/W22_B73/10.B73_W22.chain 52.sites.W22.bed -minMatch=0.1 unmap.bed

fc = file.path(dirw, '52.sites.W22.bed')
tc = read_tsv(fc, col_names=c("chrom","start",'end','mid'))

# gene mapping
f_map = '~/projects/genome/data/Zmays_B73/gene_mapping/B73_W22.tsv'
t_map = read_tsv(f_map, col_names=c('chrom1','start1','end1','gid1','chrom2','start2','end2','gid2','src','type')) %>%
    select(gid1,gid2,type)
tm = tp %>% distinct(gid) %>% left_join(t_map, by = c("gid"='gid1'))
tm %>% count(type)
tm = tm %>% select(-type)

to = tp %>%
    select(mid,sid,select_reason,
           b.chrom=chrom,b.start=start,b.end=end,b.gid=gid) %>%
    left_join(tc, by='mid') %>%
    rename(w.chrom=chrom,w.start=start,w.end=end) %>%
    left_join(tm, by=c('b.gid'='gid')) %>% rename(w.gid=gid2)
to %>% print(width=Inf)
to %>% filter(is.na(w.chrom)) %>% print(width=Inf)

fo = file.path(dirw, '55.ufmu.tsv')
write_tsv(to, fo)
#}}}

#{{{ primer design - DNA
fm = file.path(dirw, '55.ufmu.tsv')
tm = read_tsv(fm) %>%
    mutate(start.i = w.start - 500, end.i = w.end + 500)

tp = tm %>% filter(!is.na(w.chrom)) %>%
    distinct(w.chrom, start.i, end.i)
fo = file.path(dirw, '61.bed')
write_tsv(tp, fo, col_names=F)
#fasta.py extract Zmays_W22 --tsv --padding 61.bed >61.tsv

ff = file.path(dirw, '61.tsv')
tf = read_tsv(ff, col_names=c('w.chrom','start.i','end.i','seq')) %>%
    mutate(start.i = start.i-1)
tm2 = tm %>% left_join(tf, by=c('w.chrom','start.i','end.i')) %>%
    mutate(start.t = w.start-start.i, end.t = w.end-start.i)

tp = tm2 %>% filter(!is.na(seq)) %>%
    distinct(mid,start.t,end.t,seq) %>%
    rename(pid=mid,start=start.t, end=end.t)
fp = file.path(dirw, '65.seq.tsv')
write_tsv(tp, fp)

# run_primer3.py ufmu 65.seq.tsv > 66.primer3.tsv

fr = file.path(dirw, '66.primer3.tsv')
tr = read_tsv(fr) %>% rename(mid=pid) %>%
    mutate(product_size = right.start - left.start)

trs = tibble(pseq=unique(c(tr$left.seq,tr$right.seq))) %>%
    mutate(pid = sprintf("ump%03d", 1:n())) %>% select(pid, pseq)
fo1 = file.path(dirw, '68.primer.seq.tsv')
write_tsv(trs, fo1)

tr2 = tr %>% select(-n_left,-n_right,-n_pair) %>%
    inner_join(trs, by=c('left.seq'='pseq')) %>% rename(left.id=pid) %>%
    inner_join(trs, by=c('right.seq'='pseq')) %>% rename(right.id=pid) %>%
    rename(mu.start=start, mu.end=end,
           left.primer.start=left.size, right.primer.size=right.size) %>%
    mutate(left.size=mu.start-left.start, right.size=right.start-mu.end) %>%
    select(mid,mu.start,mu.end,left.start,right.start,left.size, right.size,
        product_size, left.id,left.seq, right.id, right.seq)

# create UM-row table
fi = file.path(dirw, '50.um.rows.tsv')
ti = read_tsv(fi)

to = ti %>% left_join(tr2, by='mid')

fo = file.path(dirw, '69.um.rows.primer.tsv')
write_tsv(to, fo, na='')
#}}}


#{{{ 2019 summary
fi = file.path(dirw, '55.ufmu.tsv')
ti = read_tsv(fi)

ft = file.path(dirw, '71_ufmu_genotype.xlsx')
tt = read_xlsx(ft)

gts = c('hom','het','unk','w','na')
gt_col = pal_aaas()(5); names(gt_col) = gts
tp = ti %>% select(gid=b.gid,mid,sid,reason=select_reason) %>%
    full_join(tt, by=c('sid','mid')) %>%
    count(gid, reason, gt) %>%
    mutate(gt = factor(gt, levels = gts)) %>%
    spread(gt, n) %>%
    left_join(gs$tf_fam, by='gid') %>%
    left_join(gcfg$gene[,c('gid','note1','note2')], by='gid') %>%
    arrange(-hom, -het, -unk, -w, -na) %>%
    mutate(type = ifelse(!is.na(w), 'na', 'w')) %>%
    mutate(type = ifelse(!is.na(unk), 'unk', type)) %>%
    mutate(type = ifelse(!is.na(het), 'het', type)) %>%
    mutate(type = ifelse(!is.na(hom), 'hom', type)) %>%
    mutate(type = factor(type, levels=gts)) %>%
    mutate(col = gt_col[type]) %>% select(-type) %>%
    select(gid, hom, het, unk, w, na, reason, fam, note=note2, col) %>%
    print(n=50)

to = tp %>% select(-col)
fo = file.path(dirw, '75.summary2019.tsv')
write_tsv(to, fo, na='')

options(knitr.kable.NA = '')
x = tp %>% select(-col) %>%
    replace_na(list(hom='',het='',unk='',w='',na='')) %>%
    mutate(gid = cell_spec(gid, color=tp$col),
           hom = cell_spec(hom, color=tp$col),
           het = cell_spec(het, color=tp$col),
           unk = cell_spec(unk, color=tp$col),
           w = cell_spec(w, color=tp$col),
           na = cell_spec(na, color=tp$col)) %>%
    rename(`Gene ID`=gid, Reason=reason, unknown=unk, wildtype=w,
           `TF family`=fam) %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T,
        format.args = list(big.mark=",", na.encode=F)) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
    #column_spec(1, bold=T) %>%
fo = file.path(dirw, '75.summary2019.pdf')
x %>% save_kable(fo)
#}}}

#{{{  2019 sum for Erika
res = rnaseq_cpm('rn18g')
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
tiss23 = c("radicle_root", "seedling_root", "seedling_leaf", "seedling_meristem", "coleoptile", "auricle", "blade_leaf", "internode", "tassel", "sheath", "ear", "flag_leaf", "floret", "husk", "root", "silk", "spikelet", "tassel_stem", "endosperm14D", "kernel", "embryo", "seed_imbibed", "endosperm27D")
ta = tm_m %>% inner_join(th_m, by='SampleID') %>%
    filter(Genotype == 'B73') %>% select(gid, Tissue, CPM) %>%
    mutate(Tissue=factor(Tissue, levels=tiss23)) %>%
    spread(Tissue, CPM)

fi = file.path(dirw, '55.ufmu.tsv')
ti = read_tsv(fi)
ft = file.path(dirw, '71_ufmu_genotype.xlsx')
tt = read_xlsx(ft)
fe = file.path(dirw, '21.tf.note.tsv')
te = read_tsv(fe) %>% select(b.gid=gid_B73,fam, fam_size, note1, map_type,
                             starts_with("B_"), starts_with("W_")) %>%
    inner_join(ta, by=c('b.gid'='gid'))

to = ti %>%
    inner_join(te, by='b.gid') %>%
    inner_join(tt, by=c('mid', 'sid')) %>%
    select(-select_reason) %>%
    select(rid, mid, sid, gt, fam, fam_size, note1,
           b.gid, b.chrom, b.start, b.end, map_type, w.gid, everything()) %>%
    arrange(rid)

fo = file.path(dirw, '76.tf.info.tsv')
write_tsv(to, fo)
#}}}


