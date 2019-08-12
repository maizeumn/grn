source("functions.R")
diri = '~/projects/rnaseq'
dirw = file.path(dird, '12_tables')

#{{{ tables
#{{{ GRN sum
tp = t_cfg %>% select(net_type,nid,study,note,sample_size,lgd,col) %>%
    mutate(study=str_replace(str_to_title(study), '(\\d+)$', ' et al. \\1')) %>%
    arrange(net_type, nid) %>%
    select(-nid)
x = tp %>% select(-col) %>%
    mutate(net_type = cell_spec(net_type, color=tp$col),
           lgd = cell_spec(lgd, color=tp$col)) %>%
    rename(`Network type`=net_type, Study=study, Note=note, N=sample_size, `Network label`=lgd) %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T,
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left') %>%
    #column_spec(1, bold=T) %>%
    collapse_rows(columns=c(1,2), latex_hline="major", valign="top")
fo = file.path(dirw, '01.t1.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t1.rds')
saveRDS(x, file=fo)
#}}}

#{{{ TF mutant
fi = '~/projects/barn/data/01.cfg.xlsx'
tp = read_xlsx(fi, 'mutants') %>%
    mutate(accession=str_sub(accession, 1, 11)) %>%
    replace_na(list(author='')) %>%
    mutate(author=str_to_title(author)) %>%
    select(`TF alias`=gene_alias, `TF name`=gene_name, `TF ID`=gene_id,
           Study=author, Accession=accession, Tissue=tissue, N=n)
x = tp %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirw, '01.t2.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t2.rds')
saveRDS(x, file=fo)
#}}}

#{{{ natural variation DEGs
fi = '~/projects/master.xlsx'
th = read_xlsx(fi, 'barn')
fi = file.path(dirw, '../06_deg/all.rds')
t_ds = readRDS(fi)

tp = t_ds %>% count(yid, cond, group1, group2, DE) %>% spread(DE, n) %>%
    rename(`non-DE` = non_DE) %>%
    mutate(cond = str_replace(cond, "_", " ")) %>%
    filter(!str_detect(group2,'x'))
tp = th %>% select(yid, author, study) %>%
    filter(! yid %in% c("rn15e",'rn18b')) %>%
#    mutate(author=str_replace(str_to_title(author), '(\\d+)$', ' et al. \\1')) %>%
    mutate(author=ifelse(yid=='rn14f','waters2017',author)) %>%
    mutate(study=ifelse(yid=='rn14f','stress cis-trans',study)) %>%
    inner_join(tp, by='yid') %>%
    mutate(cond=ifelse(yid=='rn14f',str_replace(cond,'^seedling','leaf3'), cond)) %>%
    select(-yid) %>%
    unite('contrast', group1, group2, sep=' vs ') %>%
    rename(condition = cond)
x = tp %>%
    mutate(author = str_to_title(author)) %>%
    kable(format='latex', escape=F, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 8, position='left') %>%
    #column_spec(1, italic = T) %>%
    collapse_rows(columns=c(1,2), latex_hline="major", valign="top")
fo = file.path(dirw, '01.t3.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t3.rds')
saveRDS(x, file=fo)
#}}}

#|{{{ trans hotspots
fi = file.path(dird, '14_eval_sum/33.hs.tsv')
tp = read_tsv(fi) %>%
    replace_na(list(reg.note='',txt='')) %>%
    mutate(reg.note = str_replace(reg.note, " *\\[.*\\]", '')) %>%
    mutate(reg.note = str_replace(reg.note, ";.*$", '')) %>%
    select(ID=reg.gid, `Support eQTL study`=qtags, `Support GRN`=studies, `TF Annotation`=reg.note, `Target enrichment`=txt)
x = tp %>%
    kable(format='latex', escape=F, longtable=T, booktabs=T,
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=T, font_size=8, position='left') %>%
    #column_spec(1:3, width = "1.85cm") %>%
    column_spec(5, width = "5cm")
fo = file.path(dirw, '01.t9.pdf')
x %>% save_kable(fo)
fo = file.path(dirw, '01.t9.rds')
saveRDS(x, file=fo)
#}}}
#}}}


