# maize.grn

This repo hosts code and data related to maize gene regulatory network study.

This repo provides a collection of published maize transcription factor functional studies, with focus on Chip-Seq and knockout mutant RNA-Seq datasets.

Raw sequencing reads were downloaded from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and mapped to the [maize B73 AGP_v4 genome](http://plants.ensembl.org/Zea_mays/Info/Index) using [STAR](https://github.com/alexdobin/STAR).  Uniquely mapped reads were assigned to and counted for the 46,117 reference gene models ([Ensembl Plants v37](http://plants.ensembl.org/Zea_mays/Info/Index)) using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).  Raw read counts were then normalized using the [TMM normalization approach](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to give CPMs (Counts Per Million reads) and then further normalized by gene CDS lengths to give FPKM (Fragments Per Kilobase of exon per Million reads) values.  Hierarchical clustering and principal component analysis were used to explore sample clustering pattern.

| Study | TF | Tissues | Targets | Note | Results |
| ----- | ------------- | ------ | ------------------- | --- | ---- |
| [Bolduc2012 KN1](https://paperpile.com/view/0820167c-9a22-0659-b253-797c30cf7ad8) | KNOTTED1 | meristem (ear, tassel, leaf) | 273, 228, 271 |  | [link]() |
| [Morohashi2012 P1](https://paperpile.com/view/7a0a1a96-7221-0fd8-a9f8-8154bb44b9d0) | Pericarp Color1 |  pericarps,  floral | 1,500 | | [link]()|
| [Eveland2014 RA1](https://paperpile.com/view/332f982e-3361-0d67-a640-57ad1fe1cd95) | RAMOSA1 | ear | 207 | | [link]() |
| [Pautler2015 FEA1](https://paperpile.com/view/1d33d29b-c855-080f-95ea-b1a6c6554b50) | fasciated ear4 | shoot meristem | 98 |  | [link]() |
| [Li2015 O2](https://paperpile.com/view/d4421338-7ca5-045b-a047-931256301428) | Opaque2 | endosperm | 35 (24) |  | [link]() |

