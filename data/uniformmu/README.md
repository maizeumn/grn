Maize UniformMu mutants characterization
================

* Raw GFF3 file mapping to B73 AGPv3 were obtained from maizeGDB
  * 79,407 insertion sites and related stocks information
* Insertion sites were lifted over to AGPv4
  * 77,085 sites
* The overlap of these sites with 40,239 AGPv4 gene models (longest alternative transcript) were checked:
  * 32,889 sites overlap w. at least one genes
  * 42,738 sites do not overlap w. any genic region (i.e., intergenic)

[Table listing all UniformMu insertions overlapping w. an AGPv4 gene](15.mu.genic.tsv)

* 15,597 genes have at least one insertion sites (CDS, UTR5, UTR3, intron)
* 13,631 have at least one exon insertion(s)
  * 7,905 have at least two exon insertions.
* 7,421 have at least one CDS insertion(s)
  * 2,681 have at least two CDS insertions.

Extra information were added to the table:
* The assignment to transcription factor family based on [PlantTFDB](http://planttfdb.cbi.pku.edu.cn/index.php?sp=Zma) was added.
* The W22 syntenic gene and mapping relationship (one-to-one, one-to-many, etc.) between B and W were added.
* Structural integrity of these B73 genes in the W22 genome were checked based on variants called by [whole genome comparison between B73 and W22](https://github.com/orionzhou/wgc/blob/master/Rmd/wgc.md).
* [Final table](16.gene.mu.tsv) lists all 13,631 genes with at least one exon insertions, with columns:
  * `gid_B73`: AGPv4 gene ID
  * `n_mu`: number exonic UniformMu insertion sites
  * `mid`, `sids`: mutant ID and stock ID
  * `fam`, `fam_size`, `fam_idx`: whether this gene is TF, the assigned TF family and member index within the family
  * `gid_W22`, `map_type`: W22 ortholog ID and mapping relationship
  * `impact`, `eff`: impact and effect of the change from B73 to W22

## Characterization of TF mutants

Starting from 470 TFs with at least one CDS insertion(s), with additional columns added:
* `TF45`: One of the TF45 involved in phenolic biosynthesis?
* `n.tgt`, `m.drc`, `tf.type`: support from the biomAP dataset
  * number of targets (`n.tgt`) supported by the biomAP dataset (in >=2 tissues)
  * `m.drc`: mean regulatory direction among all targets (-1 for consistent repression and 1 for consistent activation)
  * `tf.type`: repressor (m.drc < -0.8), activator (m.drc > 0.8), mixed
* `eQTL`: suported by at least one (out of three) previous eQTL studies (significantly co-regulate the targets of a trans-eQTL hotspot and located within 50Mb of its physical location)
* Expression in 10 tissues from B, P and W were added

[Table listing all 470 TFs with extra support information](20.tf.tsv)

## Selection of TF mutants

First, all 470 TFs were filtered by the following criteria:
* at least one exon insertion
* `One-to-One` mapping between the B73 and W22 gene model
* maximum W22 expression (RPM) > 2

Next, three different sets were obtained:
* 34 TFs with support from previous eQTL studies
* 20 TFs supported by the biomAP dataset with the most number of targets (53 - 286)
* 21 TFs in one of these families (`HSF,LBD,SBP,TCP,WRKY,MYB`) and at least 2 insertion sites

## Selection of UniformMu insertion sites and stocks

For TFs with >=3 insertions, different insertion sites were ranked by the rules below and only the top 3 insertions were selected:
* present in Erika's hand-selected list (3 TFs without good insertion sites were discarded)
* cds > utr5 > utr3 (based on insertion site location in [15.mu.genic.tsv](15.mu.genic.tsv))
* the more upstream (five prime) insertion site is given higher rank

For insertion sites with multiple stocks, the first stock is taken

[Final list of 67 selected TFs, 111 insertion sites and 106 stocks](https://docs.google.com/spreadsheets/d/1O4fHFqv-60JWQNa0E55ePWOd7gGJdj89HCXyt8e1nVA/edit?usp=sharing)

## UniformMu stock ordering

These stocks were checked for previous ordering history (from `Springer_UniformMu_orders.xlsx`):
* 11 stocks were ordered in a previous season
* 95 stocks were not ordered before and will be ordered

[95 stocks to order](https://docs.google.com/spreadsheets/d/1O4fHFqv-60JWQNa0E55ePWOd7gGJdj89HCXyt8e1nVA/edit?usp=sharing)

