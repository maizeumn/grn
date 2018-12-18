# Y1H analysis

[Original dataset](y1h.targets.tsv):
* 44 TFs
* 123 targets
* 335 links (Y1H interactions)

These TFs and target genes were checked for overlap with [49 different GRNs](../10.dataset.tsv) (top 500,000 edges):
* [Overlap table](11.support.edges.tsv)
* 34 TFs were found to regulate 98 targets in at least one GRN
* 190 interactions
  * 144 found in only 1 network
  * 29 found in 2 networks
  * 10 found in 3 networks
  * 7 found in 4 networks
* only 10 out of the 190 interactions were also supported by Y1H
  * 8 found in 1 network
  * 2 found in 2 networks

Expression in [139 different maize tissues / developmental stages](https://github.com/orionzhou/rnaseq/blob/master/data/05_read_list/mec03.tsv) (460 total samples with an average of 3 replicates per tissue):
* normalized CPM (Counts Per Million) for [44 TFs](15.tf.cpm.tsv)
* [123 targets](15.targets.cpm.tsv)
