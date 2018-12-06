# Y1H analysis

[Original dataset](y1h.targets.tsv):
* 44 TFs
* 123 targets
* 335 links (Y1H interactions)

These TFs and target genes were checked for overlap with [49 different GRNs](../10.dataset.tsv) (top 500,000 edges):
* [Overlap table](11.support.edges.tsv)
* 34 TFs were found to regulate 98 targets in at least one GRN
* 190 interactions: only 10 were also supported by Y1H
  * 144 found in only 1 network
  * 29 found in 2 networks
  * 10 found in 3 networks
  * 7 found in 4 networks
