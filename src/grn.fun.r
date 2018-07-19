#{{{ load required libraries, define common variables
require(tidyverse)
require(igraph)
require(grid)
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
options(stringsAsFactors = FALSE)
#}}}


dirg = file.path(Sys.getenv("genome"), "Zmays_v4/TF")
dirp = '~/projects/maize.grn'
dird = file.path(dirp, 'data')
dira = file.path(dirp, 'analysis')
dirr = file.path(dirp, 'Rmd')


radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
