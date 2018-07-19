require(rmarkdown) 
require(knitr) 
require(tidyverse) 
require(kableExtra) 
dirm = '~/projects/maize.grn/Rmd'
getwd()

rmarkdown::render(file.path(dirm, "grn.Rmd"))
