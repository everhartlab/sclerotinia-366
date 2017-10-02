<!-- README.md is generated from README.Rmd. Please edit that file -->
Analysis of 366 *S. sclerotiorum* isolates
==========================================

[![Last-changedate](https://img.shields.io/badge/last%20change-2017--10--02-brightgreen.svg)](https://github.com/everhartlab/sclerotinia-366/commits/master) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.1-brightgreen.svg)](https://cran.r-project.org/) [![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/) [![Circle CI](https://circleci.com/gh/everhartlab/sclerotinia-366.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/everhartlab/sclerotinia-366)

This repository contains data, code, and a manuscript for analysis of 366 isolates of *Sclerotinia sclerotiorum* from the US and various countries around the world.

> Population structure and phenotypic variation of *Sclerotinia sclerotiorum* from dry bean in the United States
>
> Z. N. Kamvar, B. S. Amaradasa, R. Jhala, S. McCoy, J. Steadman, and S. E. Everhart

TOC
===

The analyses are arranged in the following order according to the [Makefile](Makefile):

1.  [table-1.md](results/table-1.md)
2.  [MCG-virulence.md](results/MCG-virulence.md)
3.  [locus-stats.md](results/locus-stats.md)
4.  [MLG-distribution.md](results/MLG-distribution.md)
5.  [mlg-mcg.md](results/mlg-mcg.md)
6.  [RDA-analysis.md](results/RDA-analysis.md)
7.  [pop-diff.md](results/pop-diff.md)
8.  [tree.md](results/tree.md)
9.  [wmn-differentiation.md](results/wmn-differentiation.md)
10. [by-year.md](results/by-year.md)

Reproducing the analysis
========================

Locally
-------

This project is controlled via a [Makefile](Makefile) which means that everything (analyses, tables, figures, the paper itself) is controlled via one command:

    make

This will bootstrap the installation (warning: it will update packages), process the data, perform the analyses, and compile the paper.

Required software:

-   GNU Make (If you're on Windows, you can use MinGW: <http://www.mingw.org/>)
-   [R (version 3.4.1 or greater)](https://r-project.org)
-   [LaTeX](https://www.latex-project.org/get)
-   [pandoc](http://pandoc.org/) (Note: pandoc ships with Rstudio)
-   [devtools](https://github.com/hadley/devtools#readme)

Docker
------

This repository contains a [Dockerfile](Dockerfile), which specifies the instructions to build a [docker](https://www.docker.com/) container. This is designed to capture the complete development environment of the analysis so that it can be accurately reproduced. The image is ~2.71Gb, so be sure that you have enough memory on your computer to run it.

To Install Docker, go here: <https://docs.docker.com/engine/installation/#desktop>. Once you have downloaded docker, you can either pull the container or build it. Pulling is by far the quickest way to do this. The docker container is located at <https://hub.docker.com/r/zkamvar/sclerotinia-366/>

> Note: both ways assume that you are in the analysis directory

### Pulling the container

    docker run --rm -it -v $(pwd):/analysis zkamvar/sclerotinia-366:latest bash

### Building the container

If you don't want to pull from docker hub, you can build the container.

``` sh
docker build .
```

This will bootstrap the environment to build the container. Once it's finished you should see:

    Successfully built 8e1e4cd82e19

where `8e1e4cd82e19` will be replaced by your hash. Now that things are built, you can run the analysis in the image with:

    docker run -it -v $(pwd):/analysis 8e1e4cd82e19 bash

### Running the Analysis

Once you are in the container, you can run the analysis, which is mapped to `analysis/`. The `make clean` command will wipe out all derivative files and the `make` command will generate everything. Note that this took almost 2 hours to run on my machine due to several bootstrapping processes.

    cd analysis/
    make clean
    make

Packages Used
-------------

``` r
options(width = 100)
imports <- packageDescription("WorldSclerotinia")$Imports
imports <- strsplit(imports, "[^A-z]*,\n")[[1]]
for (i in imports) suppressPackageStartupMessages(library(i, character.only = TRUE))
devtools::session_info()
#> Session info --------------------------------------------------------------------------------------
#>  setting  value                       
#>  version  R version 3.4.1 (2017-06-30)
#>  system   x86_64, darwin15.6.0        
#>  ui       X11                         
#>  language (EN)                        
#>  collate  en_US.UTF-8                 
#>  tz       America/Chicago             
#>  date     2017-10-02
#> Packages ------------------------------------------------------------------------------------------
#>  package     * version    date       source                                  
#>  ade4        * 1.7-8      2017-08-09 cran (@1.7-8)                           
#>  adegenet    * 2.1.0      2017-10-02 Github (thibautjombart/adegenet@8bc0ae0)
#>  agricolae   * 1.2-8      2017-09-12 cran (@1.2-8)                           
#>  AlgDesign     1.1-7.3    2014-10-15 CRAN (R 3.4.0)                          
#>  ape           4.1        2017-02-14 CRAN (R 3.4.0)                          
#>  assertr     * 2.0.2.2    2017-06-06 CRAN (R 3.4.0)                          
#>  assertthat    0.2.0      2017-04-11 CRAN (R 3.4.0)                          
#>  backports     1.1.1      2017-09-25 CRAN (R 3.4.2)                          
#>  base        * 3.4.1      2017-07-07 local                                   
#>  bindr         0.1        2016-11-13 CRAN (R 3.4.0)                          
#>  bindrcpp      0.2        2017-06-17 CRAN (R 3.4.0)                          
#>  bookdown    * 0.5        2017-08-20 CRAN (R 3.4.1)                          
#>  boot          1.3-20     2017-07-30 CRAN (R 3.4.1)                          
#>  broom         0.4.2      2017-02-13 CRAN (R 3.4.0)                          
#>  cellranger    1.1.0      2016-07-27 CRAN (R 3.4.0)                          
#>  cluster       2.0.6      2017-03-16 CRAN (R 3.4.0)                          
#>  coda          0.19-1     2016-12-08 CRAN (R 3.4.0)                          
#>  colorspace    1.3-3      2017-08-16 R-Forge (R 3.4.1)                       
#>  combinat      0.0-8      2012-10-29 CRAN (R 3.4.0)                          
#>  compiler      3.4.1      2017-07-07 local                                   
#>  cowplot     * 0.8.0.9000 2017-08-28 Github (wilkelab/cowplot@a0b419e)       
#>  datasets    * 3.4.1      2017-07-07 local                                   
#>  DBI           0.7        2017-06-18 CRAN (R 3.4.0)                          
#>  deldir        0.1-14     2017-04-22 CRAN (R 3.4.0)                          
#>  devtools      1.13.3     2017-08-02 CRAN (R 3.4.1)                          
#>  digest        0.6.12     2017-01-27 CRAN (R 3.4.0)                          
#>  dplyr       * 0.7.4      2017-09-28 CRAN (R 3.4.1)                          
#>  evaluate      0.10.1     2017-06-24 CRAN (R 3.4.1)                          
#>  expm          0.999-2    2017-03-29 CRAN (R 3.4.0)                          
#>  ezknitr     * 0.6        2016-09-16 CRAN (R 3.4.0)                          
#>  fastmatch     1.1-0      2017-01-28 CRAN (R 3.4.0)                          
#>  forcats       0.2.0      2017-01-23 CRAN (R 3.4.0)                          
#>  foreign       0.8-69     2017-06-21 CRAN (R 3.4.0)                          
#>  gdata         2.18.0     2017-06-06 CRAN (R 3.4.0)                          
#>  ggcompoplot * 0.1.0      2017-10-02 Github (zkamvar/ggcompoplot@bcf007d)    
#>  ggforce     * 0.1.1      2016-11-28 CRAN (R 3.4.0)                          
#>  ggplot2     * 2.2.1      2016-12-30 CRAN (R 3.4.0)                          
#>  ggraph      * 1.0.0      2017-02-24 CRAN (R 3.4.0)                          
#>  ggrepel     * 0.6.12     2017-09-30 Github (slowkow/ggrepel@fd15d0a)        
#>  ggridges    * 0.4.1      2017-09-15 cran (@0.4.1)                           
#>  ggtree      * 1.9.4      2017-10-02 Github (GuangchuangYu/ggtree@07063f9)   
#>  glue          1.1.1      2017-06-21 CRAN (R 3.4.0)                          
#>  gmodels       2.16.2     2015-07-22 CRAN (R 3.4.0)                          
#>  graphics    * 3.4.1      2017-07-07 local                                   
#>  grDevices   * 3.4.1      2017-07-07 local                                   
#>  grid          3.4.1      2017-07-07 local                                   
#>  gridExtra     2.3        2017-09-09 CRAN (R 3.4.1)                          
#>  gtable        0.2.0      2016-02-26 CRAN (R 3.4.0)                          
#>  gtools        3.5.0      2015-05-29 CRAN (R 3.4.0)                          
#>  haven         1.1.0      2017-07-09 CRAN (R 3.4.1)                          
#>  hms           0.3        2016-11-22 CRAN (R 3.4.0)                          
#>  htmltools     0.3.6      2017-04-28 CRAN (R 3.4.0)                          
#>  htmlwidgets   0.9        2017-07-10 cran (@0.9)                             
#>  httpuv        1.3.5      2017-07-04 CRAN (R 3.4.1)                          
#>  httr          1.3.1      2017-08-20 cran (@1.3.1)                           
#>  huxtable    * 0.3.1.9000 2017-10-02 Github (hughjonesd/huxtable@4ee68c5)    
#>  igraph      * 1.1.2      2017-07-21 cran (@1.1.2)                           
#>  jsonlite      1.5        2017-06-01 CRAN (R 3.4.0)                          
#>  KernSmooth  * 2.23-15    2015-06-29 CRAN (R 3.4.0)                          
#>  klaR          0.6-12     2014-08-06 CRAN (R 3.4.0)                          
#>  knitr       * 1.17       2017-08-10 cran (@1.17)                            
#>  lattice     * 0.20-35    2017-03-25 CRAN (R 3.4.0)                          
#>  lazyeval      0.2.0      2016-06-12 CRAN (R 3.4.0)                          
#>  LearnBayes    2.15       2014-05-29 CRAN (R 3.4.0)                          
#>  lubridate     1.6.0      2016-09-13 CRAN (R 3.4.0)                          
#>  magrittr      1.5        2014-11-22 CRAN (R 3.4.0)                          
#>  MASS          7.3-47     2017-04-21 CRAN (R 3.4.0)                          
#>  Matrix        1.2-11     2017-08-16 CRAN (R 3.4.1)                          
#>  memoise       1.1.0      2017-04-21 CRAN (R 3.4.0)                          
#>  methods     * 3.4.1      2017-07-07 local                                   
#>  mgcv          1.8-22     2017-09-19 CRAN (R 3.4.2)                          
#>  mime          0.5        2016-07-07 CRAN (R 3.4.0)                          
#>  mnormt        1.5-5      2016-10-15 CRAN (R 3.4.0)                          
#>  modelr        0.1.1      2017-07-24 CRAN (R 3.4.1)                          
#>  munsell       0.4.3      2016-02-13 CRAN (R 3.4.0)                          
#>  nlme          3.1-131    2017-02-06 CRAN (R 3.4.0)                          
#>  parallel      3.4.1      2017-07-07 local                                   
#>  pegas         0.10       2017-05-03 CRAN (R 3.4.0)                          
#>  permute     * 0.9-4      2016-09-09 CRAN (R 3.4.0)                          
#>  phangorn      2.2.0      2017-04-03 CRAN (R 3.4.0)                          
#>  pkgconfig     2.0.1      2017-03-21 CRAN (R 3.4.0)                          
#>  plyr          1.8.4      2016-06-08 CRAN (R 3.4.0)                          
#>  poppr       * 2.5.0      2017-09-11 CRAN (R 3.4.1)                          
#>  psych         1.7.8      2017-09-09 CRAN (R 3.4.1)                          
#>  purrr       * 0.2.3      2017-08-02 CRAN (R 3.4.1)                          
#>  quadprog      1.5-5      2013-04-17 CRAN (R 3.4.0)                          
#>  R6            2.2.2      2017-06-17 cran (@2.2.2)                           
#>  Rcpp          0.12.13    2017-09-28 CRAN (R 3.4.1)                          
#>  readr       * 1.1.1      2017-05-16 CRAN (R 3.4.0)                          
#>  readxl      * 1.0.0      2017-04-18 CRAN (R 3.4.0)                          
#>  reshape2      1.4.2      2016-10-22 CRAN (R 3.4.0)                          
#>  rlang         0.1.2      2017-08-09 cran (@0.1.2)                           
#>  rmarkdown     1.6        2017-06-15 cran (@1.6)                             
#>  rprojroot     1.2        2017-01-16 CRAN (R 3.4.0)                          
#>  rticles     * 0.4.1      2017-10-02 Github (rstudio/rticles@4111a39)        
#>  rvcheck       0.0.9      2017-07-10 cran (@0.0.9)                           
#>  rvest         0.3.2      2016-06-17 CRAN (R 3.4.0)                          
#>  scales        0.5.0.9000 2017-08-28 Github (hadley/scales@d767915)          
#>  seqinr        3.4-5      2017-08-01 CRAN (R 3.4.1)                          
#>  shiny         1.0.5      2017-08-23 cran (@1.0.5)                           
#>  sp            1.2-5      2017-06-29 CRAN (R 3.4.1)                          
#>  spdep         0.6-15     2017-09-01 CRAN (R 3.4.1)                          
#>  splines       3.4.1      2017-07-07 local                                   
#>  stats       * 3.4.1      2017-07-07 local                                   
#>  stringi       1.1.5      2017-04-07 CRAN (R 3.4.0)                          
#>  stringr       1.2.0      2017-02-18 CRAN (R 3.4.0)                          
#>  tibble      * 1.3.4      2017-08-22 cran (@1.3.4)                           
#>  tidyr       * 0.7.1      2017-09-01 CRAN (R 3.4.1)                          
#>  tidyverse   * 1.1.1      2017-01-27 CRAN (R 3.4.0)                          
#>  tools         3.4.1      2017-07-07 local                                   
#>  treeio      * 1.1.2      2017-10-02 Github (GuangchuangYu/treeio@b6ae142)   
#>  tweenr        0.1.5      2016-10-10 CRAN (R 3.4.0)                          
#>  udunits2      0.13       2016-11-17 CRAN (R 3.4.0)                          
#>  units         0.4-6      2017-08-27 CRAN (R 3.4.1)                          
#>  utils       * 3.4.1      2017-07-07 local                                   
#>  vegan       * 2.4-4      2017-08-24 cran (@2.4-4)                           
#>  viridis     * 0.4.0      2017-03-27 CRAN (R 3.4.0)                          
#>  viridisLite * 0.2.0      2017-03-24 CRAN (R 3.4.0)                          
#>  visNetwork  * 2.0.1      2017-07-30 cran (@2.0.1)                           
#>  withr         2.0.0      2017-07-28 CRAN (R 3.4.1)                          
#>  xml2          1.1.1      2017-01-24 CRAN (R 3.4.0)                          
#>  xtable        1.8-2      2016-02-05 CRAN (R 3.4.0)                          
#>  yaml          2.1.14     2016-11-12 CRAN (R 3.4.0)
```
