<!-- README.md is generated from README.Rmd. Please edit that file -->
Analysis of 366 *S. sclerotiorum* isolates
==========================================

[![Last-changedate](https://img.shields.io/badge/last%20change-2017--12--07-2F4096.svg)](https://github.com/everhartlab/sclerotinia-366/commits/master) [![Licence](https://img.shields.io/badge/license-MIT%20License-2F4096.svg)](http://choosealicense.com/licenses/mit/)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.2-2F4096.svg)](https://cran.r-project.org/) [![Circle CI](https://circleci.com/gh/everhartlab/sclerotinia-366.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/everhartlab/sclerotinia-366)

This repository contains data, code, and a manuscript for analysis of 366 isolates of *Sclerotinia sclerotiorum* from the US and various countries around the world.

Citations
=========

Preprint
--------

> Kamvar ZN, Amaradasa BS, Jhala R, McCoy S, Steadman JR, Everhart SE. (2017) Population structure and phenotypic variation of *Sclerotinia sclerotiorum* from dry bean (*Phaseolus vulgaris*) in the United States. *PeerJ* 5:e4152 <https://doi.org/10.7717/peerj.4152>

Data and Code
-------------

> Kamvar, Z. N., Amaradasa, B. S., Jhala, R., McCoy, S., Steadman, J. R., & Everhart, S. E. (2017, November). Data and analysis for population structure and phenotypic variation of *Sclerotinia sclerotiorum* from dry bean (*Phaseolus vulgaris*) in the United States. <https://doi.org/10.17605/OSF.IO/K8WTM>

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
11. [compare-aldrich-wolfe.md](results/compare-aldrich-wolfe.md)

Analysis
========

Background
----------

The analysis is controlled via two docker containers:

1.  [sclerotinia-366-dependencies](https://hub.docker.com/r/everhartlab/sclerotinia-366-dependencies/): defines the complete software environment used, built on top of the [rocker/verse:3.4.2](https://hub.docker.com/r/rocker/verse) container. **(See the [Dockerfile](https://osf.io/mwv6d/))**
2.  [sclerotinia-366](https://hub.docker.com/r/everhartlab/sclerotinia-366/) is built on top of the above container and contains the results of the analysis. **(See the [Dockerfile](https://osf.io/tvfju/))**

The `sclerotinia-366-dependencies` container is regularly rebuilt on docker hub whenever `rocker/verse:3.4.2` updates and the `sclerotinia-366` container is rebuilt on [CircleCI](https://circleci.com/gh/everhartlab/sclerotinia-366) whenever the git repository is updated.

As of this writing, the containers are up to date with R version 3.4.2 and packages downloaded from the [MRAN snapshot on 2017-10-31](https://mran.microsoft.com/snapshot/2017-10-31/).

Jump to [Reproduction via Docker](#reproduction-via-docker) or [Reproduction Locally](#locally).

------------------------------------------------------------------------

Reproduction via Docker
-----------------------

This repository contains a [Dockerfile](Dockerfile), which specifies the instructions to build a [docker](https://www.docker.com/) container. This is designed to capture the complete development environment of the analysis so that it can be accurately reproduced. The image is ~3.21Gb, so be sure that you have enough memory on your computer to run it.

To Install Docker, go here: <https://docs.docker.com/engine/installation/#desktop>. Once you have downloaded docker, you can either pull the container or build it. Pulling is by far the quickest way to do this. The docker container is located at <https://hub.docker.com/r/everhartlab/sclerotinia-366/>

### RStudio Server

To run the docker container, make sure you have Docker running, open your terminal and type:

``` bash
docker run --name ssc --rm -dp 8787:8787 -e ROOT=TRUE everhartlab/sclerotinia-366:latest
```

This will first check to make sure you have the container on your machine. If you don't, Docker will automatically download it for you. It will spin up the Docker container on your machine, give it the name "ssc", and expose it to port 8787. You can open your browser and type `localhost:8787`, and an instance of Rstudio server will appear. Sign in with the following credentials:

-   username: rstudio
-   password: rstudio.

Since the files in `/analysis` are write-protected, if you wanted to explore, you should copy the directory to your current working space:

-   in the R console type: `system("cp -R /analysis .")`.
-   open `/analysis` and double click on znk\_analysis.Rproj

From here you can re-run the analyses to your heart's content. **Don't forget to stop the container when you are finished:**

``` bash
docker stop ssc
```

### Building the container locally

If you don't want to pull from docker hub, you can build the container locally. First, download the repository

``` bash
git clone https://github.com/everhartlab/sclerotinia-366.git
cd sclerotinia-366/
docker build -t sclerotinia-366 .
```

Now that things are built, you can run the analysis in the image with:

    docker run -it sclerotinia-366 bash

### Running the Analysis

Once you are in the container, you can run the analysis, which is mapped to `analysis/`. The `make clean` command will wipe out all derivative files and the `make` command will generate everything. Note that this took almost 2 hours to run on my machine due to several bootstrapping processes.

    cd analysis/
    make clean
    make

Locally
-------

This project is controlled via a [Makefile](Makefile) which means that everything (analyses, tables, figures, the paper itself) is controlled via one command:

    make

This will bootstrap the installation (warning: it will update packages), process the data, perform the analyses, and compile the paper.

> Note: This analysis is only guaranteed to work with [the stated software environment](#packages-used).

Required software:

-   GNU Make (If you're on Windows, you can use MinGW: <http://www.mingw.org/>)
-   [R (version 3.4.1 or greater)](https://r-project.org)
-   [LaTeX](https://www.latex-project.org/get)
-   [pandoc](http://pandoc.org/) (Note: pandoc ships with Rstudio)
-   [devtools](https://github.com/hadley/devtools#readme)

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
#>  version  R version 3.4.2 (2017-09-28)
#>  system   x86_64, linux-gnu           
#>  ui       X11                         
#>  language (EN)                        
#>  collate  en_US.UTF-8                 
#>  tz       UTC                         
#>  date     2017-12-07
#> Packages ------------------------------------------------------------------------------------------
#>  package     * version date       source                               
#>  ade4        * 1.7-8   2017-08-09 cran (@1.7-8)                        
#>  adegenet    * 2.1.0   2017-10-12 cran (@2.1.0)                        
#>  agricolae   * 1.2-8   2017-09-12 cran (@1.2-8)                        
#>  AlgDesign     1.1-7.3 2014-10-15 cran (@1.1-7.3)                      
#>  ape           5.0     2017-10-30 cran (@5.0)                          
#>  assertr     * 2.0.2.2 2017-06-06 cran (@2.0.2.2)                      
#>  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.2)                       
#>  backports     1.1.1   2017-09-25 CRAN (R 3.4.2)                       
#>  base        * 3.4.2   2017-11-01 local                                
#>  bindr         0.1     2016-11-13 CRAN (R 3.4.2)                       
#>  bindrcpp      0.2     2017-06-17 CRAN (R 3.4.2)                       
#>  bookdown    * 0.5     2017-08-20 CRAN (R 3.4.2)                       
#>  boot          1.3-20  2017-07-30 cran (@1.3-20)                       
#>  broom         0.4.2   2017-02-13 CRAN (R 3.4.2)                       
#>  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.2)                       
#>  cluster       2.0.6   2017-03-16 CRAN (R 3.4.2)                       
#>  coda          0.19-1  2016-12-08 cran (@0.19-1)                       
#>  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.2)                       
#>  combinat      0.0-8   2012-10-29 cran (@0.0-8)                        
#>  compiler      3.4.2   2017-11-01 local                                
#>  cowplot     * 0.8.0   2017-07-30 cran (@0.8.0)                        
#>  datasets    * 3.4.2   2017-11-01 local                                
#>  deldir        0.1-14  2017-04-22 cran (@0.1-14)                       
#>  devtools      1.13.3  2017-08-02 CRAN (R 3.4.2)                       
#>  digest        0.6.12  2017-01-27 CRAN (R 3.4.2)                       
#>  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)                       
#>  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.2)                       
#>  expm          0.999-2 2017-03-29 cran (@0.999-2)                      
#>  ezknitr     * 0.6     2016-09-16 cran (@0.6)                          
#>  fastmatch     1.1-0   2017-01-28 cran (@1.1-0)                        
#>  forcats       0.2.0   2017-01-23 CRAN (R 3.4.2)                       
#>  foreign       0.8-69  2017-06-21 CRAN (R 3.4.2)                       
#>  gdata         2.18.0  2017-06-06 cran (@2.18.0)                       
#>  ggcompoplot * 0.1.0   2017-11-15 Github (zkamvar/ggcompoplot@bcf007d) 
#>  ggforce     * 0.1.1   2016-11-28 cran (@0.1.1)                        
#>  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.2)                       
#>  ggraph      * 1.0.0   2017-02-24 cran (@1.0.0)                        
#>  ggrepel     * 0.7.0   2017-09-29 cran (@0.7.0)                        
#>  ggridges    * 0.4.1   2017-09-15 cran (@0.4.1)                        
#>  ggtree      * 1.9.4   2017-11-15 Github (GuangchuangYu/ggtree@07063f9)
#>  glue          1.2.0   2017-10-29 CRAN (R 3.4.2)                       
#>  gmodels       2.16.2  2015-07-22 cran (@2.16.2)                       
#>  graphics    * 3.4.2   2017-11-01 local                                
#>  grDevices   * 3.4.2   2017-11-01 local                                
#>  grid          3.4.2   2017-11-01 local                                
#>  gridExtra     2.3     2017-09-09 CRAN (R 3.4.2)                       
#>  gtable        0.2.0   2016-02-26 CRAN (R 3.4.2)                       
#>  gtools        3.5.0   2015-05-29 cran (@3.5.0)                        
#>  haven         1.1.0   2017-07-09 CRAN (R 3.4.2)                       
#>  hms           0.3     2016-11-22 CRAN (R 3.4.2)                       
#>  htmltools     0.3.6   2017-04-28 CRAN (R 3.4.2)                       
#>  htmlwidgets   0.9     2017-07-10 CRAN (R 3.4.2)                       
#>  httpuv        1.3.5   2017-07-04 CRAN (R 3.4.2)                       
#>  httr          1.3.1   2017-08-20 CRAN (R 3.4.2)                       
#>  huxtable    * 1.1.0   2017-10-20 cran (@1.1.0)                        
#>  igraph      * 1.1.2   2017-07-21 CRAN (R 3.4.2)                       
#>  jsonlite      1.5     2017-06-01 CRAN (R 3.4.2)                       
#>  KernSmooth  * 2.23-15 2015-06-29 cran (@2.23-15)                      
#>  klaR          0.6-12  2014-08-06 cran (@0.6-12)                       
#>  knitr       * 1.17    2017-08-10 CRAN (R 3.4.2)                       
#>  lattice     * 0.20-35 2017-03-25 CRAN (R 3.4.2)                       
#>  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2)                       
#>  LearnBayes    2.15    2014-05-29 cran (@2.15)                         
#>  lubridate     1.7.0   2017-10-29 CRAN (R 3.4.2)                       
#>  magrittr      1.5     2014-11-22 CRAN (R 3.4.2)                       
#>  MASS          7.3-47  2017-04-21 CRAN (R 3.4.2)                       
#>  Matrix        1.2-11  2017-08-16 CRAN (R 3.4.2)                       
#>  memoise       1.1.0   2017-04-21 CRAN (R 3.4.2)                       
#>  methods     * 3.4.2   2017-11-01 local                                
#>  mgcv          1.8-22  2017-09-19 CRAN (R 3.4.2)                       
#>  mime          0.5     2016-07-07 CRAN (R 3.4.2)                       
#>  mnormt        1.5-5   2016-10-15 CRAN (R 3.4.2)                       
#>  modelr        0.1.1   2017-07-24 CRAN (R 3.4.2)                       
#>  munsell       0.4.3   2016-02-13 CRAN (R 3.4.2)                       
#>  nlme          3.1-131 2017-02-06 CRAN (R 3.4.2)                       
#>  parallel      3.4.2   2017-11-01 local                                
#>  pegas         0.10    2017-05-03 cran (@0.10)                         
#>  permute     * 0.9-4   2016-09-09 cran (@0.9-4)                        
#>  phangorn      2.2.0   2017-04-03 cran (@2.2.0)                        
#>  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.2)                       
#>  plyr          1.8.4   2016-06-08 CRAN (R 3.4.2)                       
#>  poppr       * 2.5.0   2017-09-11 cran (@2.5.0)                        
#>  psych         1.7.8   2017-09-09 CRAN (R 3.4.2)                       
#>  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2)                       
#>  quadprog      1.5-5   2013-04-17 cran (@1.5-5)                        
#>  R6            2.2.2   2017-06-17 CRAN (R 3.4.2)                       
#>  Rcpp          0.12.13 2017-09-28 CRAN (R 3.4.2)                       
#>  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.2)                       
#>  readxl      * 1.0.0   2017-04-18 CRAN (R 3.4.2)                       
#>  reshape2      1.4.2   2016-10-22 CRAN (R 3.4.2)                       
#>  rlang         0.1.2   2017-08-09 CRAN (R 3.4.2)                       
#>  rmarkdown     1.6     2017-06-15 CRAN (R 3.4.2)                       
#>  rprojroot     1.2     2017-01-16 CRAN (R 3.4.2)                       
#>  rticles     * 0.4.1   2017-11-15 Github (rstudio/rticles@4111a39)     
#>  rvcheck       0.0.9   2017-07-10 cran (@0.0.9)                        
#>  rvest         0.3.2   2016-06-17 CRAN (R 3.4.2)                       
#>  scales        0.5.0   2017-08-24 CRAN (R 3.4.2)                       
#>  seqinr        3.4-5   2017-08-01 cran (@3.4-5)                        
#>  shiny         1.0.5   2017-08-23 CRAN (R 3.4.2)                       
#>  sp            1.2-5   2017-06-29 CRAN (R 3.4.2)                       
#>  spdep         0.6-15  2017-09-01 cran (@0.6-15)                       
#>  splines       3.4.2   2017-11-01 local                                
#>  stats       * 3.4.2   2017-11-01 local                                
#>  stringi       1.1.5   2017-04-07 CRAN (R 3.4.2)                       
#>  stringr       1.2.0   2017-02-18 CRAN (R 3.4.2)                       
#>  tibble      * 1.3.4   2017-08-22 CRAN (R 3.4.2)                       
#>  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2)                       
#>  tidyverse   * 1.1.1   2017-01-27 CRAN (R 3.4.2)                       
#>  tools         3.4.2   2017-11-01 local                                
#>  treeio      * 1.1.2   2017-11-15 Github (GuangchuangYu/treeio@b6ae142)
#>  tweenr        0.1.5   2016-10-10 cran (@0.1.5)                        
#>  udunits2      0.13    2016-11-17 cran (@0.13)                         
#>  units         0.4-6   2017-08-27 cran (@0.4-6)                        
#>  utils       * 3.4.2   2017-11-01 local                                
#>  vegan       * 2.4-4   2017-08-24 cran (@2.4-4)                        
#>  viridis     * 0.4.0   2017-03-27 CRAN (R 3.4.2)                       
#>  viridisLite * 0.2.0   2017-03-24 CRAN (R 3.4.2)                       
#>  visNetwork  * 2.0.1   2017-07-30 cran (@2.0.1)                        
#>  withr         2.0.0   2017-07-28 CRAN (R 3.4.2)                       
#>  xml2          1.1.1   2017-01-24 CRAN (R 3.4.2)                       
#>  xtable        1.8-2   2016-02-05 CRAN (R 3.4.2)                       
#>  yaml          2.1.14  2016-11-12 CRAN (R 3.4.2)
```
