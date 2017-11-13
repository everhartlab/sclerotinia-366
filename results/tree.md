---
title: "tree"
---






## Packages and Data


```r
library('tidyverse')
library('poppr')
library('ggtree')
```


```r
load(file.path(PROJHOME, "data", "sclerotinia_16_loci.rda"))
setPop(dat11) <- ~Host/Source/Region/Year
dat11cc <- clonecorrect(dat11, ~Host/Source/Region/Year, keep = 1:4)
dat11cc
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    165 original multilocus genotypes 
##    318 haploid individuals
##     11 codominant loci
## 
## Population information:
## 
##      5 strata - MCG, Region, Source, Year, Host
##    128 populations defined - 
## GH_unk_NE_2003, GH_unk_NY_2003, G122_wmn_MN_2003, ..., unk_pmc_ND_2010, unk_wlc_ND_2010, unk_flds_France_2012
```

```r
# Asserting that nothing messed up with the metadata.
stopifnot(identical(indNames(dat11cc), other(dat11cc)$meta$Isolate))
```


The purpose of this document is simply to calculate a bootstrapped tree for 
Bruvo's distance.


```r
set.seed(2017-08-03)
bt <- bruvo.boot(clonecorrect(dat11, strata = NA), 
                 replen = other(dat11)$REPLEN, 
                 sample = 1000, 
                 tree = "nj", 
                 showtree = FALSE)
```

```
## 
## Bootstrapping...
## (note: calculation of node labels can take a while even after the progress bar is full)
## 
## Running bootstraps:       100 / 1000Running bootstraps:       200 / 1000Running bootstraps:       300 / 1000Running bootstraps:       400 / 1000Running bootstraps:       500 / 1000Running bootstraps:       600 / 1000Running bootstraps:       700 / 1000Running bootstraps:       800 / 1000Running bootstraps:       900 / 1000Running bootstraps:       1000 / 1000
## Calculating bootstrap values... done.
```

Let's take a look at the results of the bootstrap analysis. Note, that we would
traditionally ignore results < 75.


```r
bt$node.labels
```

```
##   [1] 100   3   3   0   0   0   0   0   0   6  17   1   0  48   0  16   0
##  [18]  26   0   0   0   0   0   0  15   2   1  37  10  19   0   0  25   0
##  [35]  21  24   0   0  34   0   1   0   0   0   9   0   0   7  17   1   2
##  [52]  21   0   1   0   0  10  22   0  26   0   1   0  18   0   8  32   9
##  [69]   0   3  27   4  12   0   5  24  51  28   0  11   1   4   0   1  14
##  [86]   0  29  32  54  26   5   2   6   2  51  27   3  33   2  25   7  17
## [103]  30   2  32  25  25  16  10  26  24  13  25  24  19  30  57  36  37
## [120]  46  47  60   4  46  40  30  37  24  20  44   1  60  31  29  55   1
## [137]  35  61  42   1  51  70  24   3  12   6  58  45  40  35   9  34  21
## [154]  19  37   6  40  46  46  41  22  47  43
```

```r
summary(bt$node.labels)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     0.0     1.0    13.0    18.4    30.5   100.0
```

```r
hist(bt$node.labels)
```

<img src="./figures/tree///bsresults-1.png" title="plot of chunk bsresults" alt="plot of chunk bsresults" style="display: block; margin: auto;" />

It's clear that our results show that NONE of the clades are well supported
(except for the whole tree, which is default), so let's take a look at how the
tree looks. One of the things we want to do is see where the popualtions fit.
Since we created a clone-corrected tree here (for speed), we are going to create
a matrix to tally up the samples from different populations per MLG.

An important thing to keep track of is the fact that ggtree works from rownames
for tip labels, so we have to add them in. 


```r
data <- bind_cols(strata(dat11), other(dat11)$meta) %>% 
  add_column(MLG = mll(dat11))

otherdf <- data %>% 
  group_by(MLG, Region) %>%
  summarize(N = n()) %>%
  spread(Region, N)
otherdf <- inner_join(data %>% 
                        group_by(MLG) %>% 
                        summarize(N = n(), id = Isolate[1]), 
                      otherdf, 
                      by = "MLG")
df <- otherdf %>% 
  select(-MLG, -N) %>% 
  as.data.frame() %>% 
  column_to_rownames("id")
otherdf
```

```
## # A tibble: 165 x 17
##      MLG     N    id    NE    NY    MN    MI    OR    WA    CO    WI    ID
##    <int> <int> <chr> <int> <int> <int> <int> <int> <int> <int> <int> <int>
##  1     1     2   586    NA    NA    NA    NA     1    NA    NA    NA    NA
##  2     2     2   501    NA    NA     1    NA    NA    NA    NA    NA    NA
##  3     3     2   502    NA    NA     1    NA    NA    NA    NA    NA    NA
##  4     4     3   500    NA    NA     1    NA    NA    NA    NA    NA    NA
##  5     5     2   585    NA    NA    NA    NA     1    NA    NA    NA    NA
##  6     6     2   596    NA    NA    NA    NA    NA    NA    NA    NA    NA
##  7     7     2   714     1    NA    NA    NA    NA    NA    NA    NA    NA
##  8     8     1   456    NA    NA    NA    NA    NA     1    NA    NA    NA
##  9     9     4   703     1    NA    NA     1    NA    NA    NA    NA    NA
## 10    10     1   482    NA    NA    NA     1    NA    NA    NA    NA    NA
## # ... with 155 more rows, and 5 more variables: Australia <int>, CA <int>,
## #   France <int>, Mexico <int>, ND <int>
```

Now we can create the tree with the matrix, coloring by number of isolates.


```r
gbt <- ggtree(bt) + geom_tippoint() + theme_tree2() + xlab("Bruvo's Distance (11 loci)")
gheatmap(gbt, df, colnames = FALSE, width = 0.6) %>% 
  scale_x_ggtree() + 
  viridis::scale_fill_viridis(guide = "legend") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")) +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = 3)) +
  labs(list(fill = "N isolates"))
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill',
## which will replace the existing scale.
```

<img src="./figures/tree///the_tree-1.png" title="plot of chunk the_tree" alt="plot of chunk the_tree" style="display: block; margin: auto;" />

Of course, because of the low bootstrap support, we don't have very much faith
in this tree. One of the promising aspects, however is the fact that all the
isolates from Mexico group together, as was shown in the DAPC. Moreover, we can
see that this is not driven by private alleles:


```r
private_alleles(dat11, locus ~ Region, count.alleles = FALSE)
```

```
##           9-2(F) 12-2(H) 55-4(F)
## NE             0       0       0
## NY             0       0       0
## MN             0       0       0
## MI             0       0       0
## OR             0       0       0
## WA             0       0       1
## CO             0       0       0
## WI             0       0       0
## ID             0       0       0
## Australia      0       0       0
## CA             0       0       0
## France         0       0       0
## Mexico         0       0       0
## ND             1       1       0
```




<details>
<summary>Session Information</summary>


```
## Session info --------------------------------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.2 (2017-09-28)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       UTC                         
##  date     2017-11-13
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source                               
##  ade4        * 1.7-8   2017-08-09 cran (@1.7-8)                        
##  adegenet    * 2.1.0   2017-10-12 cran (@2.1.0)                        
##  ape           5.0     2017-10-30 cran (@5.0)                          
##  assertr       2.0.2.2 2017-06-06 cran (@2.0.2.2)                      
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.2)                       
##  base        * 3.4.2   2017-11-01 local                                
##  bindr         0.1     2016-11-13 CRAN (R 3.4.2)                       
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.2)                       
##  boot          1.3-20  2017-07-30 cran (@1.3-20)                       
##  broom         0.4.2   2017-02-13 CRAN (R 3.4.2)                       
##  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.2)                       
##  cluster       2.0.6   2017-03-16 CRAN (R 3.4.2)                       
##  coda          0.19-1  2016-12-08 cran (@0.19-1)                       
##  codetools     0.2-15  2016-10-05 CRAN (R 3.4.2)                       
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.2)                       
##  compiler      3.4.2   2017-11-01 local                                
##  datasets    * 3.4.2   2017-11-01 local                                
##  deldir        0.1-14  2017-04-22 cran (@0.1-14)                       
##  devtools      1.13.3  2017-08-02 CRAN (R 3.4.2)                       
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.2)                       
##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)                       
##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.2)                       
##  expm          0.999-2 2017-03-29 cran (@0.999-2)                      
##  ezknitr       0.6     2016-09-16 cran (@0.6)                          
##  fastmatch     1.1-0   2017-01-28 cran (@1.1-0)                        
##  forcats       0.2.0   2017-01-23 CRAN (R 3.4.2)                       
##  foreign       0.8-69  2017-06-21 CRAN (R 3.4.2)                       
##  gdata         2.18.0  2017-06-06 cran (@2.18.0)                       
##  ggcompoplot * 0.1.0   2017-11-09 Github (zkamvar/ggcompoplot@bcf007d) 
##  ggforce       0.1.1   2016-11-28 cran (@0.1.1)                        
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.2)                       
##  ggraph      * 1.0.0   2017-02-24 cran (@1.0.0)                        
##  ggrepel       0.7.0   2017-09-29 cran (@0.7.0)                        
##  ggtree      * 1.9.4   2017-11-09 Github (GuangchuangYu/ggtree@07063f9)
##  glue          1.2.0   2017-10-29 CRAN (R 3.4.2)                       
##  gmodels       2.16.2  2015-07-22 cran (@2.16.2)                       
##  graphics    * 3.4.2   2017-11-01 local                                
##  grDevices   * 3.4.2   2017-11-01 local                                
##  grid          3.4.2   2017-11-01 local                                
##  gridExtra     2.3     2017-09-09 CRAN (R 3.4.2)                       
##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.2)                       
##  gtools        3.5.0   2015-05-29 cran (@3.5.0)                        
##  haven         1.1.0   2017-07-09 CRAN (R 3.4.2)                       
##  highr         0.6     2016-05-09 CRAN (R 3.4.2)                       
##  hms           0.3     2016-11-22 CRAN (R 3.4.2)                       
##  htmltools     0.3.6   2017-04-28 CRAN (R 3.4.2)                       
##  httpuv        1.3.5   2017-07-04 CRAN (R 3.4.2)                       
##  httr          1.3.1   2017-08-20 CRAN (R 3.4.2)                       
##  igraph      * 1.1.2   2017-07-21 CRAN (R 3.4.2)                       
##  jsonlite      1.5     2017-06-01 CRAN (R 3.4.2)                       
##  knitr       * 1.17    2017-08-10 CRAN (R 3.4.2)                       
##  labeling      0.3     2014-08-23 CRAN (R 3.4.2)                       
##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.2)                       
##  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2)                       
##  LearnBayes    2.15    2014-05-29 cran (@2.15)                         
##  lubridate     1.7.0   2017-10-29 CRAN (R 3.4.2)                       
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.2)                       
##  MASS          7.3-47  2017-04-21 CRAN (R 3.4.2)                       
##  Matrix        1.2-11  2017-08-16 CRAN (R 3.4.2)                       
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.2)                       
##  methods     * 3.4.2   2017-11-01 local                                
##  mgcv          1.8-22  2017-09-19 CRAN (R 3.4.2)                       
##  mime          0.5     2016-07-07 CRAN (R 3.4.2)                       
##  mnormt        1.5-5   2016-10-15 CRAN (R 3.4.2)                       
##  modelr        0.1.1   2017-07-24 CRAN (R 3.4.2)                       
##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.2)                       
##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.2)                       
##  parallel      3.4.2   2017-11-01 local                                
##  pegas         0.10    2017-05-03 cran (@0.10)                         
##  permute       0.9-4   2016-09-09 cran (@0.9-4)                        
##  phangorn      2.2.0   2017-04-03 cran (@2.2.0)                        
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.2)                       
##  plyr          1.8.4   2016-06-08 CRAN (R 3.4.2)                       
##  poppr       * 2.5.0   2017-09-11 cran (@2.5.0)                        
##  psych         1.7.8   2017-09-09 CRAN (R 3.4.2)                       
##  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2)                       
##  quadprog      1.5-5   2013-04-17 cran (@1.5-5)                        
##  R.methodsS3   1.7.1   2016-02-16 cran (@1.7.1)                        
##  R.oo          1.21.0  2016-11-01 cran (@1.21.0)                       
##  R.utils       2.5.0   2016-11-07 cran (@2.5.0)                        
##  R6            2.2.2   2017-06-17 CRAN (R 3.4.2)                       
##  Rcpp          0.12.13 2017-09-28 CRAN (R 3.4.2)                       
##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.2)                       
##  readxl        1.0.0   2017-04-18 CRAN (R 3.4.2)                       
##  reshape2      1.4.2   2016-10-22 CRAN (R 3.4.2)                       
##  rlang         0.1.2   2017-08-09 CRAN (R 3.4.2)                       
##  rvcheck       0.0.9   2017-07-10 cran (@0.0.9)                        
##  rvest         0.3.2   2016-06-17 CRAN (R 3.4.2)                       
##  scales        0.5.0   2017-08-24 CRAN (R 3.4.2)                       
##  seqinr        3.4-5   2017-08-01 cran (@3.4-5)                        
##  shiny         1.0.5   2017-08-23 CRAN (R 3.4.2)                       
##  sp            1.2-5   2017-06-29 CRAN (R 3.4.2)                       
##  spdep         0.6-15  2017-09-01 cran (@0.6-15)                       
##  splines       3.4.2   2017-11-01 local                                
##  stats       * 3.4.2   2017-11-01 local                                
##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.2)                       
##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.2)                       
##  tibble      * 1.3.4   2017-08-22 CRAN (R 3.4.2)                       
##  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2)                       
##  tidyselect    0.2.2   2017-10-10 CRAN (R 3.4.2)                       
##  tidyverse   * 1.1.1   2017-01-27 CRAN (R 3.4.2)                       
##  tools         3.4.2   2017-11-01 local                                
##  treeio      * 1.1.2   2017-11-09 Github (GuangchuangYu/treeio@b6ae142)
##  tweenr        0.1.5   2016-10-10 cran (@0.1.5)                        
##  udunits2      0.13    2016-11-17 cran (@0.13)                         
##  units         0.4-6   2017-08-27 cran (@0.4-6)                        
##  utils       * 3.4.2   2017-11-01 local                                
##  vegan         2.4-4   2017-08-24 cran (@2.4-4)                        
##  viridis       0.4.0   2017-03-27 CRAN (R 3.4.2)                       
##  viridisLite   0.2.0   2017-03-24 CRAN (R 3.4.2)                       
##  withr         2.0.0   2017-07-28 CRAN (R 3.4.2)                       
##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.2)                       
##  xtable        1.8-2   2016-02-05 CRAN (R 3.4.2)
```

</details>
