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
##  version  R version 3.4.1 (2017-06-30)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       America/Chicago             
##  date     2017-08-17
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version    date       source                        
##  ade4        * 1.7-6      2017-03-23 CRAN (R 3.4.0)                
##  adegenet    * 2.1.0      2017-07-17 local                         
##  ape           4.1        2017-02-14 CRAN (R 3.4.0)                
##  assertr       2.0.2.2    2017-06-06 CRAN (R 3.4.0)                
##  assertthat    0.2.0      2017-04-11 CRAN (R 3.4.0)                
##  base        * 3.4.1      2017-07-07 local                         
##  bindr         0.1        2016-11-13 CRAN (R 3.4.0)                
##  bindrcpp    * 0.2        2017-06-17 CRAN (R 3.4.0)                
##  boot          1.3-20     2017-07-30 CRAN (R 3.4.1)                
##  broom         0.4.2      2017-02-13 CRAN (R 3.4.0)                
##  cellranger    1.1.0      2016-07-27 CRAN (R 3.4.0)                
##  cluster       2.0.6      2017-03-16 CRAN (R 3.4.0)                
##  coda          0.19-1     2016-12-08 CRAN (R 3.4.0)                
##  colorspace    1.3-2      2016-12-14 CRAN (R 3.4.0)                
##  compiler      3.4.1      2017-07-07 local                         
##  datasets    * 3.4.1      2017-07-07 local                         
##  deldir        0.1-14     2017-04-22 CRAN (R 3.4.0)                
##  devtools      1.13.3     2017-08-02 CRAN (R 3.4.1)                
##  digest        0.6.12     2017-01-27 CRAN (R 3.4.0)                
##  dplyr       * 0.7.2      2017-07-20 CRAN (R 3.4.1)                
##  evaluate      0.10.1     2017-06-24 CRAN (R 3.4.1)                
##  expm          0.999-2    2017-03-29 CRAN (R 3.4.0)                
##  ezknitr       0.6        2016-09-16 CRAN (R 3.4.0)                
##  fastmatch     1.1-0      2017-01-28 CRAN (R 3.4.0)                
##  forcats       0.2.0      2017-01-23 CRAN (R 3.4.0)                
##  foreign       0.8-69     2017-06-21 CRAN (R 3.4.0)                
##  gdata         2.18.0     2017-06-06 CRAN (R 3.4.0)                
##  ggplot2     * 2.2.1      2016-12-30 CRAN (R 3.4.0)                
##  ggtree      * 1.6.11     2017-08-03 Bioconductor                  
##  glue          1.1.1      2017-06-21 CRAN (R 3.4.0)                
##  gmodels       2.16.2     2015-07-22 CRAN (R 3.4.0)                
##  graphics    * 3.4.1      2017-07-07 local                         
##  grDevices   * 3.4.1      2017-07-07 local                         
##  grid          3.4.1      2017-07-07 local                         
##  gridExtra     2.2.1      2016-02-29 CRAN (R 3.4.0)                
##  gtable        0.2.0      2016-02-26 CRAN (R 3.4.0)                
##  gtools        3.5.0      2015-05-29 CRAN (R 3.4.0)                
##  haven         1.1.0      2017-07-09 CRAN (R 3.4.1)                
##  highr         0.6        2016-05-09 CRAN (R 3.4.0)                
##  hms           0.3        2016-11-22 CRAN (R 3.4.0)                
##  htmltools     0.3.6      2017-04-28 CRAN (R 3.4.0)                
##  httpuv        1.3.5      2017-07-04 CRAN (R 3.4.1)                
##  httr          1.2.1      2016-07-03 CRAN (R 3.4.0)                
##  igraph        1.1.2      2017-07-21 cran (@1.1.2)                 
##  jsonlite      1.5        2017-06-01 CRAN (R 3.4.0)                
##  knitr       * 1.16       2017-05-18 CRAN (R 3.4.0)                
##  labeling      0.3        2014-08-23 CRAN (R 3.4.0)                
##  lattice       0.20-35    2017-03-25 CRAN (R 3.4.0)                
##  lazyeval      0.2.0      2016-06-12 CRAN (R 3.4.0)                
##  LearnBayes    2.15       2014-05-29 CRAN (R 3.4.0)                
##  lubridate     1.6.0      2016-09-13 CRAN (R 3.4.0)                
##  magrittr      1.5        2014-11-22 CRAN (R 3.4.0)                
##  MASS          7.3-47     2017-04-21 CRAN (R 3.4.0)                
##  Matrix        1.2-10     2017-04-28 CRAN (R 3.4.0)                
##  memoise       1.1.0      2017-04-21 CRAN (R 3.4.0)                
##  methods     * 3.4.1      2017-07-07 local                         
##  mgcv          1.8-18     2017-07-28 CRAN (R 3.4.1)                
##  mime          0.5        2016-07-07 CRAN (R 3.4.0)                
##  mnormt        1.5-5      2016-10-15 CRAN (R 3.4.0)                
##  modelr        0.1.1      2017-07-24 CRAN (R 3.4.1)                
##  munsell       0.4.3      2016-02-13 CRAN (R 3.4.0)                
##  nlme          3.1-131    2017-02-06 CRAN (R 3.4.0)                
##  parallel      3.4.1      2017-07-07 local                         
##  pegas         0.10       2017-05-03 CRAN (R 3.4.0)                
##  permute       0.9-4      2016-09-09 CRAN (R 3.4.0)                
##  phangorn      2.2.0      2017-04-03 CRAN (R 3.4.0)                
##  pkgconfig     2.0.1      2017-03-21 CRAN (R 3.4.0)                
##  plyr          1.8.4      2016-06-08 CRAN (R 3.4.0)                
##  poppr       * 2.4.1.99-2 2017-08-13 local                         
##  psych         1.7.5      2017-05-03 CRAN (R 3.4.0)                
##  purrr       * 0.2.3      2017-08-02 CRAN (R 3.4.1)                
##  quadprog      1.5-5      2013-04-17 CRAN (R 3.4.0)                
##  R.methodsS3   1.7.1      2016-02-16 CRAN (R 3.4.0)                
##  R.oo          1.21.0     2016-11-01 CRAN (R 3.4.0)                
##  R.utils       2.5.0      2016-11-07 CRAN (R 3.4.0)                
##  R6            2.2.2      2017-06-17 cran (@2.2.2)                 
##  Rcpp          0.12.12    2017-07-15 cran (@0.12.12)               
##  readr       * 1.1.1      2017-05-16 CRAN (R 3.4.0)                
##  readxl        1.0.0      2017-04-18 CRAN (R 3.4.0)                
##  reshape2      1.4.2      2016-10-22 CRAN (R 3.4.0)                
##  rlang         0.1.1      2017-05-18 CRAN (R 3.4.0)                
##  rvest         0.3.2      2016-06-17 CRAN (R 3.4.0)                
##  scales        0.4.1.9002 2017-08-02 Github (hadley/scales@842ad87)
##  seqinr        3.4-5      2017-08-01 CRAN (R 3.4.1)                
##  shiny         1.0.3      2017-04-26 CRAN (R 3.4.0)                
##  sp            1.2-5      2017-06-29 CRAN (R 3.4.1)                
##  spdep         0.6-13     2017-04-25 CRAN (R 3.4.0)                
##  splines       3.4.1      2017-07-07 local                         
##  stats       * 3.4.1      2017-07-07 local                         
##  stats4        3.4.1      2017-07-07 local                         
##  stringi       1.1.5      2017-04-07 CRAN (R 3.4.0)                
##  stringr       1.2.0      2017-02-18 CRAN (R 3.4.0)                
##  tibble      * 1.3.3      2017-05-28 CRAN (R 3.4.0)                
##  tidyr       * 0.6.3      2017-05-15 CRAN (R 3.4.0)                
##  tidyverse   * 1.1.1      2017-01-27 CRAN (R 3.4.0)                
##  tools         3.4.1      2017-07-07 local                         
##  utils       * 3.4.1      2017-07-07 local                         
##  vegan         2.4-3      2017-04-07 CRAN (R 3.4.0)                
##  viridis       0.4.0      2017-03-27 CRAN (R 3.4.0)                
##  viridisLite   0.2.0      2017-03-24 CRAN (R 3.4.0)                
##  withr         2.0.0      2017-07-28 CRAN (R 3.4.1)                
##  xml2          1.1.1      2017-01-24 CRAN (R 3.4.0)                
##  xtable        1.8-2      2016-02-05 CRAN (R 3.4.0)
```

</details>
