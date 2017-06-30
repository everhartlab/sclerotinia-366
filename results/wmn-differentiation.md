---
title: "R Notebook"
output: 
  html_notebook:
    toc: true
editor_options: 
  chunk_output_type: inline
---







```r
library("tidyverse")
library("poppr")
(load("data/sclerotinia_16_loci.rda"))
```

```
## [1] "dat"            "dat11"          "datdf"          "keeploci"      
## [5] "corrected_loci"
```

# Purpose

The white mold nursery populations are unique because they are not fungicide
treated and have the same cultivars planted in them year after year.

The question becomes, are white mold nurseries differentiated from each other or
are they more or less homogeneous? We could use AMOVA to test for these with
location and binary source (wmn or non-wmn) as the hierarchy.


# Data Setup

First, we want to clone-correct our data down to the field level so that we
don't accidentally include non-independent samples.


```r
dat11
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    165 original multilocus genotypes 
##    366 haploid individuals
##     11 codominant loci
## 
## Population information:
## 
##      5 strata - MCG, Region, Source, Year, Host
##     14 populations defined - NE, NY, MN, ..., France, Mexico, ND
```

```r
dat11cc <- clonecorrect(dat11, ~Region/Source/Host/Year)
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
##     14 populations defined - NE, NY, MN, ..., France, Mexico, ND
```

Now that we've done that, we should make a new variable in the strata that 
separates the white mold nurseries from the others. We'll call this stratum 
"Source Type".


```r
addStrata(dat11cc) <- strata(dat11cc) %>% 
  mutate(SourceType = forcats::fct_inorder(ifelse(Source == "wmn", "wmn", "other"))) %>%
  select(SourceType)
setPop(dat11cc) <- ~SourceType
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
##      6 strata - MCG, Region, Source, Year, Host, SourceType
##      2 populations defined - other, wmn
```

I can perform AMOVA on the newly defined variable using Bruvo's distance.


```r
other(dat11cc)$REPLEN
```

```
##   5-2(F)   5-3(F)   6-2(F)   7-2(F)   8-3(H)   9-2(F)  12-2(H)  17-3(H) 
##  2.00000  4.00000  5.99999  2.00000  2.00000  2.00000  2.00000  3.00000 
##  20-3(F)  36-4(F)  50-4(F)  55-4(F)  92-4(F) 106-4(H) 110-4(H) 114-4(H) 
##  2.00000  4.00000  4.00000  4.00000  2.00000  4.00000  3.99999  4.00000
```

```r
bd <- bruvo.dist(dat11cc, replen = other(dat11cc)$REPLEN)
(ssc_amova <- poppr.amova(dat11cc, ~SourceType, dist = bd, quiet = TRUE))
```

```
## $call
## ade4::amova(samples = xtab, distances = xdist, structures = xstruct)
## 
## $results
##                  Df     Sum Sq   Mean Sq
## Between samples   1  0.8969901 0.8969901
## Within samples  316 68.7005529 0.2174068
## Total           317 69.5975430 0.2195506
## 
## $componentsofcovariance
##                                   Sigma          %
## Variations  Between samples 0.004391536   1.979968
## Variations  Within samples  0.217406813  98.020032
## Total variations            0.221798349 100.000000
## 
## $statphi
##                          Phi
## Phi-samples-total 0.01979968
```

```r
ssc_amova_test <- randtest(ssc_amova, nrepet = 999)
plot(ssc_amova_test)
```

![plot of chunk AMOVA](./figures/wmn-differentiation///AMOVA-1.png)

```r
ssc_amova_test
```

```
## Monte-Carlo test
## Call: as.randtest(sim = res, obs = sigma[1])
## 
## Observation: 0.004391536 
## 
## Based on 999 replicates
## Simulated p-value: 0.001 
## Alternative hypothesis: greater 
## 
##       Std.Obs   Expectation      Variance 
##  7.949629e+00 -3.847695e-06  3.057032e-07
```


This result is telling us that there is some subdivision between white mold
nurseries and non-white mold nurseries. Of course, from previous analyses, we
know that Mexico is differentiated from other populations, so what happens if we
account for Region? Here, we are placing region lower in the heirarchy because
we specifically want to test the effect of region on the differentiation between
white mold nurseries within different regions.


```r
(ssc_amova_region <- poppr.amova(dat11cc, ~SourceType/Region, dist = bd, quiet = TRUE))
```

```
## $call
## ade4::amova(samples = xtab, distances = xdist, structures = xstruct)
## 
## $results
##                                    Df     Sum Sq   Mean Sq
## Between SourceType                  1  0.8969901 0.8969901
## Between samples Within SourceType  20 12.0338468 0.6016923
## Within samples                    296 56.6667061 0.1914416
## Total                             317 69.5975430 0.2195506
## 
## $componentsofcovariance
##                                                       Sigma           %
## Variations  Between SourceType                -0.0003509858  -0.1583345
## Variations  Between samples Within SourceType  0.0305830685  13.7964379
## Variations  Within samples                     0.1914415747  86.3618966
## Total variations                               0.2216736574 100.0000000
## 
## $statphi
##                                 Phi
## Phi-samples-total       0.136381034
## Phi-samples-SourceType  0.137746279
## Phi-SourceType-total   -0.001583345
```

```r
ssc_amova_region_test <- randtest(ssc_amova_region, nrepet = 999)
plot(ssc_amova_region_test)
```

![plot of chunk AMOVA-Region](./figures/wmn-differentiation///AMOVA-Region-1.png)

```r
ssc_amova_region_test
```

```
## class: krandtest lightkrandtest 
## Monte-Carlo tests
## Call: randtest.amova(xtest = ssc_amova_region, nrepet = 999)
## 
## Number of tests:   3 
## 
## Adjustment method for multiple comparisons:   none 
## Permutation number:   999 
##                            Test           Obs     Std.Obs   Alter Pvalue
## 1     Variations within samples  0.1914415747 -22.5164432    less  0.001
## 2    Variations between samples  0.0305830685  22.0522452 greater  0.001
## 3 Variations between SourceType -0.0003509858   0.4199434 greater  0.252
```

Okay! This shows that when we account for Region after accounting for Source
Type, we find that the differentiation is coming mainly from the Regions. What
happens when we remove Mexico?


```r
datnomex <- setPop(dat11cc, ~Region) %>% popsub(blacklist = "Mexico")
bdnm     <- bruvo.dist(datnomex, replen = other(datnomex)$REPLEN)
(ssc_amova_nm <- poppr.amova(datnomex, ~SourceType/Region, dist = bdnm, quiet = TRUE))
```

```
## $call
## ade4::amova(samples = xtab, distances = xdist, structures = xstruct)
## 
## $results
##                                    Df     Sum Sq   Mean Sq
## Between SourceType                  1  0.7527864 0.7527864
## Between samples Within SourceType  19  8.9146924 0.4691943
## Within samples                    282 55.0854125 0.1953383
## Total                             302 64.7528912 0.2144135
## 
## $componentsofcovariance
##                                                      Sigma          %
## Variations  Between SourceType                0.0002498435   0.115560
## Variations  Between samples Within SourceType 0.0206142904   9.534715
## Variations  Within samples                    0.1953383421  90.349725
## Total variations                              0.2162024759 100.000000
## 
## $statphi
##                               Phi
## Phi-samples-total      0.09650275
## Phi-samples-SourceType 0.09545746
## Phi-SourceType-total   0.00115560
```

```r
ssc_amova_nm_test <- randtest(ssc_amova_nm, nrepet = 999)
plot(ssc_amova_nm_test)
```

![plot of chunk AMOVA-nomex](./figures/wmn-differentiation///AMOVA-nomex-1.png)

```r
ssc_amova_nm_test
```

```
## class: krandtest lightkrandtest 
## Monte-Carlo tests
## Call: randtest.amova(xtest = ssc_amova_nm, nrepet = 999)
## 
## Number of tests:   3 
## 
## Adjustment method for multiple comparisons:   none 
## Permutation number:   999 
##                            Test          Obs     Std.Obs   Alter Pvalue
## 1     Variations within samples 0.1953383421 -14.9295312    less  0.001
## 2    Variations between samples 0.0206142904  13.9977469 greater  0.001
## 3 Variations between SourceType 0.0002498435   0.4447413 greater  0.276
```

When we remove the Mexican isolates (which only contained white mold nurseries
and shared no genotypes), we see that indeed, the degree of differentiation
went down. 

<details>
<summary>Session Information</summary>


```
## Session info --------------------------------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.0 (2017-04-21)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       America/Chicago             
##  date     2017-06-30
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source                                  
##  ade4        * 1.7-6   2017-03-23 CRAN (R 3.4.0)                          
##  adegenet    * 2.1.0   2017-06-30 Github (thibautjombart/adegenet@43b4360)
##  ape           4.1     2017-02-14 CRAN (R 3.4.0)                          
##  assertr       2.0.2.2 2017-06-06 CRAN (R 3.4.0)                          
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0)                          
##  base        * 3.4.0   2017-04-21 local                                   
##  bindr         0.1     2016-11-13 CRAN (R 3.4.0)                          
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.0)                          
##  boot          1.3-19  2017-04-21 CRAN (R 3.4.0)                          
##  broom         0.4.2   2017-02-13 CRAN (R 3.4.0)                          
##  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.0)                          
##  cluster       2.0.6   2017-03-16 CRAN (R 3.4.0)                          
##  coda          0.19-1  2016-12-08 CRAN (R 3.4.0)                          
##  codetools     0.2-15  2016-10-05 CRAN (R 3.4.0)                          
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.0)                          
##  compiler      3.4.0   2017-04-21 local                                   
##  datasets    * 3.4.0   2017-04-21 local                                   
##  DBI           0.7     2017-06-18 CRAN (R 3.4.0)                          
##  deldir        0.1-14  2017-04-22 CRAN (R 3.4.0)                          
##  devtools      1.13.2  2017-06-02 CRAN (R 3.4.0)                          
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0)                          
##  dplyr       * 0.7.1   2017-06-22 CRAN (R 3.4.0)                          
##  evaluate      0.10    2016-10-11 CRAN (R 3.4.0)                          
##  expm          0.999-2 2017-03-29 CRAN (R 3.4.0)                          
##  ezknitr       0.6     2016-09-16 CRAN (R 3.4.0)                          
##  fastmatch     1.1-0   2017-01-28 CRAN (R 3.4.0)                          
##  forcats       0.2.0   2017-01-23 CRAN (R 3.4.0)                          
##  foreign       0.8-69  2017-06-21 CRAN (R 3.4.0)                          
##  gdata         2.18.0  2017-06-06 CRAN (R 3.4.0)                          
##  gdtools     * 0.1.4   2017-03-17 CRAN (R 3.4.0)                          
##  ggforce       0.1.1   2016-11-28 CRAN (R 3.4.0)                          
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.0)                          
##  ggraph      * 1.0.0   2017-02-24 CRAN (R 3.4.0)                          
##  ggrepel     * 0.6.10  2017-06-23 Github (slowkow/ggrepel@102ca39)        
##  glue          1.1.1   2017-06-21 CRAN (R 3.4.0)                          
##  gmodels       2.16.2  2015-07-22 CRAN (R 3.4.0)                          
##  graphics    * 3.4.0   2017-04-21 local                                   
##  grDevices   * 3.4.0   2017-04-21 local                                   
##  grid          3.4.0   2017-04-21 local                                   
##  gridExtra     2.2.1   2016-02-29 CRAN (R 3.4.0)                          
##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.0)                          
##  gtools        3.5.0   2015-05-29 CRAN (R 3.4.0)                          
##  haven         1.0.0   2016-09-23 CRAN (R 3.4.0)                          
##  highr         0.6     2016-05-09 CRAN (R 3.4.0)                          
##  hms           0.3     2016-11-22 CRAN (R 3.4.0)                          
##  htmltools     0.3.6   2017-04-28 CRAN (R 3.4.0)                          
##  httpuv        1.3.3   2015-08-04 CRAN (R 3.4.0)                          
##  httr          1.2.1   2016-07-03 CRAN (R 3.4.0)                          
##  igraph      * 1.0.1   2015-06-26 CRAN (R 3.4.0)                          
##  jsonlite      1.5     2017-06-01 CRAN (R 3.4.0)                          
##  knitr       * 1.16    2017-05-18 CRAN (R 3.4.0)                          
##  lattice     * 0.20-35 2017-03-25 CRAN (R 3.4.0)                          
##  lazyeval      0.2.0   2016-06-12 CRAN (R 3.4.0)                          
##  LearnBayes    2.15    2014-05-29 CRAN (R 3.4.0)                          
##  lubridate     1.6.0   2016-09-13 CRAN (R 3.4.0)                          
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0)                          
##  MASS          7.3-47  2017-04-21 CRAN (R 3.4.0)                          
##  Matrix        1.2-10  2017-04-28 CRAN (R 3.4.0)                          
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0)                          
##  methods     * 3.4.0   2017-04-21 local                                   
##  mgcv          1.8-17  2017-02-08 CRAN (R 3.4.0)                          
##  mime          0.5     2016-07-07 CRAN (R 3.4.0)                          
##  mnormt        1.5-5   2016-10-15 CRAN (R 3.4.0)                          
##  modelr        0.1.0   2016-08-31 CRAN (R 3.4.0)                          
##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.0)                          
##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.0)                          
##  parallel      3.4.0   2017-04-21 local                                   
##  pegas         0.10    2017-05-03 CRAN (R 3.4.0)                          
##  permute     * 0.9-4   2016-09-09 CRAN (R 3.4.0)                          
##  phangorn      2.2.0   2017-04-03 CRAN (R 3.4.0)                          
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0)                          
##  plyr          1.8.4   2016-06-08 CRAN (R 3.4.0)                          
##  poppr       * 2.4.1   2017-04-14 CRAN (R 3.4.0)                          
##  psych         1.7.5   2017-05-03 CRAN (R 3.4.0)                          
##  purrr       * 0.2.2.2 2017-05-11 cran (@0.2.2.2)                         
##  quadprog      1.5-5   2013-04-17 CRAN (R 3.4.0)                          
##  R.methodsS3   1.7.1   2016-02-16 CRAN (R 3.4.0)                          
##  R.oo          1.21.0  2016-11-01 CRAN (R 3.4.0)                          
##  R.utils       2.5.0   2016-11-07 CRAN (R 3.4.0)                          
##  R6            2.2.2   2017-06-17 cran (@2.2.2)                           
##  Rcpp          0.12.11 2017-05-22 cran (@0.12.11)                         
##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.0)                          
##  readxl        1.0.0   2017-04-18 CRAN (R 3.4.0)                          
##  reshape2      1.4.2   2016-10-22 CRAN (R 3.4.0)                          
##  rlang         0.1.1   2017-05-18 CRAN (R 3.4.0)                          
##  rvest         0.3.2   2016-06-17 CRAN (R 3.4.0)                          
##  scales        0.4.1   2016-11-09 CRAN (R 3.4.0)                          
##  seqinr        3.3-6   2017-04-06 CRAN (R 3.4.0)                          
##  shiny         1.0.3   2017-04-26 CRAN (R 3.4.0)                          
##  sp            1.2-4   2016-12-22 CRAN (R 3.4.0)                          
##  spdep         0.6-13  2017-04-25 CRAN (R 3.4.0)                          
##  splines       3.4.0   2017-04-21 local                                   
##  stats       * 3.4.0   2017-04-21 local                                   
##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.0)                          
##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.0)                          
##  tibble      * 1.3.3   2017-05-28 CRAN (R 3.4.0)                          
##  tidyr       * 0.6.3   2017-05-15 CRAN (R 3.4.0)                          
##  tidyverse   * 1.1.1   2017-01-27 CRAN (R 3.4.0)                          
##  tools         3.4.0   2017-04-21 local                                   
##  tweenr        0.1.5   2016-10-10 CRAN (R 3.4.0)                          
##  udunits2      0.13    2016-11-17 CRAN (R 3.4.0)                          
##  units         0.4-5   2017-06-15 CRAN (R 3.4.0)                          
##  utils       * 3.4.0   2017-04-21 local                                   
##  vegan       * 2.4-3   2017-04-07 CRAN (R 3.4.0)                          
##  viridis     * 0.4.0   2017-03-27 CRAN (R 3.4.0)                          
##  viridisLite * 0.2.0   2017-03-24 CRAN (R 3.4.0)                          
##  withr         1.0.2   2016-06-20 CRAN (R 3.4.0)                          
##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.0)                          
##  xtable        1.8-2   2016-02-05 CRAN (R 3.4.0)
```

</details>
