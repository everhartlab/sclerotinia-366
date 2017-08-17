---
title: "Locus Statistics"
output: 
  html_notebook:
    toc: true
editor_options: 
  chunk_output_type: inline
---



# Purpose

This will tabulate statistics per locus for presentation in a table.


```r
library("tidyverse")
library("poppr")
load(file.path(PROJHOME, "data", "sclerotinia_16_loci.rda"))
```


```r
makerange <- . %>% as.integer() %>% range() %>% paste(collapse = "-")

ranges <- map_chr(alleles(dat11), makerange)
ranges <- c(ranges, alleles(dat11) %>% unlist(use.names = FALSE) %>% makerange)
locus_table(dat11, information = FALSE) %>% 
  as.data.frame() %>%
  rownames_to_column("Locus") %>%
  mutate(Locus = gsub("\\([FH]\\)", "", Locus)) %>%
  add_column(Range = ranges) %>%
  select(Locus, Range, allele, Hexp, Evenness) %>%
  mutate(allele = round(allele, 2)) %>%
  rename(`*h*` = Hexp, `No. alleles` = allele) %>%
  huxtable::as_huxtable(add_colnames = TRUE) %>%
  huxtable::set_number_format(1:12, 3, 0) %>%
  huxtable::set_number_format(huxtable::everywhere, 4:5, 2) %>%
  huxtable::set_align(huxtable::everywhere, 3:5, "right") %>%
  huxtable::set_col_width(c(0.06, 0.07, 0.12, 0.06, 0.09)) %>%
  huxtable::print_md(max_width = 47)
```

```
--------------------------------------------
Locus  Range    No. alleles    *h*  Evenness 
------ ------- ------------ ------ ---------
5-2    318-324            4   0.45      0.62 

6-2    483-495            3   0.64      0.95 

7-2    158-174            7   0.73      0.76 

8-3    244-270            7   0.74      0.79 

9-2    360-382            9   0.35      0.41 

12-2   214-222            5   0.58      0.78 

17-3   342-363            7   0.55      0.53 

20-3   280-282            2   0.05      0.42 

55-4   153-216           10   0.72      0.66 

110-4  370-386            5   0.76      0.91 

114-4  339-416           10   0.83      0.80 

mean   153-495         6.27   0.58      0.69 

--------------------------------------------
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
##  glue          1.1.1      2017-06-21 CRAN (R 3.4.0)                
##  gmodels       2.16.2     2015-07-22 CRAN (R 3.4.0)                
##  graphics    * 3.4.1      2017-07-07 local                         
##  grDevices   * 3.4.1      2017-07-07 local                         
##  grid          3.4.1      2017-07-07 local                         
##  gtable        0.2.0      2016-02-26 CRAN (R 3.4.0)                
##  gtools        3.5.0      2015-05-29 CRAN (R 3.4.0)                
##  haven         1.1.0      2017-07-09 CRAN (R 3.4.1)                
##  hms           0.3        2016-11-22 CRAN (R 3.4.0)                
##  htmltools     0.3.6      2017-04-28 CRAN (R 3.4.0)                
##  httpuv        1.3.5      2017-07-04 CRAN (R 3.4.1)                
##  httr          1.2.1      2016-07-03 CRAN (R 3.4.0)                
##  huxtable      0.3.0      2017-05-18 CRAN (R 3.4.0)                
##  igraph        1.1.2      2017-07-21 cran (@1.1.2)                 
##  jsonlite      1.5        2017-06-01 CRAN (R 3.4.0)                
##  knitr       * 1.16       2017-05-18 CRAN (R 3.4.0)                
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
##  stringi       1.1.5      2017-04-07 CRAN (R 3.4.0)                
##  stringr       1.2.0      2017-02-18 CRAN (R 3.4.0)                
##  tibble      * 1.3.3      2017-05-28 CRAN (R 3.4.0)                
##  tidyr       * 0.6.3      2017-05-15 CRAN (R 3.4.0)                
##  tidyverse   * 1.1.1      2017-01-27 CRAN (R 3.4.0)                
##  tools         3.4.1      2017-07-07 local                         
##  utils       * 3.4.1      2017-07-07 local                         
##  vegan         2.4-3      2017-04-07 CRAN (R 3.4.0)                
##  withr         2.0.0      2017-07-28 CRAN (R 3.4.1)                
##  xml2          1.1.1      2017-01-24 CRAN (R 3.4.0)                
##  xtable        1.8-2      2016-02-05 CRAN (R 3.4.0)
```

</details>
