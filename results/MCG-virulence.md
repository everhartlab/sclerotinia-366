---
title: "Assessing Virulence by MCG"
output: 
  html_notebook:
    toc: true
---




```r
library("tidyverse")
```

```
## Loading tidyverse: ggplot2
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
## Loading tidyverse: dplyr
```

```
## Conflicts with tidy packages ----------------------------------------------
```

```
## filter(): dplyr, stats
## lag():    dplyr, stats
```

# Purpose

This document is to assess virulence associated with the 10 most common MCGs. We
can get this from the cleaned data set that we saved in the
`data-comparison.Rmd` file.


```r
cols(
  .default = col_integer(),
  Severity = col_double(),
  Region = col_character(),
  Source = col_character(),
  Host = col_character()
)
```

```
## cols(
##   .default = col_integer(),
##   Severity = col_double(),
##   Region = col_character(),
##   Source = col_character(),
##   Host = col_character()
## )
```

```r
dat <- read_csv(file.path(PROJHOME, "data", "clean_data.csv")) %>%
  select(MCG, Severity)
```

```
## Parsed with column specification:
## cols(
##   .default = col_integer(),
##   Severity = col_double(),
##   Region = col_character(),
##   Source = col_character(),
##   Host = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

```r
dat
```

```
## # A tibble: 366 x 2
##      MCG Severity
##    <int>    <dbl>
##  1     4      3.9
##  2    45      5.4
##  3     5      6.3
##  4     4      4.4
##  5     4      4.7
##  6     3      6.1
##  7     5      5.5
##  8     3      5.0
##  9     3      5.2
## 10     5      5.3
## # ... with 356 more rows
```

Now, I can filter out the top 10 MCGs:


```r
top_mcg <- dat %>% 
  group_by(MCG) %>%
  summarize(N = n()) %>%
  arrange(desc(N)) %>%
  slice(1:10) %>% 
  inner_join(dat, by = "MCG") %>%
  select(-N) %>%
  mutate(MCG = forcats::fct_inorder(as.character(MCG)))
top_mcg
```

```
## # A tibble: 207 x 2
##       MCG Severity
##    <fctr>    <dbl>
##  1      5      6.3
##  2      5      5.5
##  3      5      5.3
##  4      5      6.0
##  5      5      5.2
##  6      5      5.8
##  7      5      6.2
##  8      5      5.3
##  9      5      4.3
## 10      5      4.5
## # ... with 197 more rows
```

## Visualizing distributions


```r
ggplot(top_mcg, aes(x = MCG, y = Severity)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.5) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(aspect.ratio = 1/2) +
  labs(list(
    title = "Severity for the top 10 MCGs",
    caption = "MCG = Mycelial Compatibility Group"
    ))
```

![plot of chunk vis](./figures/MCG-virulence///vis-1.png)

## AMOVA

The default AMOVA in R sets contrasts as `contrast.treatment`, which compares
everything to the first factor, considered the treatment. Since we are
interested in whether or not there IS a difference between samples, this will
be sufficient.


```r
AMOVA <- aov(Severity ~ MCG, data = top_mcg)
AMOVA
```

```
## Call:
##    aov(formula = Severity ~ MCG, data = top_mcg)
## 
## Terms:
##                       MCG Residuals
## Sum of Squares   36.67876 163.83080
## Deg. of Freedom         9       197
## 
## Residual standard error: 0.9119366
## Estimated effects may be unbalanced
```

```r
summary(AMOVA)
```

```
##              Df Sum Sq Mean Sq F value   Pr(>F)    
## MCG           9  36.68   4.075   4.901 6.19e-06 ***
## Residuals   197 163.83   0.832                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

There appears to be a significant effect at p = 6.189e-06.


## Session Information


```r
options(width = 100)
devtools::session_info()
```

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
##  date     2017-06-29
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source         
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0) 
##  base        * 3.4.0   2017-04-21 local          
##  bindr         0.1     2016-11-13 CRAN (R 3.4.0) 
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.0) 
##  broom         0.4.2   2017-02-13 CRAN (R 3.4.0) 
##  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.0) 
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.0) 
##  compiler      3.4.0   2017-04-21 local          
##  datasets    * 3.4.0   2017-04-21 local          
##  devtools      1.13.2  2017-06-02 CRAN (R 3.4.0) 
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0) 
##  dplyr       * 0.7.1   2017-06-22 CRAN (R 3.4.0) 
##  evaluate      0.10    2016-10-11 CRAN (R 3.4.0) 
##  ezknitr       0.6     2016-09-16 CRAN (R 3.4.0) 
##  forcats       0.2.0   2017-01-23 CRAN (R 3.4.0) 
##  foreign       0.8-69  2017-06-21 CRAN (R 3.4.0) 
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.0) 
##  glue          1.1.1   2017-06-21 CRAN (R 3.4.0) 
##  graphics    * 3.4.0   2017-04-21 local          
##  grDevices   * 3.4.0   2017-04-21 local          
##  grid          3.4.0   2017-04-21 local          
##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.0) 
##  haven         1.0.0   2016-09-23 CRAN (R 3.4.0) 
##  highr         0.6     2016-05-09 CRAN (R 3.4.0) 
##  hms           0.3     2016-11-22 CRAN (R 3.4.0) 
##  httr          1.2.1   2016-07-03 CRAN (R 3.4.0) 
##  jsonlite      1.5     2017-06-01 CRAN (R 3.4.0) 
##  knitr       * 1.16    2017-05-18 CRAN (R 3.4.0) 
##  labeling      0.3     2014-08-23 CRAN (R 3.4.0) 
##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.0) 
##  lazyeval      0.2.0   2016-06-12 CRAN (R 3.4.0) 
##  lubridate     1.6.0   2016-09-13 CRAN (R 3.4.0) 
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0) 
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0) 
##  methods     * 3.4.0   2017-04-21 local          
##  mnormt        1.5-5   2016-10-15 CRAN (R 3.4.0) 
##  modelr        0.1.0   2016-08-31 CRAN (R 3.4.0) 
##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.0) 
##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.0) 
##  parallel      3.4.0   2017-04-21 local          
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0) 
##  plyr          1.8.4   2016-06-08 CRAN (R 3.4.0) 
##  psych         1.7.5   2017-05-03 CRAN (R 3.4.0) 
##  purrr       * 0.2.2.2 2017-05-11 cran (@0.2.2.2)
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
##  stats       * 3.4.0   2017-04-21 local          
##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.0) 
##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.0) 
##  tibble      * 1.3.3   2017-05-28 CRAN (R 3.4.0) 
##  tidyr       * 0.6.3   2017-05-15 CRAN (R 3.4.0) 
##  tidyverse   * 1.1.1   2017-01-27 CRAN (R 3.4.0) 
##  tools         3.4.0   2017-04-21 local          
##  utils       * 3.4.0   2017-04-21 local          
##  withr         1.0.2   2016-06-20 CRAN (R 3.4.0) 
##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.0)
```
