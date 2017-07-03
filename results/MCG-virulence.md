---
title: "Assessing Virulence by MCG"
output: 
  html_notebook:
    toc: true
---




```r
library("tidyverse")
```

# Purpose

This document is to assess virulence associated with the 10 most common MCGs. We
can get this from the cleaned data set that we saved in the
`data-comparison.Rmd` file. Note, we will be labelling the plot with 
"Aggressiveness" as that is the preferred term.


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
set.seed(2017-06-29)
ggplot(top_mcg, aes(x = MCG, y = Severity)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.5) +
  scale_y_continuous(limits = c(0, NA)) +
  # scale_x_discrete(position = "top") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(aspect.ratio = 1/2) +
  theme(axis.text = element_text(color = "black")) +
  # theme(axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_line(colour = "grey20")) +
  theme(panel.grid.minor = element_line(colour = "grey50")) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.border = element_blank()) +
  labs(list(
    # title = "Aggressiveness for the top 10 MCGs",
    y = "Aggressiveness",
    x = "Mycelial Compatibility Group"
    ))
```

![plot of chunk vis](./figures/MCG-virulence///vis-1.png)


```r
top_mcg %>% 
  group_by(MCG) %>%
  summarize(N = n(), 
            `Min Aggressiveness` = min(Severity),
            `Max Aggressiveness` = max(Severity),
            `Average Aggressiveness` = mean(Severity)
            ) %>%
  knitr::kable(digits = 2)
```



|MCG |  N| Min Aggressiveness| Max Aggressiveness| Average Aggressiveness|
|:---|--:|------------------:|------------------:|----------------------:|
|5   | 73|                3.6|                7.8|                   5.40|
|44  | 36|                4.3|                7.9|                   6.03|
|45  | 16|                3.9|                7.0|                   4.88|
|1   | 15|                3.9|                6.5|                   4.95|
|9   | 15|                2.8|                7.2|                   5.11|
|4   | 14|                3.9|                5.9|                   4.87|
|49  | 11|                3.3|                6.1|                   4.60|
|2   | 10|                4.1|                6.2|                   5.25|
|53  |  9|                3.6|                5.4|                   4.69|
|3   |  8|                4.7|                6.1|                   5.50|


## ANOVA

The default ANOVA in R sets contrasts as `contrast.treatment`, which compares
everything to the first factor, considered the treatment. Since we are
interested in whether or not there IS a difference between samples, this will
be sufficient.

After the ANOVA, we performed a Tukey's Honest Significant Difference test to
see exactly what groups these fell into.


```r
ANOVA <- aov(Severity ~ MCG, data = top_mcg)
ANOVA
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
summary(ANOVA)
```

```
##              Df Sum Sq Mean Sq F value   Pr(>F)    
## MCG           9  36.68   4.075   4.901 6.19e-06 ***
## Residuals   197 163.83   0.832                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
plot(TukeyHSD(ANOVA), las = 2)
```

![plot of chunk anova](./figures/MCG-virulence///anova-1.png)

```r
agricolae::HSD.test(ANOVA, "MCG")$groups
```

```
##    trt    means  M
## 1   44 6.027778  a
## 2   3  5.500000 ab
## 3   5  5.397260  b
## 4   2  5.250000  b
## 5   9  5.106667  b
## 6   1  4.953333  b
## 7   45 4.875000  b
## 8   4  4.871429  b
## 9   53 4.688889  b
## 10  49 4.600000  b
```

There appears to be a significant effect at p = 6.189e-06.


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
##  date     2017-07-03
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source         
##  agricolae     1.2-4   2016-06-12 CRAN (R 3.4.0) 
##  AlgDesign     1.1-7.3 2014-10-15 CRAN (R 3.4.0) 
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0) 
##  base        * 3.4.0   2017-04-21 local          
##  bindr         0.1     2016-11-13 CRAN (R 3.4.0) 
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.0) 
##  boot          1.3-19  2017-04-21 CRAN (R 3.4.0) 
##  broom         0.4.2   2017-02-13 CRAN (R 3.4.0) 
##  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.0) 
##  cluster       2.0.6   2017-03-16 CRAN (R 3.4.0) 
##  coda          0.19-1  2016-12-08 CRAN (R 3.4.0) 
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.0) 
##  combinat      0.0-8   2012-10-29 CRAN (R 3.4.0) 
##  compiler      3.4.0   2017-04-21 local          
##  datasets    * 3.4.0   2017-04-21 local          
##  deldir        0.1-14  2017-04-22 CRAN (R 3.4.0) 
##  devtools      1.13.2  2017-06-02 CRAN (R 3.4.0) 
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0) 
##  dplyr       * 0.7.1   2017-06-22 CRAN (R 3.4.0) 
##  evaluate      0.10    2016-10-11 CRAN (R 3.4.0) 
##  expm          0.999-2 2017-03-29 CRAN (R 3.4.0) 
##  ezknitr       0.6     2016-09-16 CRAN (R 3.4.0) 
##  forcats       0.2.0   2017-01-23 CRAN (R 3.4.0) 
##  foreign       0.8-69  2017-06-21 CRAN (R 3.4.0) 
##  gdata         2.18.0  2017-06-06 CRAN (R 3.4.0) 
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.0) 
##  glue          1.1.1   2017-06-21 CRAN (R 3.4.0) 
##  gmodels       2.16.2  2015-07-22 CRAN (R 3.4.0) 
##  graphics    * 3.4.0   2017-04-21 local          
##  grDevices   * 3.4.0   2017-04-21 local          
##  grid          3.4.0   2017-04-21 local          
##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.0) 
##  gtools        3.5.0   2015-05-29 CRAN (R 3.4.0) 
##  haven         1.0.0   2016-09-23 CRAN (R 3.4.0) 
##  highr         0.6     2016-05-09 CRAN (R 3.4.0) 
##  hms           0.3     2016-11-22 CRAN (R 3.4.0) 
##  httr          1.2.1   2016-07-03 CRAN (R 3.4.0) 
##  jsonlite      1.5     2017-06-01 CRAN (R 3.4.0) 
##  klaR          0.6-12  2014-08-06 CRAN (R 3.4.0) 
##  knitr       * 1.16    2017-05-18 CRAN (R 3.4.0) 
##  labeling      0.3     2014-08-23 CRAN (R 3.4.0) 
##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.0) 
##  lazyeval      0.2.0   2016-06-12 CRAN (R 3.4.0) 
##  LearnBayes    2.15    2014-05-29 CRAN (R 3.4.0) 
##  lubridate     1.6.0   2016-09-13 CRAN (R 3.4.0) 
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0) 
##  MASS          7.3-47  2017-04-21 CRAN (R 3.4.0) 
##  Matrix        1.2-10  2017-04-28 CRAN (R 3.4.0) 
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
##  utils       * 3.4.0   2017-04-21 local          
##  withr         1.0.2   2016-06-20 CRAN (R 3.4.0) 
##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.0)
```

</details>
