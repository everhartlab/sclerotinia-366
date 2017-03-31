---
title: "R Notebook"
output: 
  html_notebook:
    toc: true
---



In this document, I will create a graph that shows the distribution of MLGs
across populations


```r
library('tidyverse')
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

```r
library('assertr')
library('poppr')
```

```
## Loading required package: adegenet
```

```
## Loading required package: ade4
```

```
## 
##    /// adegenet 2.1.0 is loaded ////////////
## 
##    > overview: '?adegenet'
##    > tutorials/doc/questions: 'adegenetWeb()' 
##    > bug reports/feature requests: adegenetIssues()
```

```
## This is poppr version 2.3.0.99.42. To get started, type package?poppr
## OMP parallel support: available
## 
## This version of poppr is under development.
## If you find any bugs, please report them at https://github.com/grunwaldlab/poppr/issues
```

```r
library('igraph')
```

```
## 
## Attaching package: 'igraph'
```

```
## The following object is masked from 'package:poppr':
## 
##     %>%
```

```
## The following objects are masked from 'package:dplyr':
## 
##     %>%, as_data_frame, groups, union
```

```
## The following objects are masked from 'package:purrr':
## 
##     %>%, compose, simplify
```

```
## The following objects are masked from 'package:tidyr':
## 
##     %>%, crossing
```

```
## The following object is masked from 'package:tibble':
## 
##     as_data_frame
```

```
## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum
```

```
## The following object is masked from 'package:base':
## 
##     union
```

## Loading data and assertions






```r
dat <- read.genalex("../Analysis4 ForManu/A2_Copy4 EUR_AUS_forManu.csv", ploidy = 1)
splitStrata(dat) <- ~Isolate/Severity/MCG/Region/Source/Year/Host
dat
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
##      7 strata - Isolate, Severity, MCG, Region, Source, Year, Host
##    366 populations defined - 
## 152_3.9_4_NE_unk_2003_GH, 274_5.4_45_NE_unk_2003_GH, 443_6.3_5_NY_unk_2003_GH, ..., 967_5.8_34_FR_flds_2012_unk, 968_4.2_34_FR_flds_2012_unk, 970_5.2_35_FR_flds_2012_unk
```

The incoming strata includes both Severity and Isolate. Since these are not
necessary for delimiting the strata, we will place them in the "other" slot
after converting Severity to numeric. Placing this information in the "other"
slot ensures that these data will travel with the object.


```r
dat_strata <- strata(dat) %>%
  mutate_all(as.character) %>%
  mutate(Severity = as.numeric(Severity))
strata(dat)     <- select(dat_strata, -Severity, -Isolate)
indNames(dat)   <- dat_strata$Isolate
other(dat)$meta <- select(dat_strata, Severity, Isolate)
```


```r
setPop(dat) <- ~Region
dat
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
##     14 populations defined - NE, NY, MN, ..., FR, MX, ND
```

## Crossing populations


We can use `mlg.crosspop()` to tabulte which MLGs cross populations.


```r
make_from_to <- function(poplist){
  x <- combn(poplist, 2)
  tibble::data_frame(from = x[1, ], to = x[2, ])
}
crosses <- mlg.crosspop(dat, df = TRUE, quiet = TRUE)

pg <- 
  crosses %>%
  group_by(MLG) %>% 
  summarize(Population = list(make_from_to(Population)), Size = length(Count))
  
pop_graph <- pg %>%  unnest() %>%
  group_by(from, to) %>%
  mutate(weight = length(Size)) %>%
  ungroup() %>%
  select(from, to, weight, everything())

unchoose <- function(x) ceiling(sqrt(x * 2))


  
gcross           <- graph_from_data_frame(pop_graph, directed = FALSE)
E(gcross)$name   <- as.character(pop_graph$MLG)
E(gcross)$weight <- unchoose(table(E(gcross)$name)[E(gcross)$name])
colors <- viridis::magma(max(E(gcross)$weight) - 1, end = 0.8)
set.seed(50)
gcross %>%
  plot(., vertex.size = table(pop(clonecorrect(dat, ~Region)))[names(V(gcross))],
       edge.width = E(.)$weight - 1,
       layout = layout_as_star(., center = "MI"),
       # layout = layout_in_circle(.),
       edge.color = colors[E(.)$weight - 1])
```

![plot of chunk unnamed-chunk-3](./figures/MLG-distribution///unnamed-chunk-3-1.png)




```r
pops <- crosses %>%
  group_by(Population) %>%
  summarize(MLG = list(as.character(MLG))) %>%
  group_by(Population) %>%
  mutate(interactions = list(pg$Population[pg$MLG %in% unlist(MLG)] %>% setNames(unlist(MLG)) %>% bind_rows(.id = "MLG") %>% distinct)) %>%
  mutate(graph = map(interactions, ~select(.x, from, to, MLG) %>% graph_from_data_frame(directed = FALSE)))
psize <- table(pop(clonecorrect(dat, ~Region)))
par(mfrow = c(3, 4))

apply(pops, 1, function(i){
    plot(i$graph, 
         layout = layout_as_star(i$graph, center = i$Population), 
         main = i$Population, 
         vertex.size = psize[names(V(i$graph))])
  })
```

```
## NULL
```

```r
par(mfrow = c(1, 1))
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-1.png)

```r
do.call("union", pops$graph) %>% 
  plot(., vertex.size = psize[names(V(.))],
       layout = layout_as_star(., center = "MI"))
```

```
## Error in union(structure(list(4, FALSE, c(1, 1, 1, 1, 3, 3, 2), c(0, 0, : unused arguments (list(6, FALSE, c(2, 2, 5, 5, 3, 3, 3, 3, 5, 2, 5, 5, 3, 2, 3, 3, 2, 3, 2, 5, 3, 5, 5, 4, 3, 2, 4, 4, 3, 2, 3, 1, 2, 5, 3, 3, 5, 2, 5, 5, 2, 5, 5, 2, 2, 3, 3, 2, 3, 3), c(0, 1, 1, 2, 2, 2, 0, 2, 3, 0, 0, 2, 1, 1, 2, 2, 0, 1, 1, 1, 2, 3, 2, 1, 1, 1, 3, 2, 2, 0, 1, 0, 1, 1, 0, 2, 3, 0, 0, 2, 0, 0, 2, 0, 1, 2, 1, 1, 2, 2), c(31, 43, 40, 37, 29, 16, 9, 0, 47, 44, 32, 25, 18, 13, 1, 34, 6, 46, 30, 24, 17, 12, 49, 48, 45, 35, 28, 20, 15, 14, 7, 5, 4, 23, 27, 26, 41, 38, 10, 33, 19, 2, 42, 39, 22, 
## 11, 3, 36, 21, 8), c(31, 43, 40, 37, 29, 16, 9, 0, 34, 6, 41, 38, 10, 47, 44, 32, 25, 18, 13, 1, 46, 30, 24, 17, 12, 23, 33, 19, 2, 49, 48, 45, 35, 28, 20, 15, 14, 7, 5, 4, 27, 42, 39, 22, 11, 3, 26, 36, 21, 8), c(0, 0, 1, 15, 33, 36, 50), c(0, 13, 29, 46, 50, 50, 50), list(c(1, 0, 1), list(), list(name = c("WA", "NE", "CO", "MI", "MN", "ND")), list(MLG = c("MLG.15", "MLG.22", "MLG.22", "MLG.22", "MLG.25", "MLG.39", "MLG.44", "MLG.44", "MLG.44", "MLG.44", "MLG.44", "MLG.44", "MLG.56", "MLG.56", "MLG.56", 
## "MLG.60", "MLG.66", "MLG.67", "MLG.67", "MLG.67", "MLG.67", "MLG.67", "MLG.67", "MLG.106", "MLG.106", "MLG.106", "MLG.106", "MLG.106", "MLG.106", "MLG.128", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.138", "MLG.140", "MLG.140", "MLG.140", "MLG.153", "MLG.155", "MLG.162", "MLG.163", "MLG.163", "MLG.163", "MLG.165"))), <environment>), list(5, FALSE, c(4, 4, 4, 2, 4, 4, 4, 4, 4, 4), c(0, 0, 0, 1, 1, 2, 3, 3, 3, 0), c(3, 9, 2, 1, 0, 4, 5, 8, 
## 7, 6), c(9, 2, 1, 0, 3, 4, 5, 8, 7, 6), c(0, 0, 0, 1, 1, 10), c(0, 4, 6, 7, 10, 10), list(c(1, 0, 1), list(), list(name = c("MI", "MN", "AU", "WA", "FR")), list(MLG = c("MLG.4", "MLG.85", "MLG.86", "MLG.87", "MLG.87", "MLG.87", "MLG.88", "MLG.109", "MLG.126", "MLG.146"))), <environment>), list(7, FALSE, c(5, 6, 1, 1, 6, 1, 2, 2, 5, 6, 1, 2, 5, 2, 5, 5, 5, 3, 3, 2, 5, 2, 5, 3, 3, 5, 2, 5, 5, 4, 3, 3, 4, 4, 2, 3, 3, 3, 5, 1, 2, 5, 2, 5, 5, 1, 2, 3, 3, 3, 2, 2, 1, 5, 5, 5, 3, 5, 5, 1, 6, 3, 5, 5), c(0, 
## 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 3, 0, 0, 2, 3, 0, 2, 0, 2, 0, 0, 1, 2, 3, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0), c(59, 52, 45, 39, 10, 5, 3, 2, 51, 50, 46, 40, 34, 26, 21, 19, 11, 7, 6, 42, 13, 61, 56, 48, 47, 35, 30, 23, 17, 36, 49, 37, 31, 24, 18, 32, 33, 29, 63, 58, 55, 54, 53, 41, 27, 22, 20, 16, 12, 8, 0, 43, 14, 44, 28, 15, 62, 57, 38, 25, 60, 9, 4, 1), c(59, 52, 45, 39, 10, 5, 3, 2, 51, 50, 46, 40, 34, 26, 21, 19, 
## 11, 7, 6, 61, 56, 48, 47, 35, 30, 23, 17, 32, 63, 58, 55, 54, 53, 41, 27, 22, 20, 16, 12, 8, 0, 60, 9, 4, 1, 42, 13, 36, 43, 14, 49, 37, 31, 24, 18, 33, 44, 28, 15, 29, 62, 57, 38, 25), c(0, 0, 8, 21, 35, 38, 60, 64), c(0, 45, 50, 59, 64, 64, 64, 64), list(c(1, 0, 1), list(), list(name = c("MI", "WA", "CO", "NE", "MN", "ND", "FR")), list(MLG = c("MLG.9", "MLG.12", "MLG.15", "MLG.25", "MLG.35", "MLG.60", "MLG.65", "MLG.66", "MLG.67", "MLG.73", "MLG.74", "MLG.74", "MLG.74", "MLG.74", "MLG.74", "MLG.74", 
## "MLG.75", "MLG.76", "MLG.76", "MLG.76", "MLG.78", "MLG.80", "MLG.88", "MLG.102", "MLG.102", "MLG.102", "MLG.102", "MLG.102", "MLG.102", "MLG.109", "MLG.109", "MLG.109", "MLG.109", "MLG.109", "MLG.109", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.125", "MLG.126", "MLG.127", "MLG.138", "MLG.140", "MLG.140", "MLG.140", "MLG.145", "MLG.146", "MLG.152", "MLG.153", "MLG.154", "MLG.155", "MLG.155", "MLG.155", "MLG.161", "MLG.163", "MLG.165", "MLG.165", 
## "MLG.165"))), <environment>), list(7, FALSE, c(4, 3, 1, 2, 5, 2, 5, 5, 3, 3, 3, 6, 6), c(0, 0, 0, 1, 1, 0, 0, 2, 0, 0, 0, 0, 3), c(2, 5, 3, 10, 9, 8, 1, 0, 6, 4, 7, 11, 12), c(2, 5, 10, 9, 8, 1, 0, 6, 11, 3, 4, 7, 12), c(0, 0, 1, 3, 7, 8, 11, 13), c(0, 9, 11, 12, 13, 13, 13, 13), list(c(1, 0, 1), list(), list(name = c("MN", "NE", "MI", "AU", "WA", "CO", "FR")), list(MLG = c("MLG.2", "MLG.3", "MLG.4", "MLG.4", "MLG.4", "MLG.4", "MLG.4", "MLG.4", "MLG.123", "MLG.141", "MLG.165", "MLG.165", "MLG.165"
## ))), <environment>), list(5, FALSE, c(4, 4, 4, 2, 4, 4, 4, 4, 4, 3, 2, 4, 3, 4, 4, 4, 4, 4, 4, 1, 2, 4, 2, 4, 4, 1, 3, 2, 4, 3, 2, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 1, 4, 4, 1, 4, 4), c(0, 1, 1, 1, 1, 2, 1, 1, 0, 0, 0, 0, 2, 3, 2, 0, 0, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 1, 1, 0, 0, 0, 2, 3, 2, 2, 3, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0), c(45, 42, 25, 19, 30, 22, 10, 27, 20, 3, 29, 9, 26, 35, 32, 12, 47, 44, 41, 40, 39, 31, 23, 18, 16, 15, 11, 8, 0, 46, 43, 38, 28, 21, 17, 7, 6, 4, 2, 1, 37, 34, 24, 14, 5, 
## 36, 33, 13), c(45, 42, 25, 19, 30, 22, 10, 29, 9, 47, 44, 41, 40, 39, 31, 23, 18, 16, 15, 11, 8, 0, 27, 20, 3, 26, 46, 43, 38, 28, 21, 17, 7, 6, 4, 2, 1, 35, 32, 12, 37, 34, 24, 14, 5, 36, 33, 13), c(0, 0, 4, 10, 16, 48), c(0, 22, 37, 45, 48, 48), list(c(1, 0, 1), list(), list(name = c("MI", "NE", "CO", "WA", "ND")), list(MLG = c("MLG.7", "MLG.9", "MLG.15", "MLG.25", "MLG.25", "MLG.25", "MLG.39", "MLG.74", "MLG.75", "MLG.76", "MLG.76", "MLG.76", "MLG.76", "MLG.76", "MLG.76", "MLG.78", "MLG.102", 
## "MLG.104", "MLG.124", "MLG.128", "MLG.128", "MLG.128", "MLG.128", "MLG.128", "MLG.128", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.136", "MLG.144", "MLG.144", "MLG.144", "MLG.145", "MLG.152", "MLG.154", "MLG.156", "MLG.161", "MLG.161", "MLG.161", "MLG.163", "MLG.163", "MLG.163"))), <environment>), list(7, FALSE, c(5, 6, 6, 1, 6, 6, 3, 6, 6, 2, 1, 2, 6, 3, 5, 5, 2, 1, 6, 2, 6, 6, 4, 2, 1, 4, 4, 2, 3, 3, 2, 5, 1, 6, 5, 2, 6, 5, 6, 6, 5, 1, 
## 5, 5, 2, 2, 1, 2, 6, 5, 2, 6, 6, 2, 6, 6), c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 2, 1, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 2, 1, 2, 1, 5, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 2), c(46, 41, 32, 24, 17, 10, 3, 53, 50, 45, 44, 30, 23, 16, 9, 47, 35, 27, 19, 11, 29, 28, 13, 6, 22, 26, 25, 49, 43, 42, 40, 31, 14, 0, 37, 34, 15, 54, 51, 48, 33, 18, 12, 8, 7, 4, 2, 1, 39, 21, 5, 55, 52, 36, 20, 38), c(46, 41, 32, 24, 17, 10, 3, 53, 50, 45, 44, 30, 23, 16, 9, 29, 28, 13, 
## 6, 22, 49, 43, 42, 40, 31, 14, 0, 54, 51, 48, 33, 18, 12, 8, 7, 4, 2, 1, 47, 35, 27, 19, 11, 26, 37, 39, 21, 5, 25, 34, 55, 52, 36, 20, 15, 38), c(0, 0, 7, 20, 24, 27, 37, 56), c(0, 38, 48, 54, 55, 55, 56, 56), list(c(1, 0, 1), list(), list(name = c("NE", "CO", "MI", "OR", "MN", "WA", "ND")), list(MLG = c("MLG.7", "MLG.9", "MLG.17", "MLG.24", "MLG.24", "MLG.24", "MLG.25", "MLG.47", "MLG.56", "MLG.62", "MLG.62", "MLG.62", "MLG.63", "MLG.65", "MLG.65", "MLG.65", "MLG.66", "MLG.66", "MLG.66", "MLG.66", 
## "MLG.66", "MLG.66", "MLG.77", "MLG.77", "MLG.77", "MLG.77", "MLG.77", "MLG.77", "MLG.78", "MLG.103", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.104", "MLG.124", "MLG.128", "MLG.134", "MLG.136", "MLG.144", "MLG.153", "MLG.153", "MLG.153", "MLG.156", "MLG.160", "MLG.163", "MLG.163", "MLG.163", "MLG.165", "MLG.165", "MLG.165"))), <environment>), list(2, FALSE, 1, 0, 0, 0, c(0, 0, 1), c(0, 1, 1), list(c(1, 0, 1), list(), list(name = c("NY", 
## "OR")), list(MLG = "MLG.42")), <environment>), list(6, FALSE, c(3, 1, 1, 4, 4, 1, 1, 5, 2, 3, 3), c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), c(6, 5, 2, 1, 8, 10, 9, 0, 4, 3, 7), c(6, 5, 2, 1, 8, 10, 9, 0, 4, 7, 3), c(0, 0, 4, 5, 8, 10, 11), c(0, 10, 11, 11, 11, 11, 11), list(c(1, 0, 1), list(), list(name = c("OR", "NE", "NY", "CA", "WA", "AU")), list(MLG = c("MLG.1", "MLG.5", "MLG.17", "MLG.17", "MLG.17", "MLG.24", "MLG.41", "MLG.42", "MLG.83", "MLG.134", "MLG.160"))), <environment>), list(9, FALSE, c(1, 
## 4, 6, 6, 2, 6, 3, 2, 2, 2, 4, 7, 4, 7, 7, 5, 1, 5, 4, 4, 2, 1, 4, 7, 2, 4, 7, 4, 7, 7, 2, 4, 7, 7, 4, 1, 6, 1, 1, 6, 6, 2, 1, 6, 2, 6, 6, 8, 8, 8), c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 4, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 2, 2, 1, 1, 4, 1, 1, 1, 4, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1), c(42, 38, 37, 35, 21, 16, 0, 20, 44, 41, 30, 24, 9, 8, 7, 4, 6, 22, 34, 31, 27, 19, 18, 12, 1, 25, 10, 15, 17, 46, 45, 43, 40, 39, 36, 5, 3, 2, 23, 32, 28, 13, 26, 11, 33, 29, 14, 49, 48, 47), c(42, 38, 
## 37, 35, 21, 16, 0, 20, 22, 15, 23, 44, 41, 30, 24, 9, 8, 7, 4, 6, 34, 31, 27, 19, 18, 12, 1, 17, 46
```

