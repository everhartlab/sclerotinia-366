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
igraph::union(pops$graph[[1]], pops$graph[-1]) %>% 
  plot(., vertex.size = psize[names(V(.))],
       layout = layout_as_star(., center = "MI"))
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-2.png)

