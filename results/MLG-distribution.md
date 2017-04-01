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

```r
library("ggraph")
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

I realized that it's possible to use an MLG table with matrix multiplication to
get an adjency matrix.


```r
datmlg  <- mlg.table(dat, plot = FALSE) > 0 # presence/absence of MLG
crosses <- mlg.crosspop(dat, quiet = TRUE, df = TRUE) %>% tbl_df()
adjmat  <- datmlg %*% t(datmlg) 
cols    <- sort(colnames(adjmat))
adjmat  <- adjmat[cols, cols]
adjmat
```

```
##    AU CA CO FR ID MI MN MX ND NE NY OR WA WI
## AU  6  0  0  1  0  0  4  0  0  0  0  1  0  0
## CA  0 15  0  0  0  0  0  0  0  0  0  3  9  0
## CO  0  0 28  0  0 11  1  0  5  7  0  0  7  0
## FR  1  0  0 14  0  4  1  0  0  0  0  0  3  0
## ID  0  0  0  0  1  0  0  0  0  0  0  0  0  0
## MI  0  0 11  4  0 43  1  0 13  8  0  0  8  0
## MN  4  0  1  1  0  1  7  0  0  1  0  0  1  0
## MX  0  0  0  0  0  0  0  9  0  0  0  0  0  0
## ND  0  0  5  0  0 13  0  0 35 11  0  0  3  0
## NE  0  0  7  0  0  8  1  0 11 28  0  4  7  0
## NY  0  0  0  0  0  0  0  0  0  0  1  1  0  0
## OR  1  3  0  0  0  0  0  0  0  4  1 13  1  0
## WA  0  9  7  3  0  8  1  0  3  7  0  1 56  0
## WI  0  0  0  0  0  0  0  0  0  0  0  0  0  2
```

```r
crosses
```

```
## # A tibble: 169 × 3
##       MLG Population Count
##    <fctr>     <fctr> <int>
## 1   MLG.1         OR     1
## 2   MLG.1         CA     1
## 3   MLG.2         MN     1
## 4   MLG.2         AU     1
## 5   MLG.3         MN     1
## 6   MLG.3         AU     1
## 7   MLG.4         MN     1
## 8   MLG.4         AU     1
## 9   MLG.4         FR     1
## 10  MLG.5         OR     1
## # ... with 159 more rows
```


Now that we have the adjacency matrix, we can use it to construct our graph:


```r
g           <- graph_from_adjacency_matrix(adjmat, mode = "undirected", diag = FALSE)
V(g)$size   <- diag(adjmat)
g           <- delete_vertices(g, degree(g) == 0)
shared_mlg  <- (crosses %>% group_by(Population) %>% summarize(n = n()))$n
V(g)$weight <- 1 - shared_mlg/V(g)$size
el          <- as_adj_edge_list(g)
el          <- el[lengths(el) > 0]
popgraphs <- setNames(vector(mode = "list", length = length(el)), names(el))
for (v in names(el)){
  idx  <- el[[v]]
  mlgs <- crosses %>%           # How to get all MLGs from a single population:
    filter(Population == v) %>%         # Grab only the population e and then
    select(MLG) %>%                     # remove everything but the MLGs to do an
    inner_join(crosses, by = "MLG") %>% # inner join of the original list and then
    filter(Population != v) %>%         # remove the query population to give
    arrange(Population)                 # the neigboring populations in order.
  MLGS <- as.character(mlgs$MLG)
  E(g)[idx]$label  <- substr(MLGS, 5, nchar(MLGS))
  E(g)[idx]$weight <- as.integer(table(MLGS)[MLGS]) # weight == n populations visited
  popgraphs[[v]]   <- subgraph.edges(g, eids = idx)
}
par(mfrow = c(3, 4))
for (i in names(popgraphs)){
  pg <- popgraphs[[i]]
  labs <- ifelse(E(pg)$weight > 1, E(pg)$label, NA)
  labs <- ifelse(duplicated(labs), NA, labs)
  plot(pg, 
       main = i, 
       layout = layout_as_star(pg, center = i), 
       edge.width = E(pg)$weight,
       edge.label = labs)
}
par(mfrow = c(1, 1))
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-1.png)

```r
center_node <- degree(g) %>% which.max() %>% names()

new_layout <- structure(c(1, 0.809016994374947, 0.309016994374947, -0.309016994374947, 
-0.809016994374947, -1, -0.809016994374947, -0.309016994374948, 
0, 0.309016994374947, 0.809016994374947, 0, 0.587785252292473, 
0.951056516295154, 0.951056516295154, 0.587785252292473, 1.22464679914735e-16, 
-0.587785252292473, -0.951056516295154, 0, -0.951056516295154, 
-0.587785252292473), .Dim = c(11L, 2L), .Dimnames = list(c("MN", 
"AU", "OR", "FR", "WA", "CA", "NE", "CO", "MI", "ND", "NY"), 
    c("x", "y")))
new_layout <- new_layout[V(g)$name, ]
g <- add_vertices(g, length(V(g)), size = V(g)$size - shared_mlg, color = "grey90")
plot(g, 
     layout = rbind(new_layout, new_layout),
     edge.width = E(g)$weight,
     edge.label = NA)
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-2.png)

```r
glay <- create_layout(g, "manual", node.positions = as.data.frame(rbind(new_layout, new_layout)))
ggraph(glay) +
  geom_edge_fan(aes(alpha = weight + 1)) +
  geom_node_circle(aes(r = scale(size, center = FALSE)/10, fill = size, alpha = weight)) +
  geom_node_label(aes(label = name), repel = TRUE) +
  viridis::scale_fill_viridis(option = "C") +
  coord_fixed() +
  theme_void() +
  labs(list(
    title = "Shared haplotypes across regions",
    fill = "Number of\nGenotypes",
    alpha = "Fraction of\nprivate genotpes",
    edge_alpha = "Populations\nper haploytpe",
    caption = "Outer circle: Number of haplotypes in the region\nInner Circle: Number of private haplotypes in the region"
  ))
```

```
## Warning: Removed 11 rows containing missing values (geom_label_repel).
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-3.png)

