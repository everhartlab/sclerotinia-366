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
## This is poppr version 2.3.0.99.45. To get started, type package?poppr
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
library('magrittr')
```

```
## 
## Attaching package: 'magrittr'
```

```
## The following object is masked from 'package:igraph':
## 
##     %>%
```

```
## The following object is masked from 'package:purrr':
## 
##     set_names
```

```
## The following object is masked from 'package:tidyr':
## 
##     extract
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

```r
repeat_lengths <-
  c(
  `5-2` = 2.000000,
  `5-3` = 0.500000,
  `6-2` = 6.000000,
  `7-2` = 2.000000,
  `8-3` = 2.000000,
  `9-2` = 2.000000,
  `12-2` = 2.000000,
  `17-3` = 3.000000,
  `20-3` = 2.000000,
  `36-4` = 0.500000,
  `50-4` = 0.500000,
  `55-4` = 4.000000,
  `92-4` = 2.000000,
  `106-4` = 4.000000,
  `110-4` = 4.000000,
  `114-4` = 4.000000
  )
repeat_lengths <- ifelse(repeat_lengths < 1, 4, repeat_lengths)
```



```r
ex <- readxl::read_excel("../Analysis4 ForManu/A1_Copy of binned-genotypes_SE.xlsx", sheet = "GenAlex", skip = 1) %>%
  select(-1) %>%                # removing first column, which is empty
  gather(locus, allele, -1) %>% # gather all loci into tidy columns
  mutate(locus = trimws(locus)) %>% # remove (F) designator
  mutate(allele = as.integer(allele)) %>% # force alleles to integers
  spread(locus, allele) %>%
  slice(-n())
ex <- ex[!names(ex) %in% locNames(dat)]

# Function to select an adjacent allele. It will select the
# next allele if the next allele is not missing and it's distance
# is one away and the previous allele for the same conditions.
# If none of the conditions are met, it will retain the allele.
cromulent_allele <- Vectorize(function(lower, allele, higher){
  if (!is.na(higher) && abs(allele - higher) == 1){
    out <- higher
  } else if (!is.na(lower) && abs(allele - lower) == 1){
    out <- lower
  } else {
    out <- allele
  }
  out
})
ex
```

```
## # A tibble: 366 × 6
##    iso_st_mcg_org_loc_yr_hst_cult_rep `106-4(H)` `36-4(F)` `5-3(F)`
##                                 <chr>      <int>     <int>    <int>
## 1               152_3.9_4_NE__2003_GH        580       415      328
## 2              274_5.4_45_NE__2003_GH        588       415      328
## 3               443_6.3_5_NY__2003_GH        567       415      308
## 4          444_4.4_4_MN_wmn_2003_G122        580       415      328
## 5         445_4.7_4_MN_wmn_2003_Beryl        580       415      328
## 6         446_6.1_3_MI_wmn_2003_Beryl        567       415      339
## 7         447_5.5_5_MI_wmn_2003_Beryl        567       415      308
## 8           448_5_3_MI_wmn_2003_Beryl        568       414      339
## 9         449_5.2_3_MI_wmn_2003_Bunsi        568       415      339
## 10        450_5.3_5_MI_wmn_2003_Bunsi        568       415      308
## # ... with 356 more rows, and 2 more variables: `50-4(F)` <int>,
## #   `92-4(F)` <int>
```

```r
exsummary <- ex %>% 
  gather(locus, allele, -1) %>% # tidy the data
  group_by(locus, allele) %>%   
  summarize(n = n()) %>%        # summarize by count 
  ungroup() %>%
  group_by(locus) %>%           # group the loci, add the lower and upper alleles,
  mutate(lower = lag(allele), higher = lead(allele)) %>% # and then create new_alleles
  mutate(new_allele = ifelse(n < 3, cromulent_allele(lower, allele, higher), allele)) %>%
  select(locus, new_allele, allele)
exsummary
```

```
## Source: local data frame [71 x 3]
## Groups: locus [5]
## 
##       locus new_allele allele
##       <chr>      <int>  <int>
## 1  106-4(H)        502    501
## 2  106-4(H)        502    502
## 3  106-4(H)        502    503
## 4  106-4(H)        511    511
## 5  106-4(H)        532    532
## 6  106-4(H)        533    533
## 7  106-4(H)        541    540
## 8  106-4(H)        541    541
## 9  106-4(H)        541    542
## 10 106-4(H)        546    546
## # ... with 61 more rows
```

```r
corrected_loci <- ex %>% gather(locus, allele, -1) %>%
  left_join(exsummary, by = c("locus", "allele")) %>%
  mutate(allele = new_allele) %>%
  select(-new_allele) %>%
  spread(locus, allele)
datdf <- genind2df(dat, usepop = FALSE) %>% bind_cols(corrected_loci[-1])
dat  <- df2genind(datdf, ind.names = indNames(dat), strata = strata(dat), ploidy = 1) %>% as.genclone()
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
##    215 original multilocus genotypes 
##    366 haploid individuals
##     16 codominant loci
## 
## Population information:
## 
##      5 strata - MCG, Region, Source, Year, Host
##     14 populations defined - NE, NY, MN, ..., FR, MX, ND
```

```r
keeploci <- !locNames(dat) %in% colnames(corrected_loci)
genotype_curve(dat, sample = 1000)
```

![plot of chunk unnamed-chunk-2](./figures/MLG-distribution///unnamed-chunk-2-1.png)

```r
locus_table(dat)
```

```
## 
## allele = Number of observed alleles
```

```
## 
## 1-D = Simpson index
```

```
## 
## Hexp = Nei's 1978 gene diversity
```

```
## ------------------------------------------
```

```
##           summary
## locus      allele    1-D   Hexp Evenness
##   5-2(F)    4.000  0.451  0.452    0.616
##   6-2(F)    3.000  0.643  0.645    0.949
##   7-2(F)    7.000  0.727  0.729    0.764
##   8-3(H)    7.000  0.740  0.742    0.789
##   9-2(F)    9.000  0.347  0.348    0.406
##   12-2(H)   5.000  0.579  0.580    0.779
##   17-3(H)   7.000  0.551  0.553    0.526
##   20-3(F)   2.000  0.053  0.053    0.420
##   55-4(F)  10.000  0.721  0.723    0.656
##   110-4(H)  5.000  0.754  0.756    0.915
##   114-4(H) 10.000  0.828  0.831    0.801
##   106-4(H) 32.000  0.917  0.919    0.602
##   36-4(F)   4.000  0.254  0.255    0.501
##   5-3(F)   12.000  0.839  0.841    0.791
##   50-4(F)   3.000  0.256  0.257    0.629
##   92-4(F)   9.000  0.797  0.799    0.808
##   mean      8.062  0.591  0.593    0.685
```

```r
# rl <- fix_replen(dat, repeat_lengths)
# mlg.filter(dat, distance = bruvo.dist, replen = rl) <- .Machine$double.eps
dat
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    215 original multilocus genotypes 
##    366 haploid individuals
##     16 codominant loci
## 
## Population information:
## 
##      5 strata - MCG, Region, Source, Year, Host
##     14 populations defined - NE, NY, MN, ..., FR, MX, ND
```

```r
genotype_curve(dat[loc = keeploci, mlg.reset = TRUE])
```

![plot of chunk unnamed-chunk-2](./figures/MLG-distribution///unnamed-chunk-2-2.png)

```r
# dat <- dat[loc = keeploci, mlg.reset = TRUE]
dat
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    215 original multilocus genotypes 
##    366 haploid individuals
##     16 codominant loci
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
make_graph_list <- function(dat){ # dat is a genclone/snpclone object
  datmlg  <- mlg.table(dat, plot = FALSE) > 0 # presence/absence of MLG
  crosses <- mlg.crosspop(dat, quiet = TRUE, df = TRUE) %>% tbl_df()
  adjmat  <- datmlg %*% t(datmlg) 
  cols    <- sort(colnames(adjmat))
  adjmat  <- adjmat[cols, cols]
  # Creating Graph
  g           <- graph_from_adjacency_matrix(adjmat, mode = "undirected", diag = FALSE)
  V(g)$size   <- diag(adjmat)
  g           <- delete_vertices(g, degree(g) == 0)
  shared_mlg  <- (crosses %>% group_by(Population) %>% summarize(n = n()))$n
  V(g)$weight <- 1 - shared_mlg/V(g)$size
  el          <- as_adj_edge_list(g)
  el          <- el[lengths(el) > 0]
  popgraphs <- setNames(vector(mode = "list", length = length(el) + 1), c(names(el), "total"))
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
  popgraphs[["total"]] <- g
  popgraphs
}

plot_mlg_graph <- function(g, glayout = NULL){
  if (is.null(glayout)){
    glayout <- layout_nicely(g)
  } else {
    glayout <- glayout[V(g)$name, ]
  }
  shared_mlg  <- (1 - V(g)$weight) * V(g)$size
  g           <- add_vertices(g, length(V(g)), size = V(g)$size - shared_mlg, color = "grey90")
  g2          <- g
  V(g2)$label <- ifelse(!is.na(V(g2)$name), sprintf("%s [(%d/%d)]", V(g2)$name, shared_mlg, V(g2)$size), NA)
  glay        <- create_layout(g2, "manual", node.positions = as.data.frame(rbind(glayout, glayout)))
  x_nudge     <- ifelse(abs(glay$x) == 1, -glay$x/10, glay$x/10)
  ggraph(glay) +
    geom_edge_fan(aes(alpha = weight + 1), width = 0.75) +
    geom_node_circle(aes(r = drop(scale(size, center = FALSE))/10, fill = size)) +
    geom_node_label(aes(label = label), repel = TRUE, parse = TRUE, label.size = 0.75, nudge_x = x_nudge) +
    viridis::scale_fill_viridis(option = "C") +
    scale_edge_alpha_continuous(range = c(0.25, 1), breaks = c(2:5)) +
    coord_fixed() +
    theme_void() +
    theme(text = element_text(size = 14)) +
    labs(list(
      title = "Shared haplotypes across regions",
      fill = "Number of\nGenotypes",
      edge_alpha = "Populations\nper haploytpe",
      caption = "Outer circle: Number of haplotypes in the region\nInner Circle: Number of private haplotypes in the region"
    ))
}

plot_mlg_subgraph <- function(graphlist){
  for (i in names(graphlist)){
    pg   <- graphlist[[i]]
    labs <- ifelse(E(pg)$weight > 1, E(pg)$label, NA)
    labs <- ifelse(duplicated(labs), NA, labs)
    plot(pg, 
         main = i, 
         layout = layout_as_star(pg, center = i), 
         edge.width = E(pg)$weight,
         edge.label = labs)
  }
}



good_layout <- read.table(
text = 
"   x                            y
AU  0.8090169943749470071737  5.877852522924730260812e-01
CA -1.0000000000000000000000  1.224646799147350002426e-16
CO -0.3090169943749480063744 -9.510565162951539752711e-01
FR -0.3090169943749470071737  9.510565162951539752711e-01
MI  0.0000000000000000000000  0.000000000000000000000e+00
MN  1.0000000000000000000000  0.000000000000000000000e+00
ND  0.3090169943749470071737 -9.510565162951539752711e-01
NE -0.8090169943749470071737 -5.877852522924730260812e-01
NY  0.8090169943749470071737 -5.877852522924730260812e-01
OR  0.3090169943749470071737  9.510565162951539752711e-01
WA -0.8090169943749470071737  5.877852522924730260812e-01"
) %>% as.matrix()
```


## Graphs

With the fuctions above, we can create and plot the graphs


```r
# Creating the graphs
graph16loc <- make_graph_list(dat)
graph11loc <- make_graph_list(dat[loc = keeploci, mlg.reset = TRUE])

# Plotting the subgraphs
par(mfrow = c(3, 4))
plot_mlg_subgraph(graph16loc[-length(graph16loc)])

par(mfrow = c(3, 4))
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-1.png)

```r
plot_mlg_subgraph(graph11loc[-length(graph11loc)])

# Plotting the full graphs
par(mfrow = c(1, 1))
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-2.png)

```r
plot_mlg_graph(graph16loc$total, good_layout)
```

```
## Warning: Removed 10 rows containing missing values (geom_label_repel).
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-3.png)

```r
plot_mlg_graph(graph11loc$total, good_layout)
```

```
## Warning: Removed 11 rows containing missing values (geom_label_repel).
```

![plot of chunk unnamed-chunk-4](./figures/MLG-distribution///unnamed-chunk-4-4.png)


```r
null_vertices <- length(V(g)) %>% seq(from = ./2 + 1, to = ., by = 1)
sg <- g %>% 
  delete_vertices(null_vertices) %>%
  cluster_spinglass()
modularity(sg)
communities(sg)
```

