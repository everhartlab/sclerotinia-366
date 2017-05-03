---
title: "Multilocus Genotype And Mycelial Compatability Group Assessment"
output: 
  html_notebook:
    toc: true
---








# Introduction 

The purpose of this document is to assess the distribution of the MLGs and MCGs
within the data. Specifically, we want to know if there are any MLGs that 
consistently coordinate with a single MCG, or if there are anything close.

## Packages and Data


```r
library('igraph')
```

```
## 
## Attaching package: 'igraph'
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
library('ggraph')
```

```
## Loading required package: ggplot2
```

```r
library('tidyverse')
```

```
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
## as_data_frame(): dplyr, tibble, igraph
## compose():       purrr, igraph
## crossing():      tidyr, igraph
## filter():        dplyr, stats
## groups():        dplyr, igraph
## lag():           dplyr, stats
## simplify():      purrr, igraph
```

```r
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
## This is poppr version 2.4.1. To get started, type package?poppr
## OMP parallel support: available
```

```
## 
## Attaching package: 'poppr'
```

```
## The following object is masked from 'package:igraph':
## 
##     %>%
```

```r
library('viridis')
```

```
## Loading required package: viridisLite
```


```r
load(file.path(PROJHOME, "data", "sclerotinia_16_loci.rda"))
load(file.path(PROJHOME, "data", "mlg-crosspop-graph.rda"))
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
g11o <- igraph::cluster_optimal(graph11loc$total)
comm <- igraph::communities(g11o)
names(comm) <- c("International", "Costal", "Midwest")
comm
```

```
## $International
## [1] "Australia" "France"    "MN"       
## 
## $Costal
## [1] "CA" "NY" "OR" "WA"
## 
## $Midwest
## [1] "CO" "MI" "ND" "NE"
```

```r
strata(dat11) <- strata(dat11) %>%
  mutate(MLGRegion = case_when(
    .$Region %in% comm$International ~ "International",
    .$Region %in% comm$Costal ~ "Costal",
    .$Region %in% comm$Midwest ~ "Midwest",
    TRUE ~ as.character(.$Region)
  ))

setPop(dat11) <- ~MLGRegion
pal <- unlist(map2(LETTERS[2:4], comm, function(a, b) setNames(viridis::viridis(length(b), option = a), b)))
pal <- c(pal, setNames(rep("#F6F6F6FF", 3), c("Mexico", "ID", "WI")))
```


## MLG table

First, it would be nice to visualize the MLGs across populations


```r
mtab <- mlg.table(dat11, color = TRUE)
```

![plot of chunk mlg_table](./figures/mlg-mcg///mlg_table-1.png)

```r
mplot <- last_plot()
mplot + scale_fill_viridis(discrete = TRUE, direction = -1, option = "C") +
  aes(color = I("black")) 
```

![plot of chunk mlg_table](./figures/mlg-mcg///mlg_table-2.png)

Now we can take a look at the concordance of MLGs to MCGs



```r
mll.custom(dat11) <- strata(dat11)$MCG
mcgmlg <- as.data.frame(table(mll(dat11, "original"), mll(dat11, "custom"))) %>%
  setNames(c("MLG", "MCG", "Freq")) %>%
  as_data_frame() %>%
  filter(Freq > 0)
mcgs <- mcgmlg %>%
  group_by(MCG) %>%
  summarize(MLGs = sum(Freq > 0), Samples = sum(Freq), Entropy = vegan::diversity(Freq), data = list(data_frame(MLG = MLG, Freq = Freq))) %>%
  arrange(desc(MLGs))
mlgs <- mcgmlg %>%
  group_by(MLG) %>%
  summarize(MCGs = sum(Freq > 0), Samples = sum(Freq), Entropy = vegan::diversity(Freq), data = list(data_frame(MCG = MCG, Freq = Freq))) %>%
  arrange(desc(Samples), desc(MCGs))
mcgs
```

```
## # A tibble: 87 × 5
##       MCG  MLGs Samples  Entropy              data
##    <fctr> <int>   <int>    <dbl>            <list>
## 1       5    37      73 3.050136 <tibble [37 × 2]>
## 2      44    19      36 2.568269 <tibble [19 × 2]>
## 3       1    10      15 2.153532 <tibble [10 × 2]>
## 4       4     9      14 1.965237  <tibble [9 × 2]>
## 5       2     9      10 2.163956  <tibble [9 × 2]>
## 6      53     9       9 2.197225  <tibble [9 × 2]>
## 7       3     8       8 2.079442  <tibble [8 × 2]>
## 8       9     8      15 1.599015  <tibble [8 × 2]>
## 9      45     7      16 1.559581  <tibble [7 × 2]>
## 10     16     6       7 1.747868  <tibble [6 × 2]>
## # ... with 77 more rows
```

```r
mlgs
```

```
## # A tibble: 165 × 5
##       MLG  MCGs Samples   Entropy             data
##    <fctr> <int>   <int>     <dbl>           <list>
## 1      25     5      27 0.9796882 <tibble [5 × 2]>
## 2     163     2      15 0.6909233 <tibble [2 × 2]>
## 3      65     2      11 0.3046361 <tibble [2 × 2]>
## 4     140     3      10 1.0296530 <tibble [3 × 2]>
## 5      66     1       8 0.0000000 <tibble [1 × 2]>
## 6     165     3       7 0.7963116 <tibble [3 × 2]>
## 7      78     4       6 1.2424533 <tibble [4 × 2]>
## 8     160     4       5 1.3321790 <tibble [4 × 2]>
## 9     104     2       5 0.5004024 <tibble [2 × 2]>
## 10     75     3       4 1.0397208 <tibble [3 × 2]>
## # ... with 155 more rows
```

```r
any(mcgs$MLGs > 1)
```

```
## [1] TRUE
```

We can see that there are plenty of MLGs that contain multiple MCGs. This 



```r
mlg.filter(dat11, dist = bruvo.dist, replen = other(dat11)$REPLEN) <- .Machine$double.eps + (0.5/11)
dat11
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    136 contracted multilocus genotypes
##        (0.045) [t], (bruvo.dist) [d], (farthest) [a] 
##    366 haploid individuals
##     11 codominant loci
## 
## Population information:
## 
##      6 strata - MCG, Region, Source, Year, Host, MLGRegion
##      6 populations defined - 
## Midwest, Costal, International, WI, ID, Mexico
```

```r
fmcgmlg <- as.data.frame(table(mll(dat11, "contracted"), mll(dat11, "custom"))) %>%
  setNames(c("MLG", "MCG", "Freq")) %>%
  as_data_frame() %>%
  filter(Freq > 0)
fmcgs <- fmcgmlg %>%
  group_by(MCG) %>%
  summarize(MLGs = sum(Freq > 0), Samples = sum(Freq), Entropy = vegan::diversity(Freq), data = list(data_frame(MLG = MLG, Freq = Freq))) %>%
  arrange(desc(MLGs))
fmlgs <- fmcgmlg %>%
  group_by(MLG) %>%
  summarize(MCGs = sum(Freq > 0), Samples = sum(Freq), Entropy = vegan::diversity(Freq), data = list(data_frame(MCG = MCG, Freq = Freq))) %>%
  arrange(desc(Samples), desc(MCGs))
fmcgs
```

```
## # A tibble: 87 × 5
##       MCG  MLGs Samples  Entropy              data
##    <fctr> <int>   <int>    <dbl>            <list>
## 1       5    29      73 2.789216 <tibble [29 × 2]>
## 2      44    17      36 2.405295 <tibble [17 × 2]>
## 3       1    10      15 2.153532 <tibble [10 × 2]>
## 4       4     9      14 1.965237  <tibble [9 × 2]>
## 5       2     9      10 2.163956  <tibble [9 × 2]>
## 6       3     8       8 2.079442  <tibble [8 × 2]>
## 7      53     8       9 2.043192  <tibble [8 × 2]>
## 8      45     7      16 1.559581  <tibble [7 × 2]>
## 9       9     7      15 1.506595  <tibble [7 × 2]>
## 10     16     6       7 1.747868  <tibble [6 × 2]>
## # ... with 77 more rows
```

```r
fmlgs
```

```
## # A tibble: 136 × 5
##       MLG  MCGs Samples   Entropy             data
##    <fctr> <int>   <int>     <dbl>           <list>
## 1      25     5      28 0.9569788 <tibble [5 × 2]>
## 2     163     3      16 0.8815323 <tibble [3 × 2]>
## 3     140     4      11 1.2406843 <tibble [4 × 2]>
## 4      65     2      11 0.3046361 <tibble [2 × 2]>
## 5      75     6      10 1.4978661 <tibble [6 × 2]>
## 6       9     2       8 0.3767702 <tibble [2 × 2]>
## 7      66     1       8 0.0000000 <tibble [1 × 2]>
## 8     165     3       7 0.7963116 <tibble [3 × 2]>
## 9     104     2       7 0.6829081 <tibble [2 × 2]>
## 10    152     5       6 1.5607104 <tibble [5 × 2]>
## # ... with 126 more rows
```

```r
any(fmcgs$MLGs > 1)
```

```
## [1] TRUE
```


## making a graph

I believe that making a graph to visualize this might help me understand what the h\*ck is going on. 


```r
gdf <- mcgmlg %>% 
  mutate(MLG = paste0('MLG.', MLG))
MLGS <- gdf %>% 
  group_by(MLG) %>%
  summarize(size = sum(Freq)) %>%
  rename(vertex = MLG)
MCGS <- gdf %>% 
  group_by(MCG) %>%
  summarize(size = sum(Freq)) %>%
  rename(vertex = MCG)
VAT <- bind_rows(MLGS, MCGS)
g <- gdf %>% 
  select(MCG, MLG) %>%
  graph_from_data_frame(vertices = VAT)
osize <- V(g)$size
```

Because I have more control over the size and feel of the graph, I'm going to use
ggraph. Of course, since this IS a complicated data set, It's not going to be
very pretty to look at, but I'm going to save it as supplementary materials
because it's valuable to at least look this ugliness in the face and say, "Yeah,
I guess it's not so simple after all."


```r
V(g)$size <- sqrt(osize)/10
V(g)$type <- ifelse(grepl("MLG", V(g)$name), "Multilocus Genotype", "Mycelial Compatibility Group")
set.seed(2017-05-03)
lay2 <- create_layout(g, layout = "igraph", algorithm = "nicely")


mcg_mlg_graph <- ggraph(lay2) +
  geom_node_circle(aes(r = size, fill = type)) +
  geom_node_text(aes(label = gsub("MLG.", "", name), color = type, size = size/4), show.legend = FALSE) +
  geom_edge_link(aes(start_cap = circle(node1.size, unit = "native"), 
                     end_cap = circle(node2.size, unit = "native")), 
                 arrow = arrow(length = unit(0.0125, "native"))) +
  coord_fixed() +
  scale_color_manual(values = c("white", "black")) +
  scale_fill_manual(values = c("black", "white")) +
  theme_graph(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "bottom") +
  ggtitle("Relation of Multilocus Genotypes and MCGs") 
  
mcg_mlg_graph
```

![plot of chunk unnamed-chunk-4](./figures/mlg-mcg///unnamed-chunk-4-1.png)

```r
ggsave(file.path(PROJHOME, "results/figures/publication/FigureS2.svg"), 
       width = 88*3, height = 88*3, units = "mm")
```

