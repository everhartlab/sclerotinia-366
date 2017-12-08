---
title: Re-analysis of *Sclerotinia sclerotiorum* data
author: Zhian N. Kamvar
output: 
  html_notebook:
    toc: true
editor_options: 
  chunk_output_type: inline
---





# Introduction

In the PLoS paper "Genetic Variation of *Sclerotinia sclerotiorum* from
Multiple Crops in the North Central United States" (DOI:
10.1371/journal.pone.0139188), Aldrich-Wolfe *et al.* show that there is more
genetic differentiation within host than among hosts in *S. sclerotiorum*.

It is of interest here to determine if the data from the present study is similar
to the data in their study. 


# Packages

I'm going to use *poppr* to analyze the data, *dplyr* to merge the data, and
*tidyr* to unnest individuals.



```r
library("poppr")
library("tidyverse")
```


# Data Preparation

Here, I'm going to re-analyze their data. They did not provide an easily
accessible format, so I was forced to copy the data from the paper.

> A note on cleaning, The following data was corrected from soybean to sunflower:
> 
> ```
> 3;3a;2008;North Dakota;Helianthus annuus;soybean;4;129,134,180,189
> ```


```r
haplotypes <- "Haplotype	7_2	8_3	110_4	55_4	13_2	23_4	7_3	5_2	17_3	12_2	92_4	106_4
1	188	268	382	189	311	407	225	332	359	238	388	572
2a	188	270	397	189	322	407	220	332	359	232	388	572
2b	188	270	397	189	322	407	220	332	359	232	388	576
3a	186	268	390	189	332	407	225	332	353	238	388	577
3b	186	268	394	189	332	407	225	332	353	238	388	577
4	186	268	382	189	311	407	222	332	353	234	388	572
6a	186	268	397	189	332	407	222	332	353	234	388	568
6b	186	268	397	189	332	407	222	332	353	234	388	593
7	188	270	390	189	322	407	225	332	353	238	390	577
8a	188	270	390	189	322	407	220	332	359	232	390	581
8b	188	270	390	189	322	407	220	332	359	232	390	585
8c	188	270	390	189	322	407	220	332	365	232	390	605
9	188	268	390	189	332	407	225	334	359	238	391	568
11	188	268	393	189	332	407	220	332	359	232	391	554
12a	176	262	386	205	332	405	220	332	365	232	388	574
12b	176	262	386	205	332	405	220	332	365	232	388	578
13	188	272	382	193	322	407	222	332	359	234	388	565
14	188	272	393	189	332	407	225	332	361	238	390	550
15	186	272	382	173	322	407	222	332	359	234	388	554
16	186	270	397	189	332	407	222	332	353	234	388	585
17	188	270	397	189	332	407	222	332	353	234	388	574
18	186	268	382	240	332	407	225	332	359	238	388	585
19	186	270	397	181	332	407	220	334	359	232	391	550
20	188	270	382	181	322	407	222	332	353	234	391	543
21	186	270	397	189	322	407	220	332	359	232	388	572
23	186	268	397	173	311	405	220	332	359	232	388	589
24	188	270	397	189	322	407	220	332	359	232	388	572
25a	176	262	386	232	332	405	220	332	375	232	388	598
25b	176	262	386	232	332	405	220	332	375	232	388	589
28a	186	270	390	189	338	407	225	332	359	238	388	572
28b	186	270	390	189	342	407	225	332	359	238	388	572
30	186	270	382	181	322	407	225	332	353	238	388	562
31	188	268	390	181	322	407	225	334	353	238	388	572
33	188	270	390	193	322	407	220	332	359	232	390	585
34	186	268	397	177	347	407	225	332	359	238	388	554
36	186	270	394	193	322	407	222	332	353	234	388	558
37	188	268	382	193	311	407	225	332	361	238	390	572
38	186	268	382	181	311	407	222	332	359	234	390	581
43	186	270	390	185	332	407	222	332	353	234	391	585
46	188	268	382	193	322	407	225	332	359	238	391	546
53	188	268	397	181	374	407	225	332	359	238	388	572
57	176	262	386	185	332	405	220	332	367	232	388	535
60	188	268	397	181	332	407	220	334	359	232	391	568
61	188	268	397	181	332	407	225	332	359	238	388	554
62	186	268	382	173	311	405	220	332	371	232	390	577
64	186	270	397	185	332	407	222	332	353	234	388	562
66	188	270	390	189	332	407	220	332	353	232	390	581
69	188	270	382	189	311	407	225	332	378	238	395	574
71	176	262	386	228	332	405	220	332	367	232	388	535
73	186	268	382	181	311	407	220	332	353	232	390	581
77	188	268	397	181	322	407	222	332	353	234	388	565
"

haplotypes <- readr::read_table2(haplotypes)

replens <- c(
"5_2" = 2,
"7_2" = 2,
"7_3" = 2,
"8_3" = 2,
"12_2" = 2,
"13_2" = 5,
"17_3" = 3,
"23_4" = 2,
"55_4" = 4,
"92_4" = 2,
"106_4" = 4,
"110_4" = 4
)

isolates <- "MCG;Haplotype;Year;State/Province;Host;Host common name;Number of isolates;Isolate designation
1;1;2007;North Dakota;Brassica napus;canola;1;100
2;2a;2008;Minnesota;Phaseolus vulgaris;dry bean;1;137
2;2a;2008;Nebraska;Glycine max;soybean;1;191
2;2a;2008;North Dakota;Brassica napus;canola;2;175,181
2;2a;2008;North Dakota;Phaseolus vulgaris;dry bean;2;101,122
2;2a;2008;North Dakota;Glycine max;soybean;1;132
2;2a;2008;South Dakota;Glycine max;soybean;1;201
2;2a;2007;North Dakota;Brassica napus;canola;1;800
2;2a;1996;Colorado;Phaseolus vulgaris;dry bean;1;196
2;2a;2008;Manitoba, Canada;Phaseolus vulgaris;dry bean;1;103
3;3a;2008;North Dakota;Brassica napus;canola;1;159
3;3a;2008;North Dakota;Phaseolus vulgaris;dry bean;1;102
3;3a;2007;North Dakota;Brassica napus;canola;2;700,900
3;3a;-;Minnesota;Daucus carota;carrot;1;205
3;3a;2008;North Dakota;Helianthus annuus;sunflower;4;129,134,180,189
3;3b;2008;North Dakota;Phaseolus vulgaris;dry bean;1;156
3;3b;2008;North Dakota;Glycine max;soybean;1;240
3;3b;2008;North Dakota;Helianthus annuus;sunflower;1;188
4;4;2008;North Dakota;Phaseolus vulgaris;dry bean;1;104
5;2a;2008;North Dakota;Phaseolus vulgaris;dry bean;1;105
5;2a;2008;North Dakota;Glycine max;soybean;1;116
6;6a;2008;Minnesota;Helianthus annuus;sunflower;2;136,139
6;6a;2008;North Dakota;Phaseolus vulgaris;dry bean;1;106
6;6b;2008;North Dakota;Glycine max;soybean;1;176
7;7;2008;North Dakota;Glycine max;soybean;1;107
8;8a;2008;Iowa;Glycine max;soybean;1;213
8;8a;2008;North Dakota;Brassica napus;canola;2;167,186
8;8a;2008;North Dakota;Phaseolus vulgaris;dry bean;4;109,115,119,133
8;8b;2008;North Dakota;Phaseolus vulgaris;dry bean;1;120
8;8b;2007;Montana;Carthamus tinctorius;safflower;1;244
8;8b;2006;Minnesota;Glycine max;soybean;1;203
8;8c;2008;North Dakota;Glycine max;soybean;1;239
8;8c;2008;North Dakota;Helianthus annuus;sunflower;1;183
9;9;2008;Illinois;Glycine max;soybean;1;148
9;9;2008;Indiana;Glycine max;soybean;2;206,207
9;9;2008;Iowa;Glycine max;soybean;6;211,212,214,215,216,217
9;9;2008;Minnesota;Glycine max;soybean;4;145,146,197,204
9;9;2008;Minnesota;Helianthus annuus;sunflower;1;141
9;9;2008;Nebraska;Glycine max;soybean;1;153
9;9;2008;North Dakota;Brassica napus;canola;1;150
9;9;2008;North Dakota;Phaseolus vulgaris;dry bean;7;108,117,124,131,151,157,184
9;9;2008;North Dakota;Glycine max;soybean;5;113,114,125,135,178
9;9;2008;North Dakota;Helianthus annuus;sunflower;2;126,128
9;9;2007;Illinois;Glycine max;soybean;1;147
9;9;2007;North Dakota;Brassica napus;canola;1;110
9;9;2004;Missouri;Glycine max;soybean;1;234
9;9;2004;Wisconsin;Nicotiana tabacum;tobacco;1;233
9;9;2003;Wisconsin;Nicotiana tabacum;tobacco;1;232
9;9;2002;Ohio;Glycine max;soybean;1;210
9;9;2002;Wisconsin;Glycine max;soybean;1;226
9;9;2000;Wisconsin;Glycine max;soybean;4;220,221,223,224
11;11;1987;Colorado;Solanum tuberosum;potato;1;198
12;12a;2008;Minnesota;Glycine max;soybean;1;140
12;12a;2008;North Dakota;Helianthus annuus;sunflower;1;123
12;12b;2008;North Dakota;Phaseolus vulgaris;dry bean;1;118
12;12b;2008;North Dakota;Glycine max;soybean;1;127
12;12b;2006;Colorado;Helianthus annuus;sunflower;1;193
13;13;2008;North Dakota;Brassica napus;canola;1;168
14;14;2007;North Dakota;Brassica napus;canola;1;111
15;15;2008;Minnesota;Brassica napus;canola;1;190
15;15;2008;North Dakota;Brassica napus;canola;2;171,177
15;15;2007;North Dakota;Brassica napus;canola;1;112
15;15;1992;Colorado;Brassica napus;canola;1;192
16;16;2008;North Dakota;Glycine max;soybean;1;130
17;17;2008;North Dakota;Brassica napus;canola;1;185
18;18;2008;Minnesota;Phaseolus vulgaris;dry bean;1;138
19;19;2008;Minnesota;Glycine max;soybean;1;142
19;19;2003;Wisconsin;Glycine max;soybean;1;229
20;20;2008;Minnesota;Helianthus annuus;sunflower;1;143
21;21;2008;Minnesota;Glycine max;soybean;1;144
23;23;2008;North Dakota;Helianthus annuus;sunflower;1;149
24;24;2008;Nebraska;Phaseolus vulgaris;dry bean;1;154
25;25a;2008;Kansas;Helianthus annuus;sunflower;1;246
25;25a;2008;North Dakota;Phaseolus vulgaris;dry bean;1;121
25;25a;2008;North Dakota;Helianthus annuus;sunflower;1;155
25;25b;2008;North Dakota;Glycine max;soybean;1;241
28;28a;2008;North Dakota;Phaseolus vulgaris;dry bean;1;158
28;28a;2008;North Dakota;Helianthus annuus;sunflower;2;160,161
28;28a;2007;North Dakota;Brassica napus;canola;1;500
28;28b;2008;North Dakota;Helianthus annuus;sunflower;1;166
30;30;2008;North Dakota;Helianthus annuus;sunflower;1;163
31;31;2008;North Dakota;Brassica napus;canola;1;165
33;33;2008;North Dakota;Phaseolus vulgaris;dry bean;1;169
34;34;2008;North Dakota;Brassica napus;canola;1;170
36;36;2008;North Dakota;Phaseolus vulgaris;dry bean;1;172
37;37;2008;North Dakota;Brassica napus;canola;1;173
38;38;2008;North Dakota;Helianthus annuus;sunflower;2;174,182
43;43;2008;North Dakota;Glycine max;soybean;1;179
46;46;2008;North Dakota;Helianthus anuus;sunflower;1;187
49;12a;2008;Nebraska;Phaseolus vulgaris;dry bean;1;195
52;23;2008;Wyoming;Phaseolus vulgaris;dry bean;1;200
53;53;2008;South Dakota;Helianthus annuus;sunflower;1;202
56;12a;1997;Iowa;Glycine max;soybean;1;208
57;57;1998;Illinois;Glycine max;soybean;1;209
60;60;2008;Michigan;Glycine max;soybean;2;218,219
61;61;2002;North Dakota;Brassica napus;canola;1;230
61;61;2000;Wisconsin;Glycine max;soybean;1;222
62;62;2003;Wisconsin;Glycine max;soybean;2;227,228
63;24;2008;Wyoming;Phaseolus vulgaris;dry bean;1;199
64;64;2002;Wisconsin;Glycine max;soybean;1;225
66;66;2002;North Dakota;Brassica napus;canola;1;231
69;69;2005;Montana;Cynoglossum officinale;houndstonque;1;243
71;71;2008;Kansas;Phaseolus vulgaris;dry bean;1;247
73;73;2007;North Dakota;Brassica napus;canola;1;600
77;77;2008;North Dakota;Brassica napus;canola;1;162
78;19;2008;North Dakota;Brassica napus;canola;1;164"

isocol <- readr::cols(
  MCG = col_integer(),
  Haplotype = col_character(),
  Year = col_character(),
  `State/Province` = col_character(),
  Host = col_character(),
  `Host common name` = col_character(),
  `Number of isolates` = col_integer(),
  `Isolate designation` = col_character()
)

isolates <- readr::read_delim(isolates, ";", col_types = isocol)
```

Now that I have the data read in to data frames, I can create a master data
frame by joining them together.


```r
ssc <- left_join(isolates, haplotypes, by = "Haplotype") %>%
	mutate(`Isolate designation` = strsplit(`Isolate designation`, ",")) %>%
	tbl_df() %>%
	unnest() %>%
	select(`Isolate designation`, everything())
ssc
```

```
## # A tibble: 145 x 20
##    `Isolate designation`   MCG Haplotype  Year `State/Province`
##                    <chr> <int>     <chr> <chr>            <chr>
##  1                   100     1         1  2007     North Dakota
##  2                   137     2        2a  2008        Minnesota
##  3                   191     2        2a  2008         Nebraska
##  4                   175     2        2a  2008     North Dakota
##  5                   181     2        2a  2008     North Dakota
##  6                   101     2        2a  2008     North Dakota
##  7                   122     2        2a  2008     North Dakota
##  8                   132     2        2a  2008     North Dakota
##  9                   201     2        2a  2008     South Dakota
## 10                   800     2        2a  2007     North Dakota
## # ... with 135 more rows, and 15 more variables: Host <chr>, `Host common
## #   name` <chr>, `Number of isolates` <int>, `7_2` <int>, `8_3` <int>,
## #   `110_4` <int>, `55_4` <int>, `13_2` <int>, `23_4` <int>, `7_3` <int>,
## #   `5_2` <int>, `17_3` <int>, `12_2` <int>, `92_4` <int>, `106_4` <int>
```

Now we can read the data into poppr


```r
gid <- df2genind(ssc[-c(1:8)], ind.names = ssc$`Isolate designation`,
                 strata = ssc[2:8], ploidy = 1, pop = ssc$`Host common name`) %>%
  as.genclone()
gid
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##     49 original multilocus genotypes 
##    145 haploid individuals
##     12 codominant loci
## 
## Population information:
## 
##      7 strata - 
## MCG, Haplotype, Year, State.Province, Host, Host.common.name, Number.of.isolates
##      9 populations defined - 
## canola, dry bean, soybean, ..., tobacco, potato, houndstonque
```

# Comparison

First, I need to load our data:


```r
load(file.path(PROJHOME, "data/sclerotinia_16_loci.rda"))
```

Now, I need to make sure that the locus names between their data and our data are cromulent:


```r
locNames(dat11)
```

```
##  [1] "5-2(F)"   "6-2(F)"   "7-2(F)"   "8-3(H)"   "9-2(F)"   "12-2(H)" 
##  [7] "17-3(H)"  "20-3(F)"  "55-4(F)"  "110-4(H)" "114-4(H)"
```

```r
locNames(gid)
```

```
##  [1] "7_2"   "8_3"   "110_4" "55_4"  "13_2"  "23_4"  "7_3"   "5_2"  
##  [9] "17_3"  "12_2"  "92_4"  "106_4"
```

They are not, so I will have to massage them a bit to make sure that we can match up the loci and then subset them to only the common loci.


```r
ours <- dat11
locNames(gid) <- gsub("_", "-", locNames(gid))
locNames(ours) <- gsub("[(FH)]", "", locNames(ours))
locNames(gid)
```

```
##  [1] "7-2"   "8-3"   "110-4" "55-4"  "13-2"  "23-4"  "7-3"   "5-2"  
##  [9] "17-3"  "12-2"  "92-4"  "106-4"
```

```r
locNames(ours)
```

```
##  [1] "5-2"   "6-2"   "7-2"   "8-3"   "9-2"   "12-2"  "17-3"  "20-3" 
##  [9] "55-4"  "110-4" "114-4"
```

```r
(AW <- gid[loc = locNames(gid) %in% locNames(ours), mlg.reset = TRUE])
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##     43 original multilocus genotypes 
##    145 haploid individuals
##      7 codominant loci
## 
## Population information:
## 
##      7 strata - 
## MCG, Haplotype, Year, State.Province, Host, Host.common.name, Number.of.isolates
##      9 populations defined - 
## canola, dry bean, soybean, ..., tobacco, potato, houndstonque
```

```r
(ZK <- ours[loc = locNames(ours) %in% locNames(gid), mlg.reset = TRUE])
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    140 original multilocus genotypes 
##    366 haploid individuals
##      7 codominant loci
## 
## Population information:
## 
##      5 strata - MCG, Region, Source, Year, Host
##     14 populations defined - NE, NY, MN, ..., France, Mexico, ND
```

So, out of the 11 loci in this study, only 7 match up with Aldrich-Wolfe 2015. This reduces us to 140 MLGs. 
One quick way of assessing whether or not there are any shared genotypes in this data set is to look at the
allele sizes. 



```r
sort_alleles <- . %>% as.integer() %>% sort()
allist <- list(
  AW = map(alleles(AW), { . %>% sort_alleles } ),
  ZK = map(alleles(ZK), { . %>% sort_alleles } )
)
purrr::transpose(allist)
```

```
## $`7-2`
## $`7-2`$AW
## [1] 176 186 188
## 
## $`7-2`$ZK
## [1] 158 160 162 168 170 172 174
## 
## 
## $`8-3`
## $`8-3`$AW
## [1] 262 268 270 272
## 
## $`8-3`$ZK
## [1] 244 248 250 252 254 256 270
## 
## 
## $`110-4`
## $`110-4`$AW
## [1] 382 386 390 393 394 397
## 
## $`110-4`$ZK
## [1] 370 374 378 382 386
## 
## 
## $`55-4`
## $`55-4`$AW
##  [1] 173 177 181 185 189 193 205 228 232 240
## 
## $`55-4`$ZK
##  [1] 153 157 165 169 173 177 188 204 212 216
## 
## 
## $`5-2`
## $`5-2`$AW
## [1] 332 334
## 
## $`5-2`$ZK
## [1] 318 320 322 324
## 
## 
## $`17-3`
## $`17-3`$AW
## [1] 353 359 361 365 367 371 375 378
## 
## $`17-3`$ZK
## [1] 342 348 351 354 357 360 363
## 
## 
## $`12-2`
## $`12-2`$AW
## [1] 232 234 238
## 
## $`12-2`$ZK
## [1] 214 216 218 220 222
```

```r
# Do any of the alleles match?
(matching_alleles <- map_lgl(purrr::transpose(allist), ~{ any(.$AW %in% .$ZK) }))
```

```
##   7-2   8-3 110-4  55-4   5-2  17-3  12-2 
## FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE
```

Out of this, we only have 3 matches, which means that no genotypes are shared. 


<details>
<summary>Session Information</summary>


```
## Session info --------------------------------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.2 (2017-09-28)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       UTC                         
##  date     2017-12-07
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source         
##  ade4        * 1.7-8   2017-08-09 cran (@1.7-8)  
##  adegenet    * 2.1.0   2017-10-12 cran (@2.1.0)  
##  ape           5.0     2017-10-30 cran (@5.0)    
##  assertr       2.0.2.2 2017-06-06 cran (@2.0.2.2)
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.2) 
##  base        * 3.4.2   2017-12-01 local          
##  bindr         0.1     2016-11-13 CRAN (R 3.4.2) 
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.2) 
##  boot          1.3-20  2017-07-30 cran (@1.3-20) 
##  broom         0.4.3   2017-11-20 CRAN (R 3.4.2) 
##  cellranger    1.1.0   2016-07-27 CRAN (R 3.4.2) 
##  cli           1.0.0   2017-11-05 CRAN (R 3.4.2) 
##  cluster       2.0.6   2017-03-16 CRAN (R 3.4.2) 
##  coda          0.19-1  2016-12-08 cran (@0.19-1) 
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.2) 
##  compiler      3.4.2   2017-12-01 local          
##  crayon        1.3.4   2017-09-16 CRAN (R 3.4.2) 
##  datasets    * 3.4.2   2017-12-01 local          
##  deldir        0.1-14  2017-04-22 cran (@0.1-14) 
##  devtools      1.13.4  2017-11-09 CRAN (R 3.4.2) 
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.2) 
##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2) 
##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.2) 
##  expm          0.999-2 2017-03-29 cran (@0.999-2)
##  ezknitr       0.6     2016-09-16 cran (@0.6)    
##  fastmatch     1.1-0   2017-01-28 cran (@1.1-0)  
##  forcats     * 0.2.0   2017-01-23 CRAN (R 3.4.2) 
##  foreign       0.8-69  2017-06-21 CRAN (R 3.4.2) 
##  gdata         2.18.0  2017-06-06 cran (@2.18.0) 
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.4.2) 
##  glue          1.2.0   2017-10-29 CRAN (R 3.4.2) 
##  gmodels       2.16.2  2015-07-22 cran (@2.16.2) 
##  graphics    * 3.4.2   2017-12-01 local          
##  grDevices   * 3.4.2   2017-12-01 local          
##  grid          3.4.2   2017-12-01 local          
##  gtable        0.2.0   2016-02-26 CRAN (R 3.4.2) 
##  gtools        3.5.0   2015-05-29 cran (@3.5.0)  
##  haven         1.1.0   2017-07-09 CRAN (R 3.4.2) 
##  hms           0.4.0   2017-11-23 CRAN (R 3.4.2) 
##  htmltools     0.3.6   2017-04-28 CRAN (R 3.4.2) 
##  httpuv        1.3.5   2017-07-04 CRAN (R 3.4.2) 
##  httr          1.3.1   2017-08-20 CRAN (R 3.4.2) 
##  igraph        1.1.2   2017-07-21 CRAN (R 3.4.2) 
##  jsonlite      1.5     2017-06-01 CRAN (R 3.4.2) 
##  knitr       * 1.17    2017-08-10 CRAN (R 3.4.2) 
##  lattice       0.20-35 2017-03-25 CRAN (R 3.4.2) 
##  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2) 
##  LearnBayes    2.15    2014-05-29 cran (@2.15)   
##  lubridate     1.7.1   2017-11-03 CRAN (R 3.4.2) 
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.2) 
##  MASS          7.3-47  2017-04-21 CRAN (R 3.4.2) 
##  Matrix        1.2-12  2017-11-16 CRAN (R 3.4.2) 
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.2) 
##  methods     * 3.4.2   2017-12-01 local          
##  mgcv          1.8-22  2017-09-19 CRAN (R 3.4.2) 
##  mime          0.5     2016-07-07 CRAN (R 3.4.2) 
##  mnormt        1.5-5   2016-10-15 CRAN (R 3.4.2) 
##  modelr        0.1.1   2017-07-24 CRAN (R 3.4.2) 
##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.2) 
##  nlme          3.1-131 2017-02-06 CRAN (R 3.4.2) 
##  parallel      3.4.2   2017-12-01 local          
##  pegas         0.10    2017-05-03 cran (@0.10)   
##  permute       0.9-4   2016-09-09 cran (@0.9-4)  
##  phangorn      2.3.1   2017-11-01 cran (@2.3.1)  
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.2) 
##  plyr          1.8.4   2016-06-08 CRAN (R 3.4.2) 
##  poppr       * 2.5.0   2017-09-11 cran (@2.5.0)  
##  psych         1.7.8   2017-09-09 CRAN (R 3.4.2) 
##  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2) 
##  quadprog      1.5-5   2013-04-17 cran (@1.5-5)  
##  R.methodsS3   1.7.1   2016-02-16 cran (@1.7.1)  
##  R.oo          1.21.0  2016-11-01 cran (@1.21.0) 
##  R.utils       2.6.0   2017-11-05 cran (@2.6.0)  
##  R6            2.2.2   2017-06-17 CRAN (R 3.4.2) 
##  Rcpp          0.12.14 2017-11-23 CRAN (R 3.4.2) 
##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.2) 
##  readxl        1.0.0   2017-04-18 CRAN (R 3.4.2) 
##  reshape2      1.4.2   2016-10-22 CRAN (R 3.4.2) 
##  rlang         0.1.4   2017-11-05 CRAN (R 3.4.2) 
##  rstudioapi    0.7     2017-09-07 CRAN (R 3.4.2) 
##  rvest         0.3.2   2016-06-17 CRAN (R 3.4.2) 
##  scales        0.5.0   2017-08-24 CRAN (R 3.4.2) 
##  seqinr        3.4-5   2017-08-01 cran (@3.4-5)  
##  shiny         1.0.5   2017-08-23 CRAN (R 3.4.2) 
##  sp            1.2-5   2017-06-29 CRAN (R 3.4.2) 
##  spData        0.2.6.7 2017-11-28 cran (@0.2.6.7)
##  spdep         0.7-4   2017-11-22 cran (@0.7-4)  
##  splines       3.4.2   2017-12-01 local          
##  stats       * 3.4.2   2017-12-01 local          
##  stringi       1.1.6   2017-11-17 CRAN (R 3.4.2) 
##  stringr     * 1.2.0   2017-02-18 CRAN (R 3.4.2) 
##  tibble      * 1.3.4   2017-08-22 CRAN (R 3.4.2) 
##  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2) 
##  tidyverse   * 1.2.1   2017-11-14 CRAN (R 3.4.2) 
##  tools         3.4.2   2017-12-01 local          
##  utils       * 3.4.2   2017-12-01 local          
##  vegan         2.4-4   2017-08-24 cran (@2.4-4)  
##  withr         2.1.0   2017-11-01 CRAN (R 3.4.2) 
##  xml2          1.1.1   2017-01-24 CRAN (R 3.4.2) 
##  xtable        1.8-2   2016-02-05 CRAN (R 3.4.2)
```

</details>
