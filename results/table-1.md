---
title: "Table 1"
output: 
  html_notebook:
    toc: true
editor_options: 
  chunk_output_type: inline
---


# Creation of Table 1

Table one is one of those long, descriptive tables that requires a lot of
fiddling. The *huxtable* package makes working with these tables slightly less
of a pain in the neck. Because this table involves wrapping and the text
can get split up, I am manually copying and pasting the result into my document.



```r
library("readr")
library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library("huxtable")
```

```
## 
## Attaching package: 'huxtable'
```

```
## The following object is masked from 'package:dplyr':
## 
##     add_rownames
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```r
clmn <- cols(
  .default = col_integer(),
  Region = col_character(),
  Source = col_character(),
  Host = col_character()
)

# Helper function to generate positions of blanks in a column for pretty printing
blank_it <- function(x){
  res <- rep(TRUE, length(x))
  res[1] <- FALSE
  res
}

dat <- readr::read_csv(file.path(PROJHOME, "data/clean_data.csv"), col_types = clmn) %>%
  mutate(Country = ifelse(nchar(Region) == 2, "USA", Region)) %>% # Defining country column
  mutate(Region = ifelse(nchar(Region) == 2, Region, "")) %>%     # Adding blanks for international "states"
  select(Country, Region, Source, Year, Host) %>%                 # Summarizing the number of isolates
  group_by_all() %>%                                              #
  summarize(N = n()) %>%                                          #
  arrange(desc(Country), Region, Source, Year) %>%                # Rearranging rows
  group_by(Country, Region, Source) %>%                           # Collapsing Years and Hosts
  summarize(Year = paste(sort(unique(Year)), collapse = ", "),     #
            Host = paste(sort(unique(Host)), collapse = ", "),    #
            N = sum(N)) %>%                                       #
  arrange(desc(Country), Region, Source, N) %>%                   # Rearranging rows
  group_by(Country, Region) %>%                                   # Adding blanks in repeated Country and Region Names
  mutate(set_blank = blank_it(Region)) %>%                        #
  ungroup() %>%                                                   #
  mutate(Country = ifelse(set_blank, "", Country)) %>%            #
  mutate(Region = ifelse(set_blank, "", Region)) %>%              #
  select(Country, Region, Source, Year, Host, N) %>%              #
  rename(State = Region) %>%                                      # Renaming columns
  rename(`Field Code` = Source)                                   #
  
dt <- huxtable::as_huxtable(dat, add_colnames = TRUE) %>% # Creating the huxtable
  set_align(everywhere, col = 4, value = "right") %>%     # Aligning numeric columns to the right
  set_align(everywhere, col = 6, value = "right") %>%     # 
  set_wrap(TRUE) %>%                                      # Wrapping columns
  set_number_format(value = 0)                            # Print numbers without decimals
col_width(dt) <- c(0.055, 0.035, 0.06, 0.1, 0.151, 0.025) # Specify width for each column
print_md(dt)                                              # print in pandoc markdown format
```

-----------------------------------------------------------------------------
Country   State  Field Code              Year Host                          N 
--------- ------ ---------- ----------------- -------------------------- ----
USA       CA     wmn               2004, 2005 Beryl, Bunsi, G122           18 

USA       CO     eat               2007, 2010 Yellow                       13 

                 gree              2007, 2010 Pinto, Yellow                14 

                 luc               2007, 2010 Yellow                       14 

                 wmn                     2003 GH                            1 

USA       ID     unk                     2003 GH                            1 

USA       MI     eln                     2009 unk                           1 

                 hur               2008, 2009 BL, Black, SR06233, Vista     5 

                 mung                    2009 Zorro                         1 

                 sanc                    2008 Merlot                        3 

                 Sebewaing               2009 Black                         1 

                 tuc               2008, 2009 Black, Fuji, Vista            6 

                 unk               2003, 2009 GH, unk                       2 

                 wmn        2003, 2004, 2005, 11A, 37, 38, B07104, Beryl   43 
                                   2008, 2009 , Bunsi, cornell, G122, Or      
                                              ion, PO7863, WM31               

USA       MN     wmn               2003, 2004 Beryl, Bunsi, G122           11 

USA       ND     gfc               2007, 2010 unk                           9 

                 nec                     2010 unk                           1 

                 pmc                     2010 unk                           7 

                 rrv                     2007 unk                          21 

                 stc                     2010 unk                           2 

                 tlc                     2010 unk                           9 

                 wlc                     2010 unk                           4 

                 wmn                     2005 Beryl, Bunsi, G122            7 

USA       NE     mtc               2007, 2010 Beryl, Emerson, Orion        12 

                 sbf               2007, 2009 Beryl, Pinto, Weihing         6 

                 unk                     2003 GH                            2 

                 wmn        2004, 2005, 2008, Beryl, Bunsi, G122, PO7683   27 
                                         2010 , unk                           

USA       NY     unk                     2003 GH                            1 

USA       OR     corv                    2003 G122                          1 

                 unk                     2003 GH                            1 

                 wmn               2003, 2004 Beryl, Bunsi, G122           15 

USA       WA     prsr                    2007 Merlot, Pinto, redkid        22 

                 unk                     2003 GH                            1 

                 wmn        2003, 2004, 2005, 11A, 37, 38, Beryl, Bunsi,   36 
                                         2008  cornell, G122, Orion, PO7      
                                              104, PO7863, WM31               

USA       WI     dfor                    2003 GH                            1 

                 unk                     2003 GH                            1 

Mexico           wmn                     2005 Beryl, Bunsi, G122           18 

France           flds                    2012 unk                           4 

                 wmn               2004, 2005 Beryl, Bunsi, G122           18 

Australia        vic                     2004 Beryl                         2 

                 wmn                     2004 Beryl, Bunsi, G122            4 

-----------------------------------------------------------------------------



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
##  compiler      3.4.0   2017-04-21 local          
##  datasets    * 3.4.0   2017-04-21 local          
##  devtools      1.13.2  2017-06-02 CRAN (R 3.4.0) 
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.0) 
##  dplyr       * 0.7.1   2017-06-22 CRAN (R 3.4.0) 
##  evaluate      0.10    2016-10-11 CRAN (R 3.4.0) 
##  ezknitr       0.6     2016-09-16 CRAN (R 3.4.0) 
##  glue          1.1.1   2017-06-21 CRAN (R 3.4.0) 
##  graphics    * 3.4.0   2017-04-21 local          
##  grDevices   * 3.4.0   2017-04-21 local          
##  hms           0.3     2016-11-22 CRAN (R 3.4.0) 
##  huxtable    * 0.3.0   2017-05-18 CRAN (R 3.4.0) 
##  knitr       * 1.16    2017-05-18 CRAN (R 3.4.0) 
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0) 
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0) 
##  methods     * 3.4.0   2017-04-21 local          
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0) 
##  R.methodsS3   1.7.1   2016-02-16 CRAN (R 3.4.0) 
##  R.oo          1.21.0  2016-11-01 CRAN (R 3.4.0) 
##  R.utils       2.5.0   2016-11-07 CRAN (R 3.4.0) 
##  R6            2.2.2   2017-06-17 cran (@2.2.2)  
##  Rcpp          0.12.11 2017-05-22 cran (@0.12.11)
##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.0) 
##  rlang         0.1.1   2017-05-18 CRAN (R 3.4.0) 
##  stats       * 3.4.0   2017-04-21 local          
##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.0) 
##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.0) 
##  tibble        1.3.3   2017-05-28 CRAN (R 3.4.0) 
##  tools         3.4.0   2017-04-21 local          
##  utils       * 3.4.0   2017-04-21 local          
##  withr         1.0.2   2016-06-20 CRAN (R 3.4.0)
```
