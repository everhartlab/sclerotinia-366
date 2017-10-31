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
library("huxtable")
```




```r
clmn <- cols(
  .default = col_integer(),
  Severity = col_double(),
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
  mutate(Region = ifelse(nchar(Region) == 2, Region, "-")) %>%    # Adding blanks for international "states"
  mutate(Source = ifelse(Source == "wmn", "wmn", "producer")) %>% # Changing fields to binary counts
  select(Country, Region, Source, Year, Host) %>%                 # Summarizing the number of isolates
  group_by_all() %>%                                              #
  summarize(N = n()) %>%                                          #
  arrange(desc(Country), Region, Source, Year) %>%                # Rearranging rows
  group_by(Country, Region, Source) %>%                           # Collapsing Years and Hosts
  summarize(Year = paste(sort(unique(Year)), collapse = ", "),    #
            Host = paste(sort(unique(Host)), collapse = ", "),    #
            N = sum(N)) %>%                                       #
  arrange(desc(Country), Region, desc(N)) %>%                     # Rearranging rows
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
  # set_wrap(TRUE) %>%                                      # Wrapping columns
  set_number_format(value = 0) %>%                        # Print numbers without decimals
  set_col_width(c(0.055, 0.035, 0.06, 0.1, 0.151, 0.025)) %>% # Specify width for each column
  print_md(max_width = 80)                                # print in pandoc markdown format
```

```
---------------------------------------------------------------------------
 Country   State Field Code              Year Host                       N 
---------- ----- ---------- ----------------- ------------------------- ---
 USA       CA    wmn               2004, 2005 Beryl, Bunsi, G122        18 
                                                                           
 USA       CO    producer          2007, 2010 Pinto, Yellow             41 
                                                                           
                 wmn                     2003 GH                         1 
                                                                           
 USA       ID    producer                2003 GH                         1 
                                                                           
 USA       MI    wmn        2003, 2004, 2005, 11A, 37, 38, B07104,      43 
                                   2008, 2009 Beryl, Bunsi, cornell,       
                                              G122, Orion, PO7863, WM31    
                                                                           
                 producer    2003, 2008, 2009 BL, Black, Fuji, GH,      19 
                                              Merlot, SR06233, unk,        
                                              Vista, Zorro                 
                                                                           
 USA       MN    wmn               2003, 2004 Beryl, Bunsi, G122        11 
                                                                           
 USA       ND    producer          2007, 2010 unk                       53 
                                                                           
                 wmn                     2005 Beryl, Bunsi, G122         7 
                                                                           
 USA       NE    wmn        2004, 2005, 2008, Beryl, Bunsi, G122,       27 
                                         2010 PO7683, unk                  
                                                                           
                 producer   2003, 2007, 2009, Beryl, Emerson, GH,       20 
                                         2010 Orion, Pinto, Weihing        
                                                                           
 USA       NY    producer                2003 GH                         1 
                                                                           
 USA       OR    wmn               2003, 2004 Beryl, Bunsi, G122        15 
                                                                           
                 producer                2003 G122, GH                   2 
                                                                           
 USA       WA    wmn        2003, 2004, 2005, 11A, 37, 38, Beryl,       36 
                                         2008 Bunsi, cornell, G122,        
                                              Orion, PO7104, PO7863,       
                                              WM31                         
                                                                           
                 producer          2003, 2007 GH, Merlot, Pinto, redkid 23 
                                                                           
 USA       WI    producer                2003 GH                         2 
                                                                           
 Mexico    -     wmn                     2005 Beryl, Bunsi, G122        18 
                                                                           
 France    -     wmn               2004, 2005 Beryl, Bunsi, G122        18 
                                                                           
                 producer                2012 unk                        4 
                                                                           
 Australia -     wmn                     2004 Beryl, Bunsi, G122         4 
                                                                           
                 producer                2004 Beryl                      2 
---------------------------------------------------------------------------
```



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
##  date     2017-10-31
```

```
## Packages ------------------------------------------------------------------------------------------
```

```
##  package     * version date       source        
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.2)
##  base        * 3.4.2   2017-10-27 local         
##  bindr         0.1     2016-11-13 CRAN (R 3.4.2)
##  bindrcpp    * 0.2     2017-06-17 CRAN (R 3.4.2)
##  compiler      3.4.2   2017-10-27 local         
##  datasets    * 3.4.2   2017-10-27 local         
##  devtools      1.13.3  2017-08-02 CRAN (R 3.4.2)
##  digest        0.6.12  2017-01-27 CRAN (R 3.4.2)
##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)
##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.2)
##  ezknitr       0.6     2016-09-16 cran (@0.6)   
##  glue          1.1.1   2017-06-21 CRAN (R 3.4.2)
##  graphics    * 3.4.2   2017-10-27 local         
##  grDevices   * 3.4.2   2017-10-27 local         
##  hms           0.3     2016-11-22 CRAN (R 3.4.2)
##  huxtable    * 1.1.0   2017-10-20 cran (@1.1.0) 
##  knitr       * 1.17    2017-08-10 CRAN (R 3.4.2)
##  magrittr      1.5     2014-11-22 CRAN (R 3.4.2)
##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.2)
##  methods     * 3.4.2   2017-10-27 local         
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.2)
##  R.methodsS3   1.7.1   2016-02-16 cran (@1.7.1) 
##  R.oo          1.21.0  2016-11-01 cran (@1.21.0)
##  R.utils       2.5.0   2016-11-07 cran (@2.5.0) 
##  R6            2.2.2   2017-06-17 CRAN (R 3.4.2)
##  Rcpp          0.12.13 2017-09-28 CRAN (R 3.4.2)
##  readr       * 1.1.1   2017-05-16 CRAN (R 3.4.2)
##  rlang         0.1.2   2017-08-09 CRAN (R 3.4.2)
##  stats       * 3.4.2   2017-10-27 local         
##  stringi       1.1.5   2017-04-07 CRAN (R 3.4.2)
##  stringr       1.2.0   2017-02-18 CRAN (R 3.4.2)
##  tibble        1.3.4   2017-08-22 CRAN (R 3.4.2)
##  tools         3.4.2   2017-10-27 local         
##  utils       * 3.4.2   2017-10-27 local         
##  withr         2.0.0   2017-07-28 CRAN (R 3.4.2)
```

</details>
