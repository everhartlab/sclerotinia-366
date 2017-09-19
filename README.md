# Analysis of 366 *S. sclerotiorum* isolates

This is an attempt to clean up an analysis of 366 isolates of *Sclerotinia
sclerotiorum* from the US and various countries around the world. 

The manuscript has mostly been written, but it still needs some cleaning:

> Population structure and phenotypic variation of *Sclerotinia sclerotiorum*
> from dry bean in the United States
>
> Z. N. Kamvar, B. S. Amaradasa, R. Jhala, S. McCoy, J. Steadman, and S. E. Everhart

# TOC

The analyses are arranged in the following order according to the [Makefile]:

 1. [table-1.md]
 1. [MCG-virulence.md]
 1. [locus-stats.md]
 1. [MLG-distribution.md]
 1. [mlg-mcg.md]
 1. [RDA-analysis.md]
 1. [pop-diff.md]
 1. [tree.md]
 1. [wmn-differentiation.md]
 1. [by-year.md]

# Reproducing the analysis

## Locally

This project is controlled via a [Makefile] which means that everything (analyses, tables, figures, the paper itself) is controlled via one command:

```
make
```

This will bootstrap the installation (warning: it will update packages),
process the data, perform the analyses, and compile the paper. 

Required software:

 - GNU Make (If you're on Windows, you can use MinGW: http://www.mingw.org/)
 - [R (version 3.4.1 or greater)](https://r-project.org)
 - [LaTeX](https://www.latex-project.org/get)
 - [pandoc](http://pandoc.org/) (Note: pandoc ships with Rstudio)
 - [devtools](https://github.com/hadley/devtools#readme)


 ## Docker

This repository contains a [Dockerfile](Dockerfile), which specifies the
instructions to build a [docker](https://www.docker.com/) container. This is
designed to capture the complete development environment of the analysis so
that it can be accurately reproduced. The image is ~2.71Gb, so be sure that you
have enough memory on your computer to run it. 

To Install Docker, go here: https://docs.docker.com/engine/installation/#desktop

Once you've installed it, open docker. This should open a command line prompt.
Navigate to this directory and type:

```sh
docker build .
```

This will bootstrap the environment to build the container. Once it's finished
you should see: 

```
Successfully built 8e1e4cd82e19
```

where `8e1e4cd82e19` will be replaced by your hash. Now that things are built,
you can run the analysis in the image with:

```
docker run -it -v $(pwd):/analysis 8e1e4cd82e19 bash
root@a4d84d719553:/# cd analysis/
root@a4d84d719553:/analysis# make clean
root@a4d84d719553:/analysis# make
```


[Makefile]: Makefile
[table-1.md]: results/table-1.md
[MCG-virulence.md]: results/MCG-virulence.md
[locus-stats.md]: results/locus-stats.md
[MLG-distribution.md]: results/MLG-distribution.md
[mlg-mcg.md]: results/mlg-mcg.md
[RDA-analysis.md]: results/RDA-analysis.md
[pop-diff.md]: results/pop-diff.md
[tree.md]: results/tree.md
[wmn-differentiation.md]: results/wmn-differentiation.md
[by-year.md]: results/by-year.md
