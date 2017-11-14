## This Dockerfile describes the R analysis environment for:
##
## Kamvar, Z.N., Amaradasa, B.S., Jhala, R., McCoy, S., Steadman, S., and
## Everhart, S.E. (2017). Population structure and phenotypic variation of
## _Sclerotinia sclerotiorum_ from dry bean in the United States PeerJ XXX
##
## package versions here are locked to those present on 2017-09-30
##
## Note: this Dockerfile was modified from
## https://github.com/NESCent/popgen-docker/blob/193387d3f1e5484ef8a1ddf6d66cfca64ccd40d7/Rpopgen/Dockerfile
##
## It also includes logic from
## https://github.com/benmarwick/mjbtramp/blob/898ee99f17d64a41161a8b6760325572c7406b4b/Dockerfile

## Note: this pulls from rocker/verse:3.4.2
FROM zkamvar/sclerotinia-366-dependencies:latest
MAINTAINER Zhian Kamvar <zkamvar@gmail.com>

## Copy the current directory to /analysis
COPY . /analysis

RUN . /etc/environment \
&& cd /analysis \
&& make clean \
&& make -j 4

# NOTES TO SELF IN A MARWICKIAN FASION
#
## To work on this project within the container
#
# docker run --rm --name ssc -dp 8787:8787 -v $(pwd):/home/rstudio/ -e ROOT=TRUE zkamvar/sclerotinia-366
#
## This will allow me to open Rstudio from my browser and work on it within the container.
## It should be noted, however that I still need to commit from my command line as the
## git bot in Rstudio doesn't appear to update inside of the container.
## 
## When done with the container, be sure to stop it:
# 
# docker stop ssc
# 

