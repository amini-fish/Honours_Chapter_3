## CRAN dependencies:

install.packages("Rcpp")

install.packages("remotes")

install.packages("tensorA")

## support packages from Mark Bravington's repo:

install.packages(pkgs = c("atease", "mvbutils", "vecless"), repos = "https://markbravington.github.io/Rmvb-repo")

## kinference and gbasics: note that you'll need to specify your own filepaths here. The files are attached.

remotes::install_local("C:/Users/samue/Desktop/Kinference/kinference-master(8).zip", subdir = "kinference")

remotes::install_local("C:/Users/samue/Desktop/Kinference/gbasics-master/gbasics")

## Okay, I think we're good to go...

library(Rcpp)
library(remotes)
library(atease)
library(mvbutils)
library(vecless)
library(kinference)

## Packages loaded now for analysis/eda

data(dropbears)

kinPalette()

glimpse(dropbears)
head(dropbears)

head(dropbears@info)
head(dropbears@locinfo)

## a boring but necessary data-prep step:
dropbears <- kin_power(dropbears, k = 0.5)

dropbears

plod_HU <- find_HSPs(dropbears, keep_thresh = -10)
PLOD_loghisto(plod_HU)
  
## Need to work out how to convert gl to snpgeno



