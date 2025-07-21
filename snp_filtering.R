setwd("C:/Users/samue/Desktop/Honours/analysis")

#Install some packages

install_github("green-striped-gecko/dartR.sexlinked@dev")
#install.packages("hierfstat")
#devtools::install_version("ggplot2", "3.4.4")

install.packages("HardyWeinberg")

## Load all required packages
library(dartR.sexlinked)
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(HardyWeinberg)

data <- "Sawfish_SNPs_genotyped.csv"
meta <- "Sawfish_meta2.csv"

## Load all genotype data and metadata

data.gl <- dartR.base::gl.read.dart(filename = data, ind.metafile = meta); data.gl

pop(data.gl) <- data.gl@other$ind.metrics$pop

table(pop(data.gl))

gl.smearplot(data.gl) 

data.gl <- gl.keep.pop(data.gl, pop.list = "Daly")

data.gl <- gl.filter.monomorphs(data.gl, v = 5)
## Save unfiltered data 

#save(data.gl, file = "pristis_geno_raw.RData") 

###### FILTERING #####

## Get rid of sex linked loci 

## Get rid of unreliable loci

gl.report.reproducibility(data.gl)

data.gl <- dartR.base::gl.filter.reproducibility(data.gl, threshold = 0.98)

## Callrate 

dartR.base::gl.report.callrate(data.gl)

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "loc", threshold = 0.95)

#Get rid of low and super high read depth loci

dartR.base::gl.report.rdepth(data.gl)

data.gl <- dartR.base::gl.filter.rdepth(data.gl, lower = 10, upper = 75)

## Remove secondary loci - important to keep unlinked loci

data.gl <- dartR.base::gl.filter.secondaries(data.gl, method = "best") 
 
#Very low filter – this is only to get rid of your really bad individuals

dartR.base::gl.report.callrate(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "ind", threshold = 0.95)

#Always run this after removing individuals – removes loci that are no longer variable

data.gl <- dartR.base::gl.filter.hwe(data.gl)

data.gl

data.gl <- gl.filter.maf(data.gl, t = 0.05, v = 5)

data.gl <- dartR.base::gl.filter.monomorphs(data.gl)


data.gl 

## remove evidence of DNA contamination ## Important for kin finding 
dartR.base::gl.report.heterozygosity(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.heterozygosity(data.gl,t.lower = 0.2,  t.upper = 0.5) ## Visualise our cleaned data 

gl.smearplot(data.gl)

####################################################

data.gl <- gl.filter.allna(data.gl)

data.gl <- gl.filter.monomorphs(data.gl)


gl.report.ld(data.gl, probar = T, nchunks = 8)

### SAVE DALY RIVER INDS ###


save(data.gl, file = "daly_geno_clean.Rdata")

## LD??



### Check our filtering steps ###

data.gl@other$history

### SAVE CLEANED SNP DATA ###

save(data.gl, file = "pristis_geno_cleaned.Rdata")


gl <- get(load("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/analysis/daly_geno_clean.Rdata")); gl




?gl2vcf()





