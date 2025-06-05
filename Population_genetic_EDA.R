#### Set directory ####

setwd("C:/Users/samue/Desktop/Honours/analysis")

#### Install Pakcages ####
install.packages("dartRverse")
install.packages("ggplot2")
install.packages("hierfstat")
install.packages("dplyr")
install.packages("devtools")
install.packages("gplots")
install.packages("graph4lg")
install.packages("viridis")
install.packages("dartR.captive")
install.packages("dartR.base")
install.packages("dartR.sim")
installed.packages("dartR.popgen")
install.packages("ggraph")
devtools::install_version("ggplot2", "3.4.4")
install.packages("SNPRelate")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

# if needed 

devtools::install_github("green-striped-gecko/dartR.captive@dev")

#### Load Packages ####

library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.base)
library(dartR.captive)
library(gplots)
library(graph4lg)
library(viridis)
library(ggraph)
library(hierfstat)

#### LOAD IN CLEAN GENOTYPE DATA ####

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

#### EDA ####

dartR.base::gl.report.basics(gl, v = 5)

#### Calculate heterozygosity ####

dartR.base::gl.report.heterozygosity(gl,
                                     method = "pop",
                                     plot.theme = theme_bw()
                                     )

?gl.report.heterozygosity
#### Calculate allelic richness ####

?gl.report.diversity
dartR.base::gl.report.diversity(gl, 
                                plot.theme = theme_bw())


hfstat <- gl2gi(gl)

AR <- hierfstat::allelic.richness(hfstat)
print(AR)

summary(AR$Ar)
ar <- as.data.frame(AR$Ar)
mean.ar <- colMeans(ar)
sd(ar$`1`, na.rm=TRUE)/sqrt(nrow(ar) - length(which(is.na(ar$`1`)))) ## 0.002

glimpse(ar)
ggplot(data = ar, aes(x = as.factor(1),y = `1`)) +
  geom_jitter(fill = "darkolivegreen3", width = 0.4, alpha = 0.4, size = 4, shape = 21, colour = "black") + 
  scale_y_continuous(limits = c(0, 2.0)) +
  xlab("Daly") + 
  ylab("Allelic Richness") +
  theme(axis.text.x = element_blank()) +
  theme_bw()

