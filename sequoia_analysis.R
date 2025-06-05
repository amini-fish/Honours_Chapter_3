## install.packages("sequoia")

install.packages("sequoia", dependencies = TRUE)
install.packages("cli")

#### Set wd ####
setwd("C:/Users/samue/Desktop/Honours/analysis")

#### Load packages ####

library(sequoia) # analysis
library(cli)
library(SNPRelate)
library(dartRverse) # cleaning/data wrangling 
library(ggplot2) # sexy plots
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.base) 
library(dartR.captive) # becuase dartR is a pain
library(gplots) 
library(graph4lg) # more sexy plots
library(viridis) # why not
library(ggraph) # even more sexy plots

#### Import Data ####

gl <- get(load("C:/Users/samue/OneDrive/Desktop/Honours/analysis/daly_geno_clean.Rdata")) 

gl <- gl.filter.monomorphs(gl, v= 3)

gl # always good to have a check and make sure it's what you wanted...

gl.smearplot(gl) ## nice and polymorphic plus no na's YAY!

#### Convert to compatible matrix format ####

# Convert genlight to matrix
geno_matrix <- as.matrix(gl)  # Converts to allele counts

geno_matrix[is.na(geno_matrix)] <- -9 # Sequoia uses -9 for missing data

# Check dimensions
dim(geno_matrix)  # Should be (n_individuals, n_loci)

View(geno_matrix)
# Save to CSV if needed
write.csv(geno_matrix, "geno_data.csv", row.names = TRUE)

#### Read into Sequoia #### 

# Read in the formatted genotype matrix
sequoia_data <- as.matrix(read.csv("geno_data.csv", row.names = 1))

SnpStats(sequoia_data)

CheckGeno(sequoia_data)
#### Read in life history data ####

sawfish_metadata <- readxl::read_excel("sequoia_metadata.xlsx")

table(sawfish_metadata$BirthYear, useNA = "always")
str(sawfish_metadata)  # Check if BirthYear is numeric or character
sawfish_metadata$BirthYear <- as.numeric(sawfish_metadata$BirthYear)
unique(sawfish_metadata$BirthYear[is.na(sawfish_metadata$BirthYear)])
sawfish_metadata$BirthYear <- as.integer(sawfish_metadata$BirthYear)
sawfish_metadata$BirthYear[is.na(sawfish_metadata$BirthYear)] <- -999

sawfish_metadata$BirthYear <- as.integer(sawfish_metadata$BirthYear)
str(sawfish_metadata)
table(sawfish_metadata$BirthYear, useNA = "always")


View(sawfish_metadata)

# Run Sequoia

run1 <- sequoia(GenoM = sequoia_data, Err = 0.005,
        LifeHistData = sawfish_metadata, Module="par", Plot=FALSE, quiet=TRUE)


