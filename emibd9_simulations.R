install.packages("sqldf")
install.packages("stringr")


#load everything in

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(dartR.captive)
library(dartR.sim)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.base)
library(dartR.captive)
library(gplots)
library(graph4lg)
library(adegenet)

# Set our wd to load simulated data in

setwd("C:/Users/samue/Desktop/Honours/analysis/sims2")


sim_data <- readLines("sims2.ibd9")

strt <- which(grepl("^IBD", sim_data)) + 2
stp <- which(grepl("Indiv genotypes", sim_data)) - 4
linez_headings <- sim_data[strt]
linez_data <- sim_data[(strt + 1):stp]
tmp_headings <- unlist(stringr::str_split(linez_headings, " "))
tmp_data <- stringr::str_split(linez_data, " ")

#Raw data 
tmp_data_raw_1 <- lapply(tmp_data, "[", c(2:22))
tmp_data_raw_2 <- do.call("rbind", tmp_data_raw_1)
tmp_data_raw_3 <- as.data.frame(tmp_data_raw_2)

tmp_data_raw_3

tmp_data_raw_3$V3 <- lapply(tmp_data_raw_3$V3, as.numeric)

tmp_data_raw_3


colnames(tmp_data_raw_3) <- tmp_headings[2:22]

tmp_data_raw_3

df <- data.frame(ind1=tmp_data_raw_3$Indiv1, ind2=tmp_data_raw_3$Indiv2,rel= tmp_data_raw_3$`r(1,2)`)

df$rel <- as.numeric(df$rel); df

df$index <- ifelse(df$ind1 == df$ind2, "drop", "keep")

df <- df %>% 
  filter(index == "keep")%>%
  droplevels

df <- df[,-4]

df

#select all the FS 

FS <- df[str_detect(df$ind1, "FS"), ]
FS <- data.frame(FS[str_detect(FS$ind2, "FS"), ])
FS$sibtype <- c("full-sib")

## Spread of FS 
ggplot(data = FS, aes(x = rel))+
  geom_histogram()

HS <- df[str_detect(df$ind1, "HS"), ]
HS <- data.frame(HS[str_detect(HS$ind2, "HS"), ])
HS$sibtype <- c("half-sib")


FC <- df[str_detect(df$ind1, "FC"), ]
FC <- data.frame(FC[str_detect(FC$ind2, "FC"), ])
FC$sibtype <- c("first-cousin")

UP <- df[str_detect(df$ind1, "UR"), ]
UP <- data.frame(UP[str_detect(UP$ind2, "UR"), ])
UP$sibtype <- c("unrelated")

df <- rbind(FS, HS, UP)

ggplot(df, aes(x = reorder(sibtype, -rel), y = rel, fill = sibtype))+
  geom_boxplot()+
  theme_bw()

##

# Define function for confidence interval

# Compute 95% confidence interval
ci_full_sib <- quantile(FS$rel, probs = c(0.025, 0.975))

print(ci_full_sib)

threshold_full_sib <- min(FS$rel[FS$rel > max(HS$rel)])

print(threshold_full_sib)

# Print results
print(ci_full_sib)

dens_full <- density(FS$rel)
dens_half <- density(HS$rel)

# Compute probability that a given relatedness value (e.g., R = 0.4) is from a full sibling
R_value <- 0.4
p_full <- approx(dens_full$x, dens_full$y, xout = R_value)$y
p_half <- approx(dens_half$x, dens_half$y, xout = R_value)$y

# Compute probability that R is from a full sibling
prob_full_given_R <- p_full / (p_full + p_half)

print(prob_full_given_R)

###

simplot <- ggplot(df, aes(x = rel, fill = sibtype)) +
  geom_density(alpha = 0.5, position = "identity")+
  theme_bw()+
  ggtitle("Distribution of simulated relatedness coefficients") +
  xlab("Relatedness") +
  ylab("Frequency")+
  scale_fill_brewer(palette = "Dark2")

simplot <- simplot + theme(plot.title = element_text(hjust = 0.5))

print(simplot)

#res <df2#res <- matrix(df$rel, nrow = 380, ncol = 380)

#for (i in 1:nrow(df)) {
#  res[df[i, 1], df[i, 2]] <- df[i, 3]
#}


ggplot(data = df, aes(x = rel))+
  geom_histogram(bins = 500)+
  xlim(0.05, 0.3) +
  geom_vline(xintercept = 0.125, linewidth = 1.1, col = "orange") +
  geom_vline(xintercept = 0.250, linewidth = 1.1, col = "skyblue") +
  geom_vline(xintercept = 0.092, col = "black", linewidth = 1, linetype="dotted") +
  geom_vline(xintercept = 0.158, col = "black", linewidth = 1,  linetype="dotted") +
  geom_vline(xintercept = 0.204, col = "black", linewidth = 1, linetype="dotted") +
  geom_vline(xintercept = 0.296, col = "black", linewidth = 1,  linetype="dotted")
  
  ## this is inclusive of self comparisons


gl.plot.heatmap(res)

colnames(res) <- tmp_data_raw_3$Indiv1
rownames(res) <- tmp_data_raw_3$Indiv1

df
summary(df)

devtools::install_github("green-striped-gecko/dartR.captive@dev", force = TRUE)

library(dartR.captive)


gl.sim.relatedness(gl, 
                   rel = "half.sibs", 
                   emibd9.path = "C:/EMIBD9", 
                   nboots = 2)

fathers <- dartR.sim::gl.sim.ind(gl, n = 1)
mothers <- dartR.sim::gl.sim.ind(gl, n = 1)

offspring <- gl.sim.offspring(fathers, mothers, noffpermother = 2)

all <- rbind(gl, offspring[1,])

gl.run.EMIBD9(offspring, 
              Inbreed = 1, 
              emibd9.path =  "C:/EMIBD9")

afl <- gl.alf(gl)

View(afl)

write.table(afl, file = "simdata_wide.txt", row.names = F, col.names = F)

library(tidyverse)

afl <- afl %>%
  pivot_longer(cols = starts_with("alf"),
               names_to = "Allele", values_to = "Freq") 

?sample_linked_genotype_pairs

write.table(afl, file = "simdata_long.txt", row.names = F, col.names = F)


write.table(input$freqs, file = "simdata_test.txt", row.names = F, col.names = F)

afreqs_ <- data.frame(afreqs_ready[, c(6:7)])

afreqs_

write.table(afreqs_, file = "simdata_long.txt", col.names = F, row.names = F)

### Try to build a new function to simulate 100 full-siblings 

?gl.sim.relatedness

## use the entire sawfish dataset...

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/pristis_geno_cleaned.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

source_pop <- gl.sim.ind(gl, 100)

## Simulate full sibs 

papas <- source_pop[1]
mamas <- source_pop[2]

off <- gl.sim.offspring(papas, mamas, noffpermother = 100, popname = "full-sibs", v = 5)

off@ind.names <- paste0("FS_",1:nInd(off))

## Simulate unrelated/randoms - 
#need to make sure no relatives in baseline data or large dataset

random <- gl.sim.ind(source_pop, n = 100, popname = "unrelated")

random@ind.names <- paste0("UR_",1:nInd(random))

# Merge existing data 

ppoff <- rbind(random, off)

## Simulate half siblings 

parents <- sample(1:nInd(source_pop),3, replace = F)
      
mother1 <- source_pop[parents[1],]
mother2 <- source_pop[parents[2],]
        
father <- source_pop[parents[3],]
        
hs1 <- gl.sim.offspring(father, mother1, noffpermother = 100, sexratio = 0.5, popname = "half-sibs")

hs1@ind.names <- paste0("HS_",1:nInd(hs1))
        
hs2 <- gl.sim.offspring(father, mother2, noffpermother = 100, sexratio = 0.5, popname = "half-sibs")
        
hs2@ind.names <- paste0("HS_",1:nInd(hs2))

ppoff <- rbind(ppoff,hs1[c(1:50),],hs2[51:100,])

## lets try the alternative and run it in emibd9 to speed it up..

simulate.relate <- gl.run.EMIBD9(ppoff,Inbreed = 1, emibd9.path =  "C:/EMIBD9")

diag(simulate.relate$rel) <- 0 # removes self comparisons hooray

emibd.rel <- pw_mat_to_df(simulate.relate$rel) # convert the square matrix to an edge based dataframe

df <- emibd.rel[,c(1,2,4)] # remove your third column (not needed and is confusing)

df
## Check your results: 

colnames(df) <- c("ind1", "ind2", "rel")

FS <- df[str_detect(df$ind1, "FS"), ]; FS
FS <- data.frame(FS[str_detect(FS$ind2, "FS"), ])
FS$sibtype <- c("full-sib")

## Spread of FS 

dev.off()
ggplot(data = FS, aes(x = rel))+
  geom_histogram()

HS <- df[str_detect(df$ind1, "HS"), ]
HS <- data.frame(HS[str_detect(HS$ind2, "HS"), ])
HS$sibtype <- c("half-sib")

HS <- HS %>% 
  filter(rel < 0.2)

UP <- df[str_detect(df$ind1, "UR"), ]
UP <- data.frame(UP[str_detect(UP$ind2, "UR"), ])
UP$sibtype <- c("unrelated")

df <- rbind(FS, HS, UP)

ggplot(df, aes(x = reorder(sibtype, -rel), y = rel, fill = sibtype))+
  geom_boxplot()+
  theme_bw()

## Plot the simulated values where each dyad has #1000 inds

simplot <- ggplot(df, aes(x = rel, fill = sibtype)) +
  geom_density(alpha = 0.5, position = "identity")+
  theme_bw()+
  ggtitle("Distribution of simulated relatedness coefficients") +
  xlab("Relatedness") +
  ylab("Frequency")+
  scale_fill_brewer(palette = "Dark2") +
  geom_vline(xintercept = mean(FS$rel), linetype="dotdash") +
  geom_vline(xintercept = mean(HS$rel), linetype="dotdash") +
  geom_vline(xintercept = mean(UP$rel), linetype="dotdash")

simplot <- simplot + theme(plot.title = element_text(hjust = 0.5))

print(simplot)

sd(FS$rel)

install.packages("car")
library(car)

mod <- lm(rel ~ sibtype, data = df)

set.seed(42) # for reproducibility system.time

boot <- Boot(mod, R=48000)

Confint(boot, level=.9, type="norm")

## The goal is to do the following: 
#1 - load the genotype data that was simulated using COANCESTRY into a format compatible with EMIBD9
#2 - analyse it with EMIBD9 + load back into R to compare with other estimators
if (!all(c("SampleID", "MarkerID", "Genotype") %in% names(genotype_data))) {
  stop("Data does not contain required columns.")
}


## ------------------------------------------------------------
getwd()
setwd("C:/Users/samue/Desktop/coancestry_sims")

genotype_data <- read.table("GenotypeData.txt", header = F, sep = " ", stringsAsFactors = FALSE)

ind_names <- genotype_data[,1]
loc_data <- genotype_data[,-1]

# Create the pattern for the column names
num_col <- ncol(loc_data)
num_levels <- ceiling(num_col/2)  # Adjust this to the number of Lx levels you need
num_sub_levels <- 2  # Adjust this to the number of sub-levels per level
column_names <- paste0("L", rep(1:num_levels, each = num_sub_levels), "-", rep(1:num_sub_levels, num_levels))

# Print to check
colnames(loc_data)<- column_names

genotype_data <- cbind(ind_names, loc_data)

genotype_matrix <- as.matrix(genotype_data)


View(genotype_matrix)## Convert it to a genlight

gl <- adegenet::as.genlight(genotype_matrix)

# Assign individual names (if you have them)
indNames(gl) <- rownames(genotype_data)  # Or use genotype_data$IndividualID if that was your column name

# Inspect the genlight object
print(gl)

gl.run.EMIBD9(gl, 
              Inbreed = 1, 
              emibd9.path =  "C:/EMIBD9", 
              outpath = getwd(), 
              )


gl2EMIBD9(gl, 
          Inbreed = 1, 
          ISeed = 42,
          GtypeFile = "C:/Users/samue/Desktop/coancestry_sims/simulations_4/EMIBD9_Gen.dat",
          OutFileName_par = "C:/Users/samue/Desktop/coancestry_sims/simulations_4/MyData.par",
          OutFileName = "C:/Users/samue/Desktop/coancestry_sims/simulations_4/EMIBD9_Res.ibd9")

