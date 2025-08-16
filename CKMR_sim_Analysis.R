#### Install packages ####

install.packages("vcfR")

install.packages("remotes")
remotes::install_github("eriqande/CKMRsim")

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}


install_github("https://github.com/eriqande/CKMRsim/tree/master") 

#### Load packages ####

library(devtools)
library(vcfR)
library(CKMRsim)
library(dartRverse)
library(tidyverse)
library(devEMF)
library(stats)
library(jtools)
library(bbmle)
library(car)
library(interactions)
library(lsmeans)
library(ggpubr)
library(dplyr)
library(caret)
library(GGally)
library(DHARMa)
library(pwr)
library(WebPower)
library(simr)
library(splines)
library(brglm2)
library(brms)
library(mclust)
library(MASS)     
library(ggeffects)
library(ggplot2)
library(devEMF)
library(MuMIn)
library(jtools)

#### Load the data ####

setwd("C:/Users/samue/Desktop")

gl <- get(load("C:/Users/samue/Desktop/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

gl.write.csv(gl, outfile = "new_outfile.csv",
             outpath = "C:/Users/samue/Desktop", 
             verbose = 5)

data <- read.csv("new_outfile.csv")

glimpse(data)

gl
#### Prepare the data ####

## SKIP TO LINE 107 TO START ANALYSIS ## 

gl2vcf(gl, 
      plink.bin.path = "C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/analysis/plink", 
     outfile = "outfile", 
       outpath = getwd())

geno_dummy <- read.vcfR("outfile.vcf")

geno_dummy

genotypes <- geno_dummy@gt

genotypes <- data.frame(genotypes)

genotypes

drop <- c("FORMAT")
genotypes = genotypes[,!(names(genotypes) %in% drop)]
genotypes

references <- data.frame(geno_dummy@fix)
references
SNPS <- references$ID 
SNPS <- gsub("_", "", SNPS)

genotypes$Locus <- SNPS
genotypes <- genotypes %>% relocate(Locus)

rownames(genotypes) <- SNPS

genotypes

### I want to transpose this genotype data 

genotypes_2 <- t(genotypes)

View(genotypes_2)

genotypes_3 <- data.frame(genotypes_2)

View(genotypes_3)

genotypes_3 <- genotypes_3[-1,]

genotypes_3 <- genotypes_3 %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = "_")

View(genotypes_3)

ind_names <- gl@ind.names

genotypes_3$Indiv <- ind_names

genotypes_3 <- genotypes_3 %>% relocate(Indiv)

genotypes_3 ## FUCK YES!!!

colnames(genotypes_3)

long_genos <- genotypes_3 %>%
  pivot_longer(
    cols = -Indiv, 
    names_to = c("Locus", "gene_copy"), 
    names_sep = "\\_", 
    values_to = "Allele"
  )

## Write the data to a csv file so we can just load it in as a tibble for next time and bypass all the nonesense

write.csv(long_genos, file = "Pristis_genofor_CKMRsim.csv" )

#### START HERE - Long genos ####

long_genos <- read.csv("Pristis_genofor_CKMRsim.csv"); long_genos

long_genos ## Looking MINT!!! We got there! 

long_genos <- long_genos[,-1]
locus_names <- unique(long_genos$Locus)

long_genos

afreqs_ready <- long_genos %>%
  count(Locus, Allele) %>%  
  group_by(Locus) %>%
  mutate(
    Freq = n / sum(n),
    Chrom = "Unk",
    Pos = as.integer(factor(Locus, levels = locus_names))
  ) %>%
  ungroup() %>%
  dplyr::select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele)) %>%
  reindex_markers()

afreqs_ready 

#### Create CKMR object ####

### Creating more dataframes for our analysis 

ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("PO","FS", "HS", "HAN", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

ckmr

#### Subset Kappas - IBD probs for desired relats ####

kappas # all possible

kappas[c("PO","FS", "HS", "HAN", "HFC", "U"), ]

#### Simulate dyads for each to approximate the LLRs and therefore the fpos/fneg ####

##  Gives us our true relatedness log likelihoods based on the assumption of no linkage 
##  Good filtering should minimise the effect - can't fully escape it with SNPs

### Desired fpos rate is 0.00002849002 i.e. 2.849 x 10-5

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("PO","FS", "HS", "HAN", "HFC", "U"),
  sim_relats = c("PO","FS", "HS", "HAN", "HFC", "U") 
)

#### Linkage Model ####

## Lets use Pristis pectinata genome as our model for linkage model

n_chromosomes = 46 #chromosome number

L_genome = 2744*0.001 #Gb genome size

sl_diff = 2723055/149968914	#smallest chromosome / largest chromosome


fake_chromo_lengths <- geometric_chromo_lengths(
  n = n_chromosomes,
  L = L_genome,
  sl = sl_diff
)

fake_chromo_lengths$chrom_length_plot

set.seed(42)

afreqs_link <- sprinkle_markers_into_genome(afreqs_ready, fake_chromo_lengths$chrom_lengths)

## Our CKMR model that allows for some linkage in our data - i.e. more overlap between LLRs and more stringent kin-assignmnets

ckmr_link <- create_ckmr(
  D = afreqs_link,
  kappa_matrix = kappas[c("PO","FS", "HS", "HAN", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## Simulate our linkage simulated true kinship values 
## This may take some time as uses MENDEL outside of Rstudio...

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("PO","FS", "HS", "HAN", "HFC", "U"),
  sim_relats = c("PO","FS", "HS", "HAN", "HFC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

#### Check the data for matching genotypes ####

# We allow for up to 500 matching loci here...

matchers <- find_close_matching_genotypes(
  LG = long_genos,
  CK = ckmr,
  max_mismatch = 50
)

matchers # have a geeze

#### Run pairwsie estimates of LLRs for our empirical data #### 

## Make sure that all relationship combinations are consisent with CKMR models

pw_4_LRTs <- lapply(
  X = list(
    POU = c("PO", "U"), 
    POFS = c("PO", "FS"),
    FSU = c("FS", "U"),
    FSHS = c("FS", "HS"),
    HSU = c("HS", "U"),
    HSHAN = c("HS", "HAN"),
    HSHFC = c("HS", "HFC"),
    HFCU = c("HFC", "U")
  ),
  FUN = function(x) {
    pairwise_kin_logl_ratios(
      D1 = long_genos, 
      D2 = long_genos, 
      CK = ckmr,
      numer = x[1],
      denom = x[2],
      num_cores = 8
    )
  }
) %>%
  bind_rows(
    .id = "lr_type"
  ) %>%
  pivot_wider(names_from = lr_type, values_from = logl_ratio)

pw_4_LRTs

##### Extract this for FSP/Unrelated ####

POU_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "PO", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

print(POU_thresholds)

## Visualise the dist of LLRs for non and linked models

PO_U_gg <- Qs %>%
  extract_logls(numer = c(PO = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

PO_U_gg

PO_U_link_gg <- Qs_link_BIG %>% 
  extract_logls(numer = c(PO = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  #ggtitle("FSU/U Linkage Logl Ratio") +
  scale_fill_brewer(palette = "RdYlGn")

PO_U_link_gg + theme_bw()

PO_U_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(PO = 1),
  denom = c(U = 1)
)

PO_U_logls

set.seed(42) # for the jittering

po_fpos <- 4.76e-182
po_fneg <-  0.0001

POU_link_plot <- PO_U_link_gg + 
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = POU, y = -0.002, colour = POU > 290.0195),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3.5
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 650,
           y = 0.013,  
           label = paste("FPR:", po_fpos, "\nFNR:", po_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

print(POU_link_plot)

### PO/FS ###

POFS_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "PO", 
  de = "FS", 
  method = "vanilla", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

print(POFS_thresholds)

## Visualise the dist of LLRs for non and linked models

PO_FS_gg <- Qs %>%
  extract_logls(numer = c(PO = 1), denom = c(FS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

PO_FS_gg + theme_bw()

PO_FS_link_gg <- Qs_link_BIG %>% 
  extract_logls(numer = c(PO = 1), denom = c(FS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  #ggtitle("FSU/U Linkage Logl Ratio") +
  scale_fill_brewer(palette = "RdYlGn") +
  theme_bw()

PO_FS_link_gg + theme_bw()

PO_FS_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(PO = 1),
  denom = c(FS = 1)
)

PO_FS_logls

set.seed(42) # for the jittering

po_fpos <- 4.76e-182
po_fneg <-  0.0001

POFS_link_plot <- PO_FS_link_gg + 
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = POFS, y = -0.002, colour = POFS > 43.78201),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3.5
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 650,
           y = 0.013,  
           label = paste("FPR:", po_fpos, "\nFNR:", po_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

print(POU_link_plot)

## Compute the thresholds we need using MCMC method

FSU_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

print(FSU_thresholds)

# What we need - FPR ~ Lamda* = 161

## Visualise the dist of LLRs for non and linked models

FS_U_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

FS_U_gg

FS_U_link_gg <- Qs_link_BIG %>% 
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  #ggtitle("FSU/U Linkage Logl Ratio") +
  scale_fill_brewer(palette = "RdYlGn")

FS_U_link_gg

FS_U_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(U = 1)
)


set.seed(42) # for the jittering

fs_fpos <- 1.07e-157
fs_fneg <-  0.0001

FSU_link_plot <- FS_U_link_gg + 
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x =FSU, y = -0.002, colour = FSU > 140),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3.5
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 650,
           y = 0.013,  
           label = paste("FPR:", fs_fpos, "\nFNR:", fs_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

print(FSU_link_plot)

emf("C:/Users/samue/Desktop/Honours/FS_UP_LLR_plot.emf", width = 10, height = 8) 
print(FSU_link_plot)
dev.off()

write.csv(FSU_thresholds, "FSU_thresholds_linked.csv")

## Extract Full Sibling Pairs ##

topFS_U <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSU)) %>%
  filter(FSU > 140)

topFS_U

topFS_U$rel <- rep("full-sibs")

set.seed(42)# for the jittering

FS_U_link_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSU, y = -0.002, colour = FSU > 140),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 2
  )

#### Full sib Half sib ####

#Linked
FS_HS_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(FS = 1), denom = c(HS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  scale_fill_brewer(palette = "RdYlGn")

FS_HS_linked_gg + theme_bw()

FS_HS_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(HS = 1)
)

FSHS_MCMC <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "HS", 
  method = "both",
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

print(FSHS_MCMC)

#Quick visual inspection

topFS_HS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > -28.69)

topFS_HS

topFS_HS$rel <- rep("full-sibs")

set.seed(54) # for the jittering

FSHS_plot <- FS_HS_linked_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSHS, y = -0.002, colour = FSHS > -28.69),
    width = 0, 
    height = 0.001,
    shape = 21, 
    fill = NA,
    size = 4
  ) +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  xlab("Log Likelihood Ratio") +
  ylab("Density")+
  labs(fill = "True Relationship") +
  annotate("text", 
           x = 250,
           y = 0.021,  
           label = paste("FPR: 5.86e-143", "\nFNR: 0.001"), 
           hjust = 0, 
           size = 3, 
           color = "black")

print(FSHS_plot)

emf("C:/Users/samue/Desktop/Honours/FS_HS_LLR_plot.emf", width = 10, height = 8) 
print(FSHS_plot)
dev.off()

## So we have our FSPs but we need to keep going with the rest of the individuals 

remaining_pairs <- pw_4_LRTs %>%
  anti_join(bind_rows(topFS_HS), by = c("D2_indiv", "D1_indiv"))

remaining_pairs # work from this now

#### Half sib Unrelated model ####

HS_U_gg <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("HS / UP Logl Ratio")

HS_U_gg

HS_U_logls <- extract_logls(
  Qs,
  numer = c(HS = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
mc_sample_simple(
  Qs,
  nu = "HS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001))

topHS_U <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSU)) %>%
  filter(HSU > 31.3)

topHS_U

set.seed(54) # for the jittering

# Unlinked model

HS_U_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 50.3),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )

### What we really care about is Linked Model

HS_UP_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

HS_UP_linked_gg

HS_UP_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HS = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 

HS_U_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HS_U_thresholds)

topHS_U <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSU)) %>%
  filter(HSU > 16.56556)

topHS_U$rel <- rep("half-sibs")

topHS_U

HSU_fpos <- 1.40e-26
HSU_fneg <- 0.0001 

HSU_plot <- HS_UP_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 16.56556),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 4
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 400,
           y = 0.027,  
           label = paste("FPR:", HSU_fpos, "\nFNR:", HSU_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HSU_plot)

emf("C:/Users/samue/Desktop/Honours/HSU_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSU_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(topHS_U), by = c("D2_indiv", "D1_indiv"))

#### Need to check HS/HAN ####

HS_HAN_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(HAN = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

HS_HAN_linked_gg

HS_HAN_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HS = 1),
  denom = c(HAN = 1)
)

### Get our false positive and false negative rates to guide selection of T 
HS_HAN_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HS", 
  de = "HAN", 
  method = "vanilla", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HS_HAN_thresholds)

## Find the true half sibs 


top_HS_HAN <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSHAN)) %>%
  filter(HSHAN > 5.3567302)

top_HS_HAN
topHS_U

HSFC_fpos <- 4.68e-25
HSFC_fneg <- 0.0001 

HSHAN_plot <- HS_HAN_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSHAN, y = -0.002, colour = HSHAN > 0.544),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 4
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HSHAN_plot)

emf("C:/Users/samue/Desktop/Honours/HSFC_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSFC_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(topHS_U), by = c("D2_indiv", "D1_indiv"))


#### Need to check HS/HFC #### 

HS_HFC_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(HFC = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

HS_HFC_linked_gg

HS_HFC_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HS = 1),
  denom = c(HFC = 1)
)

### Get our false positive and false negative rates to guide selection of T 
HS_HFC_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HS", 
  de = "HFC", 
  method = "both", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HS_HFC_thresholds)

## Find the true half sibs 

top_HS_HFC <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSHFC)) %>%
  filter(HSHFC > 22.070405)

top_HS_HFC
topHS_U

HSFC_fpos <- 0
HSFC_fneg <- 0.0001

HSHFC_plot <- HS_HFC_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSHFC, y = -0.002, colour = HSHFC > 19.3),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 4
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 160,
           y = 0.05,  
           label = paste("FPR:", HSFC_fpos, "\nFNR:", HSFC_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HSHFC_plot)

emf("C:/Users/samue/Desktop/Honours/HSFC_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSFC_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(topHS_U), by = c("D2_indiv", "D1_indiv"))


#### Half First Cousins and Unrelated ####

HFC_U_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HFC= 1),
  denom = c(U = 1)
)


### Get our false positive and false negative rates to guide selection of T 
HFC_U_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HFC", 
  de = "U", 
  method = "vanilla", 
  FNRs = c(0.5, 0.25, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HFC_U_thresholds)

## Select our true HSP

top_HFC_U <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HFCU)) %>%
  filter(HFCU > 4.7)

top_HFC_U$rel <- rep("half-first cousin")
top_HFC_U

HFCU_fpos <- 0.0005
HFCU_fneg <- 0.3

HFC_U_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HFC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)


HFCU_plot <- HFC_U_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HFCU, y = -0.014, colour = HFCU > 3.14),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21, 
    size = 4
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 132,
           y = 0.1,  
           label = paste("FPR:", HFCU_fpos, "\nFNR:", HFCU_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HFCU_plot)

emf("C:/Users/samue/Desktop/Honours/HFCU_LLR_plot.emf", width = 10, height = 8) 
print(HFCU_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(top_HFC_U), by = c("D2_indiv", "D1_indiv"))

#### Merge all relatives ####

library(igraph)
library(ggplot2)
library(dplyr)

all_kin <- bind_rows(topFS_HS, topHS_U)

print(all_kin)

all_kin <- all_kin %>%
  dplyr::select(D2_indiv, D1_indiv, num_loc, rel)

table(all_kin$rel)

colnames(all_kin) <- c("id_1", "id_2", "loc", "rel"); all_kin

## We need to a) check the billabongs b) check birth year 

# first we can deal with selecting the billabong based on the meta file 

meta <- read.csv("C:/Users/samue/Desktop/Honours/analysis/Daly_meta.csv")

all_kin <- all_kin %>%
  left_join(meta, by = c("id_1" = "id")) %>%
rename(birth_year_1 = birthyear, billabong_1 = billabong) %>%
  left_join(meta, by = c("id_2" = "id")) %>%
  rename(birth_year_2 = birthyear, billabong_2 = billabong) %>%
  mutate(
    birth_year_diff = abs(birth_year_1 - birth_year_2),  # Absolute age difference
    same_billabong = billabong_1 == billabong_2) %>%
  select(id_1, id_2, rel, birth_year_diff, same_billabong)

# Last step is to export it 

View(all_kin)

ggplot() +
  geom_bar(data = all_kin, aes(x = ))

write.csv(all_kin, "C:/Users/samue/Desktop/CKMR_kin.csv")

#### Visualise as a network plot with extra relatives ####

## write our sibling results as csv 

sibs <- read.csv("C:/Users/samue/Desktop/CKMRsim_kin_metadata.csv")
meta <- read.csv("C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv")

kinNWdata <- sibs %>%
  dplyr::select(id_1, id_2, rel, birth_yr_diff, pw_loc, pw_loc_word)

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE); print(network) # works

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, Sex, birth_year, Site_ID)

vertices <- as_tibble(vertices, what = "vertices")

vertices

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)

## Plot the network ##

kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = kinNWdata$birth_yr_diff,
        edge_colour = factor(kinNWdata$rel), 
        edge_linetype = kinNWdata$pw_loc_word))+
  ggraph::scale_edge_width(range = c(3, 1), breaks = seq(0,12, 1), name = "Cohort Gap") +
  ggraph::scale_edge_linetype_manual(values = c("dashed", "solid"), 
                                     name = "Capture Location", 
                                     aesthetics = "edge_linetype") +
  ggraph::scale_edge_colour_manual(values = c("skyblue", "orange"),
                                   name = "Kin Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(shape = Sex),
                          size = 5) +
  ggplot2::scale_color_manual(values = adegenet::funky(9), na.value = "grey50") +
  labs(shape = "Sex") +
  ggplot2::theme_bw() +
  #guides(edge_width = "none") +
  ggplot2::theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    plot.margin = unit(rep(1,5), "cm"), 
    margin = margin(20, 0, 40, 10))

print(kin_network1)

## Save it as an EMF

devEMF::emf("C:/Users/samue/Desktop/CKMRsim_Network.emf", width = 10, height = 8)  
print(kin_network1)
dev.off()

ggraph::geom_node_text( aes(label = df$id), repel = TRUE, 
                        size = 5, color = "black") +

#### Permutational Test for differences between sample year and billabong #### 


## We need to make a seperate file for all pairwise comparisons with the same format as CKMR_kin

ind_metadata <- read.csv("C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv")

ind_metadata <- ind_metadata %>%
  mutate(year = year(dmy(Date)))

remaining_pairs$rel <- rep("unrelated")

all_pairs <- bind_rows(topFS_U, topHS_U, top_HFC_U, remaining_pairs)

print(all_pairs)

all_pairs <- all_pairs %>%
  dplyr::select(D2_indiv, D1_indiv, rel)

colnames(all_pairs) <- c("id_1", "id_2", "rel"); all_pairs

write.csv(all_pairs, "C:/Users/samue/Desktop/CKMRsim_all.csv")

all_pairs2 <- all_pairs %>%
  left_join(meta, by = c("id_1" = "id")) %>%
  rename(birth_year_1 = birthyear, billabong_1 = billabong, capture_year1 = Year_caught) %>%
  left_join(meta, by = c("id_2" = "id")) %>%
  rename(birth_year_2 = birthyear, billabong_2 = billabong, capture_year2 = Year_caught) %>%
  mutate(
    birth_year_diff = abs(birth_year_1 - birth_year_2),  # Absolute age difference
    Within_Billabong = as.integer(billabong_1 == billabong_2), # returns true false as 1,0
    Within_Cohort = as.integer(capture_year1 == capture_year2)) %>% #returns true false as 1,0
  select(id_1, id_2, rel,billabong_1, billabong_2, birth_year_1, birth_year_2, capture_year1, capture_year2, birth_year_diff, Within_Billabong, Within_Cohort) # returns desired metadata

all_pairs2 <- all_pairs %>%
  dplyr::left_join(ind_metadata %>%
                     dplyr::select(id, Sex, Site_ID, year), 
                   by = c("id_1" = "id")) %>%
  dplyr::mutate(year_1 = year, loc_1 = Site_ID, .keep = "unused") %>% 
  dplyr::left_join(ind_metadata %>%
                     dplyr::select(id, Sex, Site_ID, year), 
                   by = c("id_2" = "id")) %>%
  dplyr::mutate(year_2 = year, loc_2 = Site_ID, .keep = "unused")


all_pairs2

# need to add kin/underlated as 1,0 

all_pairs2$Relatives <- ifelse(all_pairs2$rel == c("unrelated"),
                                yes = 0, 
                                no = 1)

all_pairs2$pw_loc <- ifelse(all_pairs2$loc_1 == all_pairs2$loc_2, "1", "0")

all_pairs2$pw_yr <- ifelse(all_pairs2$year_1 == all_pairs2$year_2, "1", "0")

all_pairs2 <- all_pairs2 %>%
  mutate(across(c(Within_Billabong, Within_Cohort, Relatives), as.integer))

View(all_pairs2)
# looks good 

## lets save this for statistical testing 

write.csv(all_pairs2, "C:/Users/samue/Desktop/CKMRsim_all_pairs2.csv")

# Load in parallel package for parallel computing

all_data <- read.csv("C:/Users/samue/Desktop/Honours/analysis/CKMRsim_all_pairs.csv")

library(parallel)

set.seed(42)  # For reproducibility


