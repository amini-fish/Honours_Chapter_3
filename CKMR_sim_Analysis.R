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

setwd("C:/Users/samue/Desktop/Honours/analysis")

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

gl.write.csv(gl, outfile = "outfile.csv", outpath = "C:/Users/samue/Desktop/Honours/analysis", verbose = NULL)

data <- read.csv("outfile.csv")

glimpse(data)

#### Prepare the data ####

## SKIP TO LINE 107 TO START ANALYSIS ## 

gl2vcf(gl, 
      plink.bin.path = "C:/Users/samue/Desktop/Honours/analysis/plink", 
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

genotypes_2

genotypes_3 <- data.frame(genotypes_2)

genotypes_3 <- genotypes_3 %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = "_")

length(gl@ind.names)

genotypes_3 <- genotypes_3[-1,]

genotypes_3

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

#write.csv(long_genos, file = "Pristis_genofor_CKMRsim.csv" )

#### START HERE - Long genos ####

long_genos <- read.csv("Pristis_genofor_CKMRsim.csv"); long_genos

long_genos ## Looking MINT!!! We got there! 

long_genos <- long_genos[,-1]
locus_names <- unique(long_genos$Locus)

afreqs_ready <- long_genos %>%
  count(Locus, Allele) %>%  
  group_by(Locus) %>%
  mutate(
    Freq = n / sum(n),
    Chrom = "Unk",
    Pos = as.integer(factor(Locus, levels = locus_names))
  ) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
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
  kappa_matrix = kappas[c("FS", "HS","HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

ckmr

#### Subset Kappas - IBD probs for desired relats ####

kappas # all possible

kappas[c("FS", "HS", "HFC", "U"), ]

#### Simulate dyads for each to approximate the LLRs and therefore the fpos/fneg ####

##  Gives us our true relatedness log likelihoods based on the assumption of no linkage 
##  Good filtering should minimise the effect - can't fully escape it with SNPs

### Desired fpos rate is 0.00002849002 i.e. 2.849 x 10-5

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS", "HFC", "U") 
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
  kappa_matrix = kappas[c("FS", "HS", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## Simulate our linkage simulated true kinship values 
## This may take some time as uses MENDEL outside of Rstudio...

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS",  "HFC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

#### Check the data for matching genotypes ####

# We allow for up to 500 matching loci here...

matchers <- find_close_matching_genotypes(
  LG = long_genos,
  CK = ckmr,
  max_mismatch = 500
)

matchers # have a geeze

#### Run pairwsie estimates of LLRs for our empirical data #### 

## Make sure that all relationship combinations are consisent with CKMR models

pw_4_LRTs <- lapply(
  X = list(
    FSU = c("FS", "U"),
    FSHS = c("FS", "HS"),
    HSU = c("HS", "U"),
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

fs_fpos <- 3.49e-137
fs_fneg <-  0.0001

FSU_link_plot <- FS_U_link_gg + 
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x =FSU, y = -0.002, colour = FSU > 161),
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
  filter(FSU > 161)

topFS_U

topFS_U$rel <- rep("full-sibs")

set.seed(42)# for the jittering

FS_U_link_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSU, y = -0.002, colour = FSU > 190),
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

FS_HS_linked_gg

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
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001)
)

print(FSHS_MCMC)

#Quick visual inspection

topFS_HS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > -7.99)

topFS_HS

topFS_HS$rel <- rep("full-sibs")

set.seed(54) # for the jittering

FSHS_plot <- FS_HS_linked_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSHS, y = -0.002, colour = FSHS > -7.99),
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
  filter(HSU > 54)

topHS_U

set.seed(54) # for the jittering

# Unlinked model

HS_UP_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 54),
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
  filter(HSU > 12.3)

topHS_U$rel <- rep("half-sibs")

topHS_U

HSU_fpos <- 1.40e-26
HSU_fneg <- 0.0001 

HSU_plot <- HS_UP_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 12.3),
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
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HS_HFC_thresholds)

## Find the true half sibs 

top_HS_HFC <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSHFC)) %>%
  filter(HSHFC > -11.5)

top_HS_HFC
topHS_U

HSFC_fpos <- 4.68e-25
HSFC_fneg <- 0.0001 

HSFC_plot <- HS_FC_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSFC, y = -0.002, colour = HSFC > -20.8),
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

plot(HSFC_plot)

emf("C:/Users/samue/Desktop/Honours/HSFC_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSFC_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(topHS_U), by = c("D2_indiv", "D1_indiv"))


#### Half First Cousins and Unrelated ####


HFC_U_linked_gg

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
  method = "IS", 
  FNRs = c(0.5, 0.25, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)) 

print(HFC_U_thresholds)

## Select our true HSP

top_HFC_U <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HFCU)) %>%
  filter(HFCU > 3.14)

top_HFC_U$rel <- rep("half-first cousin")
top_HFC_U

HFCU_fpos <- 0.0033
HFCU_fneg <- 0.25

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

all_kin <- bind_rows(topFS_U, topHS_U, top_HFC_U)

print(all_kin)

all_kin <- all_kin %>%
  select(D2_indiv, D1_indiv, num_loc, rel)

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

write.csv(all_kin, "CKMR_kin.csv")

#### Visualise as a network plot with extra relatives ####

## write our sibling results as csv 

sibs <- read.csv("C:/Users/samue/Desktop/Honours/analysis/CKMR_kin.csv")
meta <- read.csv("C:/Users/samue/Desktop/Honours/analysis/Daly_meta.csv")

meta <- meta[meta$id %in% c(sibs$id_1, sibs$id_2),]
meta

kinNWdata <- sibs %>%
  select(id_1, id_2, rel, birth_year_diff, same_billabong)

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE); print(network) # works

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, birthyear, billabong)

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
    aes(width = kinNWdata$birth_year_diff,
        edge_colour = factor(kinNWdata$rel), 
        edge_linetype = kinNWdata$same_billabong, 
        edge_alpha = 1))+
  ggraph::scale_edge_width(range = c(2, 3), breaks = c(0,1, 6), name = "Cohort Gap") +
  ggraph::scale_edge_linetype_manual(values = c("dashed", "solid"), 
                                     name = "Capture Location", 
                                     aesthetics = "edge_linetype") +
  ggraph::scale_edge_colour_manual(values = c("skyblue", "red3", "orange"),
                                   name = "Kin Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(shape = sex),
                          size = 5) +
  ggraph::geom_node_text( aes(label = df$id), repel = TRUE, 
                          size = 5, color = "black") +
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

devEMF::emf("C:/Users/samue/Desktop/Honours/CKMRsim_Network.emf", width = 10, height = 8)  
print(kin_network1)
dev.off()

#### Permutational Test for differences between sample year and billabong #### 


## We need to make a seperate file for all pairwise comparisons with the same format as CKMR_kin

meta <- read.csv("C:/Users/samue/Desktop/Honours/analysis/Daly_meta.csv")

remaining_pairs$rel <- rep("unrelated")

all_pairs <- bind_rows(topFS_U, topHS_U, top_HFC_U, remaining_pairs)

print(all_pairs)

all_pairs <- all_pairs %>%
  select(D2_indiv, D1_indiv, rel)

colnames(all_pairs) <- c("id_1", "id_2", "rel"); all_pairs

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

# need to add kin/underlated as 1,0 

all_pairs2$Relatives <- ifelse(all_pairs2$rel == c("unrelated"),
                                yes = 0, 
                                no = 1)


all_pairs2 <- all_pairs2 %>%
  mutate(across(c(Within_Billabong, Within_Cohort, Relatives), as.integer))

View(all_pairs2)
# looks good 

## lets save this for statistical testing 

write.csv(all_pairs2, "CKMRsim_all_pairs.csv")

# Load in parallel package for parallel computing

all_data <- read.csv("C:/Users/samue/Desktop/Honours/analysis/CKMRsim_all_pairs.csv")

library(parallel)

set.seed(42)  # For reproducibility

# Function to calculate the difference in proportions
calc_diff <- function(data) {
  within_prop <- mean(data$Relatives[data$Within_Billabong == 1 & data$Within_Cohort == 1])
  between_prop <- mean(data$Relatives[data$Within_Billabong == 0 | data$Within_Cohort == 0])
  return(within_prop - between_prop)
}

# Compute observed difference
observed_diff <- calc_diff(all_data)

# Number of permutations
n_perm <- 10000
n_cores <- detectCores() - 1  # Use all available cores except 1

all_data

cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c("all_data", "calc_diff"))

# Parallelized permutation function
permute_function <- function(i, data) {
  permuted_data <- data  # Work on a copy
  permuted_data$billabong_1 <- sample(permuted_data$billabong_1)
  permuted_data$billabong_2 <- sample(permuted_data$billabong_2)
  permuted_data$capture_year1 <- sample(permuted_data$capture_year1)
  permuted_data$capture_year2 <- sample(permuted_data$capture_year2)
  
  # Recalculate within-billabong and within-cohort
  permuted_data$Within_Billabong <- as.integer(permuted_data$billabong_1 == permuted_data$billabong_2)
  permuted_data$Within_Cohort <- as.integer(permuted_data$capture_year1 == permuted_data$capture_year2)
  
  calc_diff(permuted_data)
}

# Run parallel permutations with data explicitly passed
perm_diffs <- parSapply(cl, 1:n_perm, permute_function, data = all_data)
# Run permutations in parallel

stopCluster(cl)  # Stop the cluster when done

# Compute p-value
p_value <- mean(perm_diffs >= observed_diff)

# Results
cat("Observed Difference:", observed_diff, "\n")
cat("P-value:", p_value, "\n")

# Convert permuted differences to a dataframe for plotting
perm_diffs_df <- data.frame(perm_diffs)

# Create the plot
ggplot(perm_diffs_df, aes(x = perm_diffs)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black", alpha = 0.7) +  # Histogram of permutations
  geom_vline(aes(xintercept = observed_diff), color = "red", linetype = "dashed", size = 1) +  # Observed difference
  labs(
    title = "Permutation Test: Observed vs. Null Distribution",
    x = "Difference in Kinship Proportions (Within - Between)",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(perm_diffs_df, aes(x = perm_diffs)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(aes(xintercept = observed_diff), color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Permutation Test Density Plot",
    x = "Difference in Kinship Proportions",
    y = "Density"
  ) +
  theme_minimal()


#### GLMM ##### 

library(lme4)

glimpse(all_data)

# Fit GLMM model
glmm_model <- glmer(Relatives ~ Within_Billabong + Within_Cohort + 
                      (1 | id_1) + (1 | id_2), 
                    data = all_data,
                    family = binomial, 
                    control = glmerControl(optimizer = "bobyqa"))

summary(glmm_model)

all_data$predicted <- predict(glmm_model, type = "response")

library(ggplot2)

# Visualizing predicted probabilities based on Within_Billabong and Within_Cohort
glmm_plot <- ggplot(all_data, aes(x = factor(Within_Billabong), y = predicted, color = factor(Within_Cohort))) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 4) + 
  stat_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +
  labs(x = "Same Location", y = "Predicted Probability of Kinship", color = "Same Cohort") +
  theme_bw() +
  scale_colour_manual(values = c("darkolivegreen2", "grey")) + 
  theme(legend.position = "top")

print(glmm_plot)
#### Troubleshooting ####
ggplot(all_data, aes(x = Within_Billabong, y = Within_Cohort, color = Relatives)) + 
  geom_jitter(size = 3) + 
  labs(title = "Pairwise Relationships")

# Correlation matrix
cor(all_data[, c("Within_Billabong", "Within_Cohort", "birth_year_diff")])


## Try with interaction term 

glmer_model_interaction <- glmer(Relatives ~ Within_Billabong * Within_Cohort + 
                                   (1 | id_1) + (1 | id_2), 
                                 data = all_data, family = binomial)
summary(glmer_model_interaction)


AICctab(glmm_model, glmer_model_interaction)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/kinship_glmm.emf", width = 10, height = 8)  # Set the width and height in inches
print(glmm_plot)
dev.off()


## Using ggpredict 

preds <- ggpredict(glmm_model, terms = c("Within_Billabong", "Within_Cohort"), bias_correction = T, type = "fixed", back.transform = FALSE)
glimpse(preds)

preds$odds <- exp(preds$predicted)
preds$odds_low <- exp(preds$conf.low)
preds$odds_high <- exp(preds$conf.high)

# Plotting odds
glmm_plot2 <- ggplot(preds, aes(x = x, y = odds, color = group)) +
  geom_point(size = 4.5) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = odds_low, ymax = odds_high, fill = group), alpha = 0.2, color = NA) +
  scale_y_log10() +  # Plot odds on log10 scale
  labs(
    x = "Within Billabong (0 = no, 1 = yes)",
    y = "Odds of Being Relatives (log10 scale)",
    color = "Within Cohort",
    fill = "Within Cohort"
  ) +
  scale_colour_manual(values = c("darkolivegreen2", "grey")) + 
  scale_fill_manual(values = c("darkolivegreen3", "grey")) + 
  scale_x_continuous(breaks = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14), 
        legend.position = "top")


print(glmm_plot2)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/kinship_glmm_2.emf", width = 10, height = 8)  # Set the width and height in inches
print(glmm_plot2)
dev.off()

## Updated plot 
library(ggplot2)

# Plot with jittered empirical points and predicted odds
ggplot() +
  # Jitter the empirical data points
  geom_jitter(data = all_data, 
              aes(x = Within_Billabong, y = Relatives, color = as.factor(Within_Cohort)), 
              width = 0.1, height = 0.1, alpha = 0.6, size = 3) +
  
  # Add predicted odds from the model
  geom_line(data = preds, aes(x = x, y = odds, color = group), size = 1) +
  geom_ribbon(data = preds, aes(x = x, ymin = odds_low, ymax = odds_high, fill = group), alpha = 0.2) +
  
  # Customize the plot
  labs(
    x = "Within Billabong (0 = no, 1 = yes)",
    y = "Predicted Odds of Being Relatives",
    color = "Within Cohort",
    fill = "Within Cohort"
  ) +
  theme_minimal() +
  scale_y_log10() +  # Log10 scale for the odds
  scale_x_continuous(breaks = c(0, 1)) +  # Set x-axis breaks to 0 and 1
  theme(legend.position = "top")
