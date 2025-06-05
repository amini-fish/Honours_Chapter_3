library(dartRverse)
library(related)
library(tidyverse)
library(dplyr)

## Set the WD

setwd("C:/EMIBD/for_Jinliang")

sim_data <- readLines("simulated_long.ibd9")

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
tmp_data_raw_3$V3 <- lapply(tmp_data_raw_3$V3, as.numeric)

colnames(tmp_data_raw_3) <- tmp_headings[2:22]

df <- data.frame(ind1=tmp_data_raw_3$Indiv1, ind2=tmp_data_raw_3$Indiv2,rel= tmp_data_raw_3$`r(1,2)`)

df$rel <- as.numeric(df$rel); df

df$index <- ifelse(df$ind1 == df$ind2, "drop", "keep")

df <- df %>% 
  filter(index == "keep")%>%
  droplevels

df <- df[,-4]

write.csv(df, "C:/EMIBD/for_Jinliang/EM1_emibd9_simdata.csv") # so we don't have to worry about reading it in the long way when it crashes...

## Okay now that we have our EMIBD9 simulations loaded into R (what a fucking task that was)

# We need to remove all the pairwise comparisons that don't match up with those from coancestry for consistency (and computation time)

## OUr template data

####

##### START HERE #####

####

df <- read.csv("C:/EMIBD/for_Jinliang/EM1_emibd9_simdata.csv")

view(df)

dim(df)

co_sims <- read.csv("C:/Users/samue/Desktop/Relatedness_simulations/coancestry_sims/coancestry_sims_results_2.csv", stringsAsFactors = T)

pairwise_index <- co_sims %>%
  unite(merger, Ind1, Ind2)

pairwise_index

df_x <- df %>%
  unite(merger2, ind1, ind2)

df_x

filtered_df <- filter(df_x, merger2 %in% pairwise_index$merger)

dim(filtered_df)

# Whilst the previous code I wrote works for the EMIBD9 output, I'd rather just do it a new way LOL

#select all the FS 

#NOW WE NEED TO ASSIGN THE KNOWN SIBTYPES TO EACH DYAD FROM com_sims

pairwise_index$Sibtype

## need to get in the same order 

# try dplyr::arrange()

EM1_emibd9_df <- arrange(filtered_df, merger)

co_sims_df <- arrange(pairwise_index, merger)

head(EM1_emibd9_df)

dim(EM1_emibd9_df)
head(co_sims_df)

EM1_emibd9_df$Sibtype <- co_sims_df$Sibtype

##----------------------- Now rinse and repeat for EM2 -----------------------##

setwd("C:/Users/samue/Desktop/Relatedness_simulations/coancestry_sims/sims_6")

sim_data <- readLines("sims_6.ibd9")

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
tmp_data_raw_3$V3 <- lapply(tmp_data_raw_3$V3, as.numeric)

colnames(tmp_data_raw_3) <- tmp_headings[2:22]

df <- data.frame(ind1=tmp_data_raw_3$Indiv1, ind2=tmp_data_raw_3$Indiv2,rel= tmp_data_raw_3$`r(1,2)`)

df$rel <- as.numeric(df$rel); df

df$index <- ifelse(df$ind1 == df$ind2, "drop", "keep")

df <- df %>% 
  filter(index == "keep")%>%
  droplevels

df <- df[,-4]

write.csv(df, "EM2_emibd9_simdata.csv")

df_2 <- read.csv("EM2_emibd9_simdata.csv")

co_sims <- read.csv("C:/Users/samue/Desktop/Relatedness_simulations/coancestry_sims/coancestry_sims_results_2.csv", stringsAsFactors = T)

pairwise_index <- co_sims %>%
  unite(merger, Ind1, Ind2)

pairwise_index

df_2
df_2 <- df_2 %>%
  unite(merger, ind1, ind2)

filtered_df_2 <- filter(df_2, merger %in% pairwise_index$merger)

dim(filtered_df_2)

pairwise_index$Sibtype

## need to get in the same order 

# try dplyr::arrange()

Em2_emibd9_df <- arrange(filtered_df_2, merger)

co_sims_df <- arrange(pairwise_index, merger)

head(Em2_emibd9_df)
head(co_sims_df)

Em2_emibd9_df$Sibtype <- co_sims_df$Sibtype


## --------------------------------------------------------------------------


## Hell yeah it's sorted

write.csv(emibd9_df, "C:/Users/samue/Desktop/Relatedness_simulations/simulated_emibd9_EM1_rel.csv")

write.csv(Em2_emibd9_df, "C:/Users/samue/Desktop/Relatedness_simulations/simulated_emibd9_EM2_rel.csv")

EM1_EMIBD9 <- emibd9_df$rel; EMIBD9

EM2_EMIBD9 <- Em2_emibd9_df$rel

simulations_combined <- cbind(co_sims_df, EM1_EMIBD9, EM2_EMIBD9)

simulations_combined

simrel <- read.csv("C:/Users/samue/Desktop/Relatedness_simulations/simulated_rel.csv")

simrel <- cbind(simrel, EM2_EMIBD9)

write.csv(simrel, "C:/Users/samue/Desktop/Relatedness_simulations/simulated_rel.csv") ## This is what we load into "relatedness_simulations.Rdata"

## Plot the simulated values where each dyad has #1000 inds

simplot <- ggplot(emibd9_df, aes(x = rel, fill = Sibtype)) +
  geom_density(alpha = 0.5, position = "identity")+
  theme_bw()+
  ggtitle("Distribution of simulated relatedness coefficients") +
  xlab("Relatedness") +
  ylab("Frequency")+
  scale_fill_brewer(palette = "Dark2")

simplot <- simplot + theme(plot.title = element_text(hjust = 0.5))

print(simplot)
