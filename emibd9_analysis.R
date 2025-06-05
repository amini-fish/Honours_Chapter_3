#### Welcome brave traveler...this script will be used to run the EMIBD9 program through dartR.captive ####

#### Install packages ####
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

devtools::install_github("green-striped-gecko/dartR.captive@dev")

#### Set working directory ####

setwd("C:/Users/samue/Desktop/Honours/analysis")

#### Load Packages ####

library(SNPRelate)
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

#### Load the filtered genotype data in - for more info on how the data has been filtered use basics ####

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

#gl.report.basics(gl)

#### Run our analysis - using EMIBD9 implementation #### 

## Calculate the number of independant pairwise comps where K(K-1)/2

nInd(gl)

29*(29-1)/2

## Run it via existing function - or can run from CMD 

daly.rel <- dartR.captive::gl.run.EMIBD9(gl, 
                          Inbreed = 1, 
                          emibd9.path =  "C:/EMIBD")

#### Results ####

diag(daly.rel$rel) <- 0 # removes self comparisons hooray

emibd.rel <- pw_mat_to_df(daly.rel$rel) # convert the square matrix to an edge based dataframe

emibd.rel <- emibd.rel[,c(1,2,4)] # remove your third column (not needed and is confusing)

## Check your results: 

View(emibd.rel)

save(emibd.rel, file = "emibd_results.csv") # Not essential but helps as a just in case

#### Summary statistics ####

mean.rel <- mean(as.numeric(emibd.rel$value))  
sd.rel <- sd(as.numeric(emibd.rel$value)) 
med <- median(as.numeric(emibd.rel$value))

mean.rel #0.013
sd.rel #0.033
med #0.0052

#### Visualise the distribition of results as a histogram ####

rel.hist <- ggplot(data = emibd.rel, aes(x = as.numeric(value))) +
  geom_histogram(bins = 80, col = I("black")) +
  xlab("Relatedness") +
  ylab("Count") +
  theme_bw()

rel.hist <- rel.hist + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 11)) +
  geom_vline(xintercept = med, linewidth = 1.1, col = "green") +
  geom_vline(xintercept = 0.125, linewidth = 1.1, col = "orange") +
  geom_vline(xintercept = 0.250, linewidth = 1.1, col = "skyblue") +
  geom_vline(xintercept = 0.092, col = "black", linewidth = 1, linetype="dotted") +
  geom_vline(xintercept = 0.158, col = "black", linewidth = 1,  linetype="dotted") +
  geom_vline(xintercept = 0.204, col = "black", linewidth = 1, linetype="dotted") +
  geom_vline(xintercept = 0.296, col = "black", linewidth = 1,  linetype="dotted") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(n.breaks = 10) 
  
print(rel.hist)

emf("C:/Users/samue/Desktop/Honours/EMIBD9_Histogram.emf", width = 10, height = 8)
print(rel.hist)
dev.off()

#### NETWORK PLOT FOR KIN ONLY ####

## write our sibling results as csv 

#write.csv(emibd.siblings, "emibd_sibs_daly.csv")

sibs  <- read.csv("C:/Users/samue/Desktop/Honours/analysis/emibd_sibs_daly.csv")
meta <- read.csv("C:/Users/samue/Desktop/Honours/analysis/Daly_meta.csv")

meta <- meta[meta$id %in% c(sibs$id_1, sibs$id_2),]

data <- sibs

kinNWdata <- data %>% 
  #dplyr::mutate(billabong_both, paste(df$n, "-", df$s))%>%
  #dplyr::mutate(Cohort_gap = cohort1 - cohort2) %>%
  #dplyr::select(id_1, id_2, relatedness, Billabong_ID2) #, Cohort_gap

meta

#Calculate the gap between cohorts, and add it into a new dataframe

#kinNWdata$Cohort_gap[kinNWdata$Cohort_gap<0] <- 1 # this replaces all of our "negative" cohort gaps so that we have a tidy plot

billabong <- c("different", "same", "same", "different", "same", "different", "different", "same", "same", "same")

kinNWdata$billabong <- billabong

kinNWdata

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE); network

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, Cohort)

vertices <- as_tibble(vertices, what = "vertices")

vertices

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)


kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = 2,
        edge_colour = factor(relatedness), 
        edge_linetype = billabong),
    #arrow = arrow(length = unit(3, 'mm')), 
    #end_cap = ggraph::circle(2, 'mm'),
    edge_alpha = 1) +
  #ggraph::scale_edge_width(range = c(2), breaks = c(0,1), name = "Cohort Gap") +
  ggraph::scale_edge_linetype_manual(values = c("dotted", "solid"), 
                                     name = "Capture Location", 
                                     aesthetics = "edge_linetype") +
  ggraph::scale_edge_colour_manual(values = c("skyblue", "orange"),
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
  guides(edge_width = "none") +
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

emf("C:/Users/samue/Desktop/Honours/EMIBD9_Network.emf", width = 10, height = 8)  
print(kin_network1)
dev.off()

#### Misc code ####

# Combine together
ibd9DT <- as.matrix(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

# Combine together
ibd9DT <- data.frame(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

mean(as.numeric(ibd9DT$r.1.2.)) ## Mean is 0.0129 
sd(as.numeric(ibd9DT$r.1.2.)) ## SD is +- 0.0323

# Test for within and between cohort means, note this isn't perfect
group_mean<- aggregate(x= as.numeric(ibd9DT$r.1.2.),
                       # Specify group indicator
                       by = list(ibd9DT$CC),      
                       # Specify function (i.e. mean)
                       FUN = mean)

group_mean

group_sd <- aggregate(x= as.numeric(ibd9DT$r.1.2.),
                       # Specify group indicator
                       by = list(ibd9DT$CC),      
                       # Specify function (i.e. mean)
                       FUN = sd)

#A quick look at the spread of rel values...

hist(as.numeric(ibd9DT$r.1.2.), breaks = 250)

## Putative siblings - we use the 95% CI from the simulations to groundtruth our sibling relationships as mean rel is 0.01

emibd.sibs <- ifelse(emibd.rel$value >= 0.092 & emibd.rel$value <= 0.180, 
               yes = "hsp", 
               no = ifelse(emibd.rel$value >=0.204 & emibd.rel$value <= 0.296, 
                           yes = "fsp",
                           "unelated"))

#Add the sibling assignments to a new data frame 

emibd.results <- data.frame(cbind(emibd.rel, emibd.sibs)); emibd.results

####  Assign kin to sibling network  #####

#hsps - extract

half.sibs <- subset(emibd.results, emibd.sibs == "hsp"); half.sibs

#fsps - extract 

full.sibs <- subset(emibd.results, emibd.sibs == "fsp"); full.sibs

# ALLLLL stitched together now 

emibd.siblings <- rbind(half.sibs, full.sibs); sibs.all

# Remove duplicated pairs (i.e., Ab & BA)

sibs.all[!duplicated(sibs.all$r.1.2.), ] 
sibs.all

## Fucking around

emibd.rel

write.csv(emibd.rel, "daly_rel_all.csv")

relate <- read.csv("daly_rel_all.csv"); relate

meta2 <- read.csv("Daly_meta.csv")

## Distance 

library(sf)

df = read.csv('Daly_meta.csv') %>%
st_as_sf(coords=c("Long","Lat"), crs=4326)

df <- df%>%st_transform(3857)

df

distance <- st_distance(x = df, by_element = F, which = "Euclidean")

View(distance)

dist_km <- distance/1000 #as km

dist_km

meta2$id

rownames(dist_km) <- meta2$id
colnames(dist_km) <- meta2$id

dist_km <- as.data.frame(dist_km)

dist_km <- pw_mat_to_df(as.matrix(dist_km))

View(dist_km)

## now we need to merge the distance matrix as a column to the other data - so we can run corMLPE

dist_data <- dist_km %>% arrange(id_link)

dist_data <- rownames_to_column(dist_km, "X")

dist_data

dim(dist_data)
dist <- dist_data[, 5]

## Both datasets in same oder...

meta3 <- read.csv("pairs_meta.csv")

meta3 <- meta3 %>% arrange(X)

head(meta3)

## get pairwise year of capture comparisons (thanks Mirjam)

meta3$year_caught_diff <- ifelse(meta3$Year_caught_id1 == meta3$Year_caught_id2, "Same", "Different")

meta3 <- meta3 %>%
  unite("year_caught_both",
      as.character("Year_caught_id1"),as.character("Year_caught_id2"), 
      sep = ":",
      remove = FALSE, na.rm = FALSE) %>%
  mutate(year_caught_both=factor(year_caught_both))

## Do the same for Catch Set

meta3$catch_set_diff <- ifelse(meta3$Catch.Set_ID1 == meta3$Catch.Set_ID2, "Same", "Different")

meta3 <- meta3  %>%
unite("catch_set_both",
      as.character("Catch.Set_ID1"),as.character("Catch.Set_ID2"), 
      sep = ":",
      remove = FALSE, na.rm = FALSE) %>%
  mutate(catch_set_both=factor(catch_set_both))

## Do the same for Billabong

meta3$billabong_diff <- ifelse(meta3$Billabong_ID1 == meta3$Billabong_ID2, "Same", "Different")

meta3 <- meta3 %>% 
  unite("billabong_both",
        as.character("Billabong_ID1"),as.character("Billabong_ID2"), 
        sep = ":",
        remove = FALSE, na.rm = FALSE) %>%
  mutate(billabong_both=factor(billabong_both))

meta4 <- cbind(meta3, dist)

write.csv(meta4, "corMLPE_data.csv")


