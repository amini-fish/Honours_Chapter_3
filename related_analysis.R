setwd("C:/Users/samue/Desktop/Honours/analysis")

### INSTALL WHERE NEEDED (assumed starting from new PC for reporducing) ### 

install.packages("hierfstat")
devtools::install_version("ggplot2", "3.4.4")
#install.packages("related")
install.packages("related", repos="http://R-Forge.R-project.org")
install.packages("parrallel")

### LOAD DESIRED PACKAGES INTO R ###

library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(viridis)
library(related)

## Not run: 
#---Read data into R---#

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl


input <- dartR.base::gl2related(gl, save = F)

#---Calculate Relatedness---#

output <- coancestry(genotype.data = input, 
                     trioml = 1,
                     allow.inbreeding = T,
                     rng.seed = 42)

#---View Point Estimates---#
output$relatedness

## End(Not run)	

### Save trioML results into a new dataframe

trioML <- data.frame(output$relatedness[, c(2,3,5)])

mean <- mean(trioML$trioml)
sd <- sd(trioML$trioml)

mean
sd


## Simulate 95% CI for relatives from file - relatedness_simulations.R  

print(ci_FS)
print(ci_HS)

## Extract the indiviudals with relatedness above x and y -  need to find a solid guideline for this assignment not arbitrary numbers.

pull.kin <- ifelse(trioML$trioml >= 0.211495 & trioML$trioml <= 0.279125, 
               yes = "hsp", 
               no = ifelse(trioML$trioml >=0.4574700  & trioML$trioml <= 0.5203725, 
                           yes = "fsp",
                           "unelated"))

## Combine the two dataframes so we have all putative hsps and fsps identified.

trioML <- data.frame(cbind(trioML, pull.kin))

## from this data we can extract our siblings and put them in a new data frame

trioML.fsp <- trioML[ which(trioML$pull.kin == "fsp"),]
trioML.hsp <- trioML[ which(trioML$pull.kin == "hsp"),]

trioML.sibs<- rbind(trioML.fsp, trioML.hsp)

trioML.sibs

###############
## Visualise ##
###############

rel.hist <- ggplot(data = trioML, aes(x = trioml)) +
  geom_histogram(bins = 100, col = I("black")) +
  labs(title = "Histogram of kinship coefficients (θ) between 2012 & 2013 sawfish") +
  xlab("Relatedness (R)") +
  ylab("Count")
  
print(rel.hist)


############# EXTRA NONESENSE ###################
theme(plot.title = element_text(hjust = 0.5, size = 22), 
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      axis.text = element_text(size = 13)) +
  geom_vline(xintercept = mean, size = 1.1, col = "green")
#geom_vline(xintercept = 0.125, size = 1.1, col = "orange") +
#geom_vline(xintercept = 0.250, size = 1.1, col = "skyblue") +
geom_vline(xintercept = (-sd), col = "black", size = 1, linetype="dotted") +
  geom_vline(xintercept = (sd), col = "black", size = 1,  linetype="dotted") +
  #geom_vline(xintercept = 0.204, col = "black", size = 1, linetype="dotted") +
  #geom_vline(xintercept = 0.296, col = "black", size = 1,  linetype="dotted") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(n.breaks = 10)


## Make a network plot for this analysis type too 

#write.csv(trioML.sibs, "tioML_sibs_daly.csv")

#we now need to manually add birth cohorts in (important to get right)

#reload baack into R

sibs  <- read.csv("tioML_sibs_daly.csv")
meta <- read.csv("Daly_meta.csv")
sibs

meta <- meta[meta$id %in% c(sibs$id_1, sibs$id_2),]
data <- sibs
kinNWdata <- data %>% 
  dplyr::mutate(Cohort_gap = cohort1 - cohort2) %>%
  dplyr::select(id_1, id_2, sib, Cohort_gap)

meta

#Calculate the gap between cohorts, and add it into a new dataframe

kinNWdata$Cohort_gap[kinNWdata$Cohort_gap<0] <- 1 # this replaces all of our "negative" cohort gaps so that we have a tidy plot

kinNWdata

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE) 

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, Cohort)
vertices <- as_data_frame(vertices, what = "vertices")

vertices

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)

kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = Cohort_gap,
        edge_colour = factor(sib)),
    #arrow = arrow(length = unit(3, 'mm')), 
    #end_cap = ggraph::circle(2, 'mm'),
    edge_alpha = 1) +
  ggraph::scale_edge_width(range = c(1,2), breaks = c(0,1), name = "Cohort Gap") +
  ggraph::scale_edge_colour_manual(values = c("orange", "skyblue"),
                                   name = "Kin Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(shape = sex),
                          size = 4) +
  ggraph::geom_node_text( aes(label = df$id), repel = TRUE, 
                          size = 5, color = "black") +
  ggplot2::scale_color_manual(values = adegenet::funky(9), na.value = "grey50") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = unit(rep(1,4), "cm")) 


print(kin_network1)


#### TRYING THE SIMULATIONS ### 

freqs <- input$freqs

simfreqs <- sample(freqs, 100, replace=F)
View(simfreqs)

View(input)

### Simulate relatedness values ### 

sim <- familysim(simfreqs, 100)

output <- coancestry(sim, wang=1)

#sim_rel <- coancestry(genotype.data = sim, 
                     trioml = 2, 
                     wang = 2, 
                     allow.inbreeding = T,
                     ci95.num.bootstrap = 100,
                     rng.seed = 42)

sim_rel <- cleanuprvals(output, 100)

## First we need to subsample our larger dataset to 100 loci (due to an internal error associated with related package)

subsample <- gl.subsample.loc(gl, 100, replace = F, v = 3, error.check = T)

## Then we convert this subset into a related compatible text file 

n100_geno <- gl2related(subsample, 
                   outfile = "subsample_related.txt", , 
                   outpath = "C:/Users/samue/OneDrive/Desktop/Honours/analysis",
                   v = 5)

## Load it back into R

sim_input <- readgenotypedata("subsample_related.txt")

## And compare the relatedness estimators to see how they go

comp_est <- compareestimators(sim_input, 100)

comp_est


sim <- familysim(sim_input$freqs, 100)


