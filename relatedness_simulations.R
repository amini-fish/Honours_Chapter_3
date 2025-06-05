## Analysing the simulated data 

## Load in the relevant packages

library(dartRverse)
library(related)
library(tidyverse)
library(dplyr)

## Set the WD

setwd("C:/Users/samue/Desktop/Honours/Relatedness_simulations")

## Load in simulated data from COANCESTRY

sims_df <- read.csv("simulated_rel.csv", stringsAsFactors = T) ## We use this as it has EMIBD9 EM1 results in it too 

sims_df

sims2 <- read.csv("simulated_rel.csv", stringsAsFactors = T)

sims2$True_rel <- sims2$True_rel/2

Em1 <- cor(sims2$EMIBD9, sims2$True_rel)
round(Em1, 4)
Em2 <- cor(sims2$EM2_EMIBD9, sims2$True_rel)

round(Em2, 3)

## NOTE - remember that you can just substitute co_sims and sims_df :)) 

sims_df$EMIBD9 <- sims_df$EMIBD9*2

sims_df$EM2_EMIBD9 <- sims_df$EM2_EMIBD9*2

sims_df <- sims_df %>%
  dplyr::rename("EM1" = "EMIBD9", "EM2"= "EM2_EMIBD9")

sims_df

sims_df <- sims_df %>% 
  pivot_longer(cols = c(TrioML, LynchRd, Ritland, QuellerGt, EM1, EM2), 
               names_to = "Estimator", 
               values_to = "rel")

sims_df <- sims_df %>%
  filter(Estimator != c("DyadML","LynchLi", "Wang"))


sims_df$theta <- sims_df$rel/2
  
sims_df$Sibtype <- factor(sims_df$Sibtype, levels =  c("FSP", "HSP", "FCP", "HFCP", "UP"))


sims_df$Estimator <- factor(sims_df$Estimator, levels = c("LynchRd", "QuellerGt", "Ritland", "TrioML", "EM1", "EM2"))

levels(sims_df$Sibtype)

True_rel <- data.frame(
  Sibtype = factor(c("FSP", "HSP", "FCP", "HFCP", "UP"), levels = c("FSP", "HSP", "FCP", "HFCP", "UP")),
  yintercept = c(0.25, 0.125, 0.0625, 0.03125, 0)
  )


## now in approrpriate format to plot...

plot1 <- ggplot(data = sims_df, aes(x = Estimator, y = theta, fill = Estimator)) +
  geom_boxplot(alpha = 0.5)+
  theme_bw()+
  scale_fill_brewer(palette = "RdYlGn") 

p1 <- plot1 + 
  facet_wrap(~Sibtype)

p1 <- p1 + 
  ggtitle("Simulated Relatedness Estimates for Relevant Kin-Types") +
  ylab("Relatedness")+
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        strip.text = element_text(size = 12), 
        axis.text = element_text(size = 11, angle = 90), 
        axis.title = element_text(size = 13)) + 
  geom_hline(data = True_rel, aes(yintercept = yintercept), linetype = "dotted", linewidth = 0.8, alpha = 0.7, col = "black")

print(p1)  

dat_text <- data.frame(
  Estimator = factor(c("LynchRd", "QuellerGt", "Ritland", "TrioML", "EM1", "EM2")),
  label = c(0.9941, 0.9937, 0.9918, 0.9958, 0.9955, 0.9946), 
  x = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0))  # X position for text

levels(dat_text$Estimator)
levels(sims_df$Estimator)

base_plot <- ggplot(sims_df, aes(x = theta, fill = Sibtype)) +
  geom_density(alpha = 1, position = "identity")+
  theme_bw()+
  ggtitle("Distribution of simulated relatedness coefficients") +
  xlab("Kinship") +
  ylab("Frequency")+
  facet_wrap(~Estimator) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 20, 0))) +
  geom_vline(data = True_rel, aes(xintercept = yintercept), linetype = "dotted", linewidth = 0.7, alpha = 0.9, col = "black")

print(base_plot)

tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}

my_tag <- c("r = 0.994", "r = 0.994", "r = 0.992", "r = 0.996", "r = 0.992", "r = 0.996", "r = 0.995")

base_plot <- tag_facet(base_plot,
    x = 0.18, y = 110, 
          vjust = 1, hjust = -0.25,
          open = "", close = "",
          size = 3,
          tag_pool = my_tag)
 
my_tag2 <- c("RMSE = 0.020", "RMSE = 0.020", "RMSE = 0.023", "RMSE = 0.0170", "RMSE = 0.030", "RMSE = 0.029", "RMSE = 0.018")

plot2 <- tag_facet(base_plot,
          x = 0.17, y = 120, 
          vjust = 1, hjust = -0.25,
          open = "", close = "",
          size = 3,
          tag_pool = my_tag2)

print(plot2)

#tiff(filename = "plot2.tiff",
     width = 25, height = 15, units = "cm", 
     compression = c("lzw"),
     bg = "white", res = 1500)

plot2

dev.off()
##------------------------------------------------------------------------------

## Calculating RMSE 

#rmse(data$actual, data$predicted) for structure of commands


my_rmse <- function(obs, pred) {
  sqrt(mean((obs-pred)^2))
}

rmse_em1 <- my_rmse(sims_df$EM1, sims_df$True_rel); rmse_em1

rmse_em2 <- my_rmse(sims_df$EM2, sims_df$True_rel); rmse_em2

rmse_trioml <- my_rmse(sims_df$TrioML, sims_df$True_rel); rmse_trioml

rmse_wang <- my_rmse(sims_df$Wang, sims_df$True_rel); rmse_wang

rmse_LyLi <- my_rmse(sims_df$LynchLi, sims_df$True_rel); rmse_LyLi

rmse_LyRd <- my_rmse(sims_df$LynchRd, sims_df$True_rel); rmse_LyRd

rmse_Qgt <- my_rmse(sims_df$QuellerGt, sims_df$True_rel); rmse_Qgt

rmse_rit <- my_rmse(sims_df$Ritland, sims_df$True_rel); rmse_rit

RMSE_all <- rbind(rmse_em1, rmse_em2, rmse_LyRd, rmse_Qgt, rmse_rit, rmse_trioml, rmse_wang)

RMSE_all <- data.frame(RMSE_all[order(RMSE_all), ])

colnames(RMSE_all) <- c("RMSE"); RMSE_all


#### EMIBD Em2 estimated #### 

glimpse(sims_df)

FS <- sims_df %>%
  filter(Sibtype == "FSP", Estimator == "TrioML") %>%
  select(rel) %>%
  droplevels()


HS <- sims_df %>%
  filter(Sibtype == "HSP", Estimator == "TrioML") %>%
  select(rel) %>%
  droplevels()


# Compute 95% confidence interval
ci_FS <- quantile(FS$rel, probs = c(0.01, 0.99))
ci_HS <- quantile(HS$rel, probs = c(0.01, 0.99))
# Print results
print(ci_FS)
print(ci_HS)

# Estimate density functions
# Find the minimum FS value that doesn't overlap with HS
threshold <- (max(HS$rel) + min(FS$rel)) / 2

print(threshold)
