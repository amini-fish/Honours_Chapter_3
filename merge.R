library(tidyverse)

pairs <- read.csv("daly_rel_all.csv", stringsAsFactors=TRUE)
meta <- read.csv('Daly_meta.csv', stringsAsFactors=TRUE)

names(pairs)
names(meta)

# merge meta data for Birth_Cohort - need to do separately for id1 and id2
# same approach for other ID specific variables

## id1
pairs <- pairs %>%
  rename(id = id_1)  %>%
  distinct()

meta <- meta %>%
  select(id, Year_caught, Catch.Set_ID, billabong) %>%
  rename(Year_caught_id1 = Year_caught, 
         Catch.Set_ID1 = Catch.Set_ID, 
         Billabong_ID1 = billabong)

names(pairs)
names(meta)
      
pairs_meta_id1 <- left_join(pairs, meta, by=("id"))

head(pairs_meta_id1)

## id2
pairs_meta_id1 <- pairs_meta_id1 %>%
  rename(id_1 = id) %>%
  rename(id = id_2) 
names(pairs_meta_id1)

meta <- meta %>%
  select(id, 
         Year_caught_id1, 
         Catch.Set_ID1, 
         Billabong_ID1) %>%
  rename(Year_caught_id2 = Year_caught_id1, 
         Catch.Set_ID2 = Catch.Set_ID1, 
         Billabong_ID2 = Billabong_ID1)

names(meta)

pairs_meta_id12 <- left_join(pairs_meta_id1, meta, by=("id"))

head(pairs_meta_id12)

write.csv(pairs_meta_id12, "pairs_meta.csv")

# 8 levels if including NA's
summary(pairs_meta_id12$year_caught_both)

pairs_meta_id12



pairs_meta_id12 <- pairs_meta_id12 %>%
  unite("year_caught_both",
        as.character("Year_caught_id1"), as.character("Year_caught_id2"), 
        sep = ".",
        remove = FALSE, na.rm = FALSE) %>%
    mutate(year_caught_both=factor(year_caught_both))

pairs_meta_id12

pairs_meta_id12 <- pairs_meta_id12 %>%
  unite("billabong_both",
        as.character("Billabong_ID1"), as.character("Billabong_ID2"), 
        sep = "_",
        remove = FALSE, na.rm = FALSE) %>%
  mutate(billabong_both=factor(billabong_both))

pairs_meta_id12  
