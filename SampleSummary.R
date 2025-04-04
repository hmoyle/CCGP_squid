### Creating Sample Summary
### All CCGP samples that have been sequenced
### plus the 14 that have been prepped
# 2025-04-04

library(tidyverse)
metadata <- read_csv("data/SquidSequenced_plus_Prepped20250404.csv")

ind_depth1 <- read_delim("data/vcftools_stats/hc-pools/hc-CCGP-squid_stats.idepth", 
                         delim = "\t", 
                         col_names = c("ind", "nsites", 
                                       "depth"), 
                         skip = 1)

ind_depth2 <- read_delim("data/vcftools_stats/lc-pools/lc-CCGP-squid_stats.idepth", 
                         delim = "\t", 
                         col_names = c("ind", "nsites", 
                                       "depth"), 
                         skip = 1)

ind_depth <- ind_depth1 %>%
  add_row(., ind_depth2)

data <- metadata %>% 
  rename(NMFS_DNA_ID = "NMFS_DNA_ID...6") %>%
  select(NMFS_DNA_ID, PoolName, LibraryName, Sequencer, 
         STATE_M, COLLECTION_DATE, CRUISE, 
         HAUL, LATITUDE_M, LONGITUDE_M) %>%
  left_join(., 
            ind_depth, 
            by = c("NMFS_DNA_ID" = "ind")) %>%
  select(-nsites)

write_csv(data, "data/SquidSamplesSummary_20250404.csv")



