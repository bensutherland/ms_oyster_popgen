# Fix population map to avoid small groups being separated

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}


# Load pop map
pop_map.FN <- "01-info_files/population_map_retained.txt"
print(paste0("Loading population map from: ", pop_map.FN))
pop.map <- read.table(file = pop_map.FN)
colnames(pop.map) <- c("sample","pop")
head(pop.map)
dim(pop.map) # this should only include the samples retained in the analysis

# Rename all BC as pop 1 (not currently used)
# pop.map$pop[pop.map$pop==1] <- 1
# pop.map$pop[pop.map$pop==2] <- 1
# pop.map$pop[pop.map$pop==6] <- 1
# pop.map$pop[pop.map$pop==7] <- 1
# pop.map$pop[pop.map$pop==8] <- 1

# Rename all France as pop 14
pop.map$pop[pop.map$pop==14] <- 14
pop.map$pop[pop.map$pop==15] <- 14

# Keep Rosewall (4)

# Keep Guernsey (9)

# Rename the two wild china as 12
pop.map$pop[pop.map$pop==12] <- 12
pop.map$pop[pop.map$pop==13] <- 12
# Keep ChinaRSC (10)
# Keep ChinaQDC (11)

# Rename Pendrell and PendrellFarm as 2
pop.map$pop[pop.map$pop==6] <- 2
pop.map$pop[pop.map$pop==2] <- 2

table(pop.map$pop)

# total number alleles
sum(table(pop.map$pop)) * 2 * 0.01 
# 0.01 is totally reasonable as a maf global. For a sample size of 352, 
# ...means need to see allele 7 times (or at least in four individuals)

write.table(x = pop.map, file = "01-info_files/population_map_retained_renamed.txt", sep = "\t", col.names = F, row.names = F, quote = F)
