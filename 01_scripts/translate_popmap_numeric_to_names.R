# Script to transform the numeric code to str code to produce the population_map_retained file 

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


# Load pop_map
pop_map.FN <- "01-info_files/population_map_retained_numeric.txt"
print(paste0("Loading population map from: ", pop_map.FN))
pop_map <- read.table(file = pop_map.FN)
colnames(pop_map) <- c("sample","numeric.id")
head(pop_map)

# Load pop codes
pop_codes.FN <- "./../ms_oyster_popgen/00_archive/pop_codes.txt"
print(paste0("Loading pop codes from: ", pop_codes.FN))
pop_codes <- read.delim2(file = pop_codes.FN, header = F, sep = "\t")
colnames(pop_codes) <- c("pop.id", "numeric.id")
pop_codes

# Join pop_map and pop_codes
print("Combining map with pop codes")
pop_map_renamed <- merge(x = pop_map, y = pop_codes, all.x = T, by = "numeric.id")
head(pop_map_renamed)

# Prep to use as the new population map
pop_map_renamed <- pop_map_renamed[, c("sample", "pop.id")]

head(pop_map_renamed)

write.table(x = pop_map_renamed, file = "01-info_files/population_map_retained_names.txt"
            , quote = F, sep = "\t", col.names = F, row.names = F)
