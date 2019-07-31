# This is to generate summary statistics (e.g. samples per population)

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory
setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/")

pop_map_ret <- read.table(file = "01-info_files/population_map_retained.txt", stringsAsFactors = F)
colnames(pop_map_ret) <- c("indiv", "pop.code")
head(pop_map_ret)
str(pop_map_ret)

pop_map_ret$pop <- gsub(pattern = "_.*", replacement = "", x = pop_map_ret$indiv)

print("This is the number of samples in each population")
table(pop_map_ret$pop)


#### Set pop definitions ####
pop_types <- list()
pop_types[["naturalized.pops"]] <- c("CHN", "FRA", "HIS", "JPN", "PEN", "PIP", "SER")
pop_types[["naturalized_bc.pops"]] <- c("HIS", "PEN", "PIP", "SER")
pop_types[["nat_to_farm.pops"]] <- c("CHNF", "FRAF", "PENF")
pop_types[["hatch.pops"]] <- c("DPB", "GUR", "QDC", "ROS", "RSC")

##### How many samples in each pop type?#####
for(i in 1:length(pop_types)){
  # Select a pop type
  pop_types_selected <- pop_types[[i]]
  
  # Report on the number of samples in this pop type
  print(paste("The number of samples in", names(pop_types)[i]))
  print(table(pop_map_ret$pop %in% pop_types_selected)["TRUE"])
}


#### how many reads in each sample ? ####
reads_per_sample <- read.table(file = "04-all_samples/reads_per_sample_table.txt", stringsAsFactors = F)
colnames(reads_per_sample) <- c("indiv", "reads")

summary(reads_per_sample$reads)

#### NEXT TO DO: Why isn't the unique mappings present?


#### how many mappings in each sample? ####
# NOT CURRENTLY WORKING, NEED TO CALCULATE
mappings.data <- read.table(file = "04-all_samples/"
                            , col.names = c("sample", "mappings"))
mappings.data
