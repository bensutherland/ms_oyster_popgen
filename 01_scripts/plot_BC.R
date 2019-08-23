# This is a script to plot a map with points on it, then add labels to the plot, and connect points with lines
# specifically here for the Moore Oyster Project

#### 00. Front Matter ####
# Clean space
# rm(list=ls())

# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("/hdd/01_moore_oyster_project/stacks_workflow/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}

# Install libraries
library(devtools)
devtools::install_github("dkahle/ggmap")
devtools::install_github("hadley/ggplot2")
library("ggmap")
library("ggplot2")


#### 1. Manual Input Locations ####
# May need to redo the boundary bay site (probably do all by hand..)
pipestem <- c(-125.291863, 49.022070)
hisnit <- c(-126.503, 49.7329)
loc.df <- geocode(c("Pendrell Sound", "Boundary Bay"))
loc.df
loc.df <- rbind(loc.df, pipestem, hisnit)
rownames(loc.df) <- c("PEN", "HIS", "BDB", "PIP")
loc.df

visit.x <- loc.df$lon
visit.y <- loc.df$lat



#### 2. Plotting ###
myLocation <- c(lon = -125, lat = 50) # set the point of reference
# Collect google map
myMap <- get_map(location = myLocation, source = "google", maptype = "terrain", crop = FALSE, zoom = 7, color = "bw")
# Collect stamen map
myMap <- get_map(location = myLocation, source = "stamen", maptype = "terrain", crop = FALSE, zoom = 7, color = "bw")
# Plot using ggmap
ma <- ggmap(myMap)
ma

# Insert the plot points
ma <- ma + geom_point(aes(x=visit.x, y=visit.y), color="black", size=3) 
ma

# use ggrepel to add sample sizes per location
# possibly sample sizes and fish sizes
library(ggrepel)
