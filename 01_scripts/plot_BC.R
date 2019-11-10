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
# devtools::install_github("dkahle/ggmap")
# devtools::install_github("hadley/ggplot2")
library("ggmap")
library("ggplot2")


#### 1. Manually Input BC GPS Locations ####
## Alternate automated method:
# loc.df <- geocode(c("Pendrell Sound", "Boundary Bay")) # no longer working without registration

PIP <- c(-125.291863, 49.022070)
HIS <- c(-126.503, 49.7329)
PEN <- c(-124.723141, 50.251223)
SER <- c(-122.922485, 49.018080)
ROS_DPB <- c(-124.766629, 49.466708)
# DPB <- c(-124.743098, 49.466778)

loc.df <- rbind(PIP, HIS, PEN, SER, ROS_DPB)
rownames(loc.df) <- c("PIP", "HIS", "PEN", "SER", "ROS/DPB")
loc.df <- as.data.frame(x = loc.df, stringsAsFactors = F)

colnames(loc.df) <- c("lon", "lat")

# Extract relevant info
visit.x <- loc.df$lon
visit.y <- loc.df$lat


#### 2. Plotting ###
# Plot a map of BC
# install.packages("maps")
library(maps)
# install.packages("mapdata")
library(mapdata)

# Plot locations
pdf(file = "BC_map.pdf")
map("worldHires","Canada", xlim=c(-130,-122), ylim=c(47.5,52.7), col="gray70", fill=TRUE)  #plot the region of Canada I want
map("worldHires","usa", xlim=c(-130,-122),ylim=c(47.5,52.7), col="gray70", fill=TRUE, add=TRUE)  #add the adjacent parts of the US; can't forget my homeland
points(loc.df$lon, loc.df$lat, pch=19, col="red", cex=1)  #plot my sample sites
text(x = loc.df$lon, y = loc.df$lat, labels =  rownames(loc.df), adj = c(-0.02,0.3), cex = 1.2)
text(x = -126.4, y = 50.2, labels = "Vancouver Island", cex = 0.9)
text(x = -124.7, y = 51.5, labels = "Mainland BC", cex = 0.9)
text(x = -129, y = 48.5, labels = "Pacific Ocean", cex = 0.9)
text(x = -123.1, y = 49.25, labels = "Vancouver", cex = 0.9)


dev.off()
