# This is a script to plot a map with points on it, then add labels to the plot, and connect points with lines
# specifically here for the Moore Oyster Project

# Clean space
# rm(list=ls())

setwd(dir = "~/Desktop")

# Install libraries
library("ggmap")
# install.packages("maptools")
library(maptools)
library(maps)
# install.packages("ggrepel")
library(ggrepel)

# Make Fst dataframe
comparisons <- c("PEN-ROS", "PEN-DPB", "PEN-GUR"
                 , "PEN-FRA", "FRA-JPN", "FRA-CHN"
                 , "CHN-JPN", "PEN-JPN"
                 , "CHN-QDC", "CHN-RSC")
length(comparisons)
fst.vals <- c(0.008, 0.058, 0.030
              , 0.003, 0.003, 0.017
              , 0.024, 0.006
              , 0.017, 0.007)
length(fst.vals)

comparisons_fst.df <- as.data.frame(cbind(comparisons, fst.vals))
comparisons_fst.df
rownames(comparisons_fst.df) <- comparisons_fst.df$comparisons
comparisons_fst.df

# todo:
#make second map, smaller one, with BC locations

#### 1. Identify locations ####
# Identify locations to put points (actual sampling sites)
visited <- c("Pendrell Sound", "France", "Qingdao", "Japan")
ll.visited <- geocode(visited)
rownames(ll.visited) <- visited # rename the rows
ll.visited
# pull out vectors for plotting purposes
visit.x <- ll.visited$lon
visit.y <- ll.visited$lat

#### 2. Starting mapping ####
# Using GGPLOT, plot the Base World Map
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mapNA <- borders("world", xlim = c(-200,200), ylim = c(-25, 80), colour="gray50", fill="gray50")
mapNA2 <- borders("world", xlim = c(-200,200), ylim = c(10, 80), colour="gray50", fill="gray50")
mp <- ggplot() +   mapWorld
ma <- ggplot() +   mapNA
ma2 <- ggplot() + mapNA2

test <- ggplot() + mapWorld
test
test <- test + coord_cartesian(xlim = c(-200, 200), ylim = c(-15, 80))

# NOT THIS WAY!
# test2 <- ggplot() + mapWorld
# test2
# test2 + scale_x_continuous(limits = c(-200, 200)) +
#   scale_y_continuous(limits = c(10, 80))

# use NA map only
#mp <- ma 
# mp <- ma2
mp <- test
mp

# Add location points on top
mp <- mp + geom_point(aes(x=visit.x, y=visit.y), color="black", size=3) 
mp

# ADD LABELS TO PLACES !
mp <- mp + geom_text_repel(aes(x = visit.x[1], y = visit.y[1], label = "PEN"), fill = "white", direction = "y", vjust = 1)
# vjust puts above, direction = y adjusts only on the vertical
mp <- mp + geom_text_repel(aes(x = visit.x[2], y = visit.y[2], label = "FRA"), fill = "white", direction = "y", vjust = 1)
mp <- mp + geom_text_repel(aes(x = visit.x[3], y = visit.y[3], label = "CHN"), fill = "white", direction = "y", vjust = 1)
mp <- mp + geom_text_repel(aes(x = visit.x[4], y = visit.y[4], label = "JPN"), fill = "white", direction = "x", hjust = 1)
mp

# Connect locations with single line through all sites 
mp <- mp + geom_line(aes(x=visit.x[1:3], y = visit.y[1:3]), col = "black")
mp

# Add additional connectors to connect 1) BC to Japan; 2) France to Japan; 3) China to Japan w/ loop
# France to Japan
dip.fr.x <- c(2.213749, 130, 133, 136, 138.252924)
dip.fr.y <- c(46.22764, 45, 44, 43, 36.20482)
mp <- mp + geom_line(aes(x=dip.fr.x, y = dip.fr.y), col = "black")
mp

# BC to Japan
dip.bc.x <- c(-124.723670, 130, 133, 136, 138.252924)
dip.bc.y <- c(50.29198, 15, 16, 17, 36.20482)
mp <- mp + geom_line(aes(x=dip.bc.x, y = dip.bc.y), col = "black")
mp

# China to Japan
dip.ch.x <- c(120.382609, 121, 123, 124, 138.252924)
dip.ch.y <- c(36.06711, 30, 30, 30, 36.20482)
mp <- mp + geom_line(aes(x=dip.ch.x, y = dip.ch.y), col = "black")
mp

#### 3. Draw Rectangle for farms
mp <- mp + annotate("rect", xmin=-200, xmax=-150, ymin = 25, ymax=50, alpha=0.6)
mp <- mp + annotate("rect", xmin=150, xmax=200, ymin = 25, ymax=50, alpha=0.6)
mp


# BC-based farms 
# Add geom_labels for farms
mp <- mp + geom_label(aes(x = -185, y = 45, label = "ROS"), fill = "white")
mp <- mp + geom_label(aes(x = -185, y = 37.5, label = "DPB"), fill = "white")
mp <- mp + geom_label(aes(x = -185, y = 30, label = "GUR"), fill = "white")
mp

# Farms box label
mp <- mp + annotate("text", x = -175, y = 52, label = "Farms")
mp

# China-based farms
# Add geom_labels for farms
mp <- mp + geom_label(aes(x = 185, y = 45, label = "QDC"), fill = "white")
mp <- mp + geom_label(aes(x = 185, y = 30, label = "RSC"), fill = "white")
mp

# Add label for farms
mp <- mp + annotate("text", x = 175, y = 52, label = "Farms")
mp


# Add connectors from farms to Pendrell
# Connect locations with single line through all sites
# BC farms
mp <- mp + geom_line(aes(x=c(-124.723670, -170), y = c(50.29198, 45)), col = "blue") # ROS
mp <- mp + geom_line(aes(x=c(-124.723670, -170), y = c(50.29198, 37.5)), col = "blue") # DPB
mp <- mp + geom_line(aes(x=c(-124.723670, -170), y = c(50.29198, 30)), col = "blue") # GUR
mp

# China farms
mp <- mp + geom_line(aes(x=c(120.382609, 170), y = c(36.06711, 45)), col = "blue")
# mp <- mp + geom_line(aes(x=c(120.382609, 170), y = c(36.06711, 37.5)), col = "blue")
mp <- mp + geom_line(aes(x=c(120.382609, 170), y = c(36.06711, 30)), col = "blue")
mp


#### Add wild fst labels manually ####
mp <- mp + geom_label(aes(x = -61.25496, y = 48.25981, label =  comparisons_fst.df["PEN-FRA", "fst.vals"]), fill = "lightgrey", size = 3) #bc-fr
mp <- mp + geom_label(aes(x = 70, y = 46, label = comparisons_fst.df["FRA-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #fr-jpn
mp <- mp + geom_label(aes(x = 6.8, y = 32.6, label = comparisons_fst.df["PEN-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #bc-jpn
mp <- mp + geom_label(aes(x = 61.3, y = 41.15, label = comparisons_fst.df["FRA-CHN", "fst.vals"]), fill = "lightgrey", size = 3) #fr-chn
mp <- mp + geom_label(aes(x = 121, y = 28, label = comparisons_fst.df["CHN-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #chn-jpn
mp

#### Add farm Fst labels manually ####
# BC FARMS
bc.farm.fst <- as.numeric(c(as.character(comparisons_fst.df["PEN-ROS", "fst.vals"]), as.character(comparisons_fst.df["PEN-DPB", "fst.vals"]), as.character(comparisons_fst.df["PEN-GUR", "fst.vals"])) )
mp <- mp + geom_label(aes(x = c(-160, -160, -160)
                          , y = c(46.5, 39, 31.5)
                          , label = bc.farm.fst), fill = "lightgrey", size = 3) #ROS, DPB, GUR
mp


# CHN FARMS
chn.farm.fst <- as.numeric( 
  c(as.character(comparisons_fst.df["CHN-QDC", "fst.vals"])
    , as.character(comparisons_fst.df["CHN-RSC", "fst.vals"])
    ))
mp <- mp + geom_label(aes(x = c(160,  160)
                          , y = c(45, 30)
                          , label = chn.farm.fst), fill = "lightgrey", size = 3) #QDC, RSC
mp

# Save out
pdf(file = "oyster_fst_world_map_figure2.pdf", width = 8, height = 5
    )
mp
dev.off()

# HERE BE DRAGONS

# #### Automated location identification: add Fst data ####
# # test <- mp + annotate("text", x = 61, y = 50.3, label = "Fst=0.02")
# # test
# ll.visited # is the object with the positions of the sample sites
# 
# # Write for loop to find the location to put the fst values
# location.lon <- NULL
# location.lat <- NULL
# location <- NULL
# fst.locs <- NULL
# all.fst.locs <- NULL
# 
# for(i in 2:nrow(ll.visited)){
#   location.lon <- ( (ll.visited[i, 1] - ll.visited[(i-1), 1]) / 2 ) + ll.visited[(i-1), 1]
#   location.lat <- ( (ll.visited[i, 2] - ll.visited[(i-1), 2]) / 2 ) + ll.visited[(i-1), 2]
#   fst.locs <- c(location.lon, location.lat)
#   
#   all.fst.locs <- rbind(fst.locs, all.fst.locs )
# }
# 
# colnames(all.fst.locs) <- c("lon","lat")
# all.fst.locs
# 
# # add hacked one
# all.fst.locs <- rbind(all.fst.locs, c(65, 10))
# 
# # add text?
# test <- mp + annotate("text", x = all.fst.locs[,1], y = (all.fst.locs[,2] + 5 )
#                       , label = c("Fst=0.024", "Fst=0.017", "Fst=0.003", "Fst=0.005"))
# test