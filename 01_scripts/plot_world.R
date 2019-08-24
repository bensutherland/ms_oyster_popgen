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
# install.packages("ggrepel")
library(ggrepel)
## NEW
library(maps)
library(mapdata)


#### 0. Make FST dataframe manually ####
# Make Fst dataframe
comparisons <- c("PEN-ROS", "PEN-DPB", "PEN-GUR"
                 , "PEN-FRA", "FRA-JPN", "FRA-CHN"
                 , "CHN-JPN", "PEN-JPN"
                 , "CHN-QDC", "CHN-RSC")

fst.vals <- c(0.0085, 0.0604, 0.0288
              , 0.0024, 0.0048, 0.0191
              , 0.0241, 0.0060
              , 0.0196, 0.0087)

comparisons_fst.df <- as.data.frame(cbind(comparisons, fst.vals))
comparisons_fst.df
rownames(comparisons_fst.df) <- comparisons_fst.df$comparisons
comparisons_fst.df


#### 1. Identify locations ####
# Manually enter GPS
PEN <- c(-124.723141, 50.251223)
FRA <- c(-2.994716, 48.765167)
CHN <- c(120.267280, 36.038603)
JPN <- c(139.601050, 35.088636)

loc_global.df <- rbind(PEN, FRA, CHN, JPN)
rownames(loc_global.df) <- c("PEN", "FRA", "CHN", "JPN")
loc_global.df <- as.data.frame(x = loc_global.df, stringsAsFactors = F)

colnames(loc_global.df) <- c("lon", "lat")
loc_global.df

# Extract relevant info
visit.x <- loc_global.df$lon
visit.y <- loc_global.df$lat





#### 2. Starting mapping ####
# Option 1. maps package
## < probably too detailed for global, use ggplot below >
# pdf(file = "maps_version_global.pdf")
# map("worldHires"
#     #, regions = "."
#     #, xlim=c(-130,-122)
#     #, wrap = c(0,360)
#     , ylim=c(15, 75)
#     , col="gray70"
#     , fill=TRUE
# )  
# 
# points(loc_global.df$lon, loc_global.df$lat, pch=19, col="red", cex=1)  #plot my sample sites
# text(x = loc_global.df$lon, y = loc_global.df$lat, labels =  rownames(loc_global.df), adj = -0.2, cex = 1.2)
# dev.off()


# Option 2. ggplot
# Set variables
ylim.upper <- 80
ylim.lower <- 10
xlim.left <- -180
xlim.right <- 180

map_trimmed <- NULL ; map_plot <- NULL
map_trimmed <- borders(database = "world"
                 , xlim = c(xlim.left, xlim.right)
                 , ylim=c(ylim.lower, ylim.upper)
                 , colour="grey50"
                 , fill="gray70"
                 )

map_plot <- ggplot() + map_trimmed
map_plot

# Add location points
map_plot <- map_plot + geom_point(aes(x=visit.x, y=visit.y), color="black", size=2.5) 
map_plot


loc_global.df

map_plot

# Add location labels
map_plot <- map_plot + 
            geom_label_repel(
              # Add parameters
                direction = "x"
              , vjust = -2
              , hjust = 0
              
              # Add aesthetics
              , aes(x = loc_global.df$lon
                  , y = loc_global.df$lat
                  , label = rownames(loc_global.df)
                    )
              
            )

map_plot

## Add connection lines

# Connect locations with single line through all sites 
map_plot <- map_plot + geom_line(
                    aes(x=loc_global.df$lon[1:4]
                        , y = loc_global.df$lat[1:4]
                        )
                    , col = "black"
                    )
map_plot

# Connect different countries via indirect lines
contrasts.starts <- c("FRA", "PEN", "CHN")
contrasts.ends <- c("JPN", "JPN", "JPN")
contrasts.df <- as.data.frame(cbind(contrasts.starts, contrasts.ends), strinsAsFactors = FALSE)
contrasts.df$contrasts.starts <- as.character(contrasts.df$contrasts.starts)
contrasts.df$contrasts.ends <- as.character(contrasts.df$contrasts.ends)

# Add midpoints specific to each contrasts
contrasts.df$mid1.x <- c(130, 130, 121)
contrasts.df$mid2.x <- c(133, 133, 123)
contrasts.df$mid3.x <- c(136, 136, 124)
contrasts.df$mid1.y <- c(45, 15, 30)
contrasts.df$mid2.y <- c(44, 16, 30)
contrasts.df$mid3.y <- c(43, 17, 30)
contrasts.df
loc_global.df

# Add lines
START <- NULL ; END <- NULL; dip.x <- NULL; dip.y <- NULL; map_plots.temp <- NULL
dips.list <- list()
for(i in 1:nrow(contrasts.df)){
  
  
  print(i)
  
  # Set contrast
  START <- contrasts.df$contrasts.starts[i] 
  END <- contrasts.df$contrasts.ends[i] 
  
  dip.x <- c(loc_global.df$lon[rownames(loc_global.df)==START]
              
              # Take from the midpoints
              , contrasts.df$mid1.x[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
              , contrasts.df$mid2.x[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
              , contrasts.df$mid3.x[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
              
              , loc_global.df$lon[rownames(loc_global.df)==END]
              )
  dip.y <- c(loc_global.df$lat[rownames(loc_global.df)==START]
              
             # Take from the midpoints
             , contrasts.df$mid1.y[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
             , contrasts.df$mid2.y[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
             , contrasts.df$mid3.y[contrasts.df$contrasts.starts==START & contrasts.df$contrasts.ends==END]
             
             , loc_global.df$lat[rownames(loc_global.df)==END]
             
             )
  
  dips.list[[i]] <- cbind(dip.x, dip.y)
  
  }

# Add connectors to plot
map_plot <- map_plot + geom_line(aes(x = dips.list[[1]][,"dip.x"], y = dips.list[[1]][,"dip.y"]), col = "black")
map_plot <- map_plot + geom_line(aes(x = dips.list[[2]][,"dip.x"], y = dips.list[[2]][,"dip.y"]), col = "black")
map_plot <- map_plot + geom_line(aes(x = dips.list[[3]][,"dip.x"], y = dips.list[[3]][,"dip.y"]), col = "black")

map_plot





#### 3. Draw Rectangle for farms #####
mp <- map_plot # for troubleshooting
mp <- mp + annotate("rect", xmin=-210, xmax=-160, ymin = 15, ymax=40, alpha=0.6)
mp <- mp + annotate("rect", xmin=160, xmax=210, ymin = 15, ymax=40, alpha=0.6)
mp

# Add connectors from farms to Pendrell
# Connect locations with single line through all sites
# BC farms
mp <- mp + geom_line(aes(x=c(loc_global.df[1,1], -185), y = c(loc_global.df[1,2], 35)), col = "blue") # ROS
mp <- mp + geom_line(aes(x=c(loc_global.df[1,1], -185), y = c(loc_global.df[1,2], 27.5)), col = "blue") # DPB
mp <- mp + geom_line(aes(x=c(loc_global.df[1,1], -185), y = c(loc_global.df[1,2], 20)), col = "blue") # GUR
mp

# BC-based farms 
# Add geom_labels for farms
mp <- mp + geom_label(aes(x = -185, y = 35, label = "ROS"), fill = "white")
mp <- mp + geom_label(aes(x = -185, y = 27.5, label = "DPB"), fill = "white")
mp <- mp + geom_label(aes(x = -185, y = 20, label = "GUR"), fill = "white")
mp

# Add connectors from farms to CHN
mp <- mp + geom_line(aes(x=c(loc_global.df[3,1], 185), y = c(loc_global.df[3,2], 35)), col = "blue")
mp <- mp + geom_line(aes(x=c(loc_global.df[3,1], 185), y = c(loc_global.df[3,2], 20)), col = "blue")
mp

# China-based farms
# Add geom_labels for farms
mp <- mp + geom_label(aes(x = 185, y = 35, label = "QDC"), fill = "white")
mp <- mp + geom_label(aes(x = 185, y = 20, label = "RSC"), fill = "white")
mp

# Add label for farms
mp <- mp + annotate("text", x =  185, y = 42, label = "Farms")
mp <- mp + annotate("text", x = -185, y = 42, label = "Farms")
mp


#### Add wild fst labels manually ####
mp <- mp + geom_label(aes(x = -61.25496, y = 48.25981, label =  comparisons_fst.df["PEN-FRA", "fst.vals"]), fill = "lightgrey", size = 3) #bc-fr
mp <- mp + geom_label(aes(x = 70, y = 46, label = comparisons_fst.df["FRA-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #fr-jpn
mp <- mp + geom_label(aes(x = 6.8, y = 32.6, label = comparisons_fst.df["PEN-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #bc-jpn
mp <- mp + geom_label(aes(x = 61.3, y = 40, label = comparisons_fst.df["FRA-CHN", "fst.vals"]), fill = "lightgrey", size = 3) #fr-chn
mp <- mp + geom_label(aes(x = 121, y = 28, label = comparisons_fst.df["CHN-JPN", "fst.vals"]), fill = "lightgrey", size = 3) #chn-jpn
mp

#### Add farm Fst labels manually ####
# BC FARMS
bc.farm.fst <- as.numeric(c(as.character(comparisons_fst.df["PEN-ROS", "fst.vals"]), as.character(comparisons_fst.df["PEN-DPB", "fst.vals"]), as.character(comparisons_fst.df["PEN-GUR", "fst.vals"])) )
mp <- mp + geom_label(aes(x = c(-152, -152, -152)
                          , y = c(44.5, 37, 29.5)
                          , label = bc.farm.fst), fill = "lightgrey", size = 3) #ROS, DPB, GUR
mp

# CHN FARMS
chn.farm.fst <- as.numeric(c(as.character(comparisons_fst.df["CHN-QDC", "fst.vals"])
                , as.character(comparisons_fst.df["CHN-RSC", "fst.vals"])
                ))
mp <- mp + geom_label(aes(x = c(155,  155)
                          , y = c(35, 27)
                          , label = chn.farm.fst), fill = "lightgrey", size = 3) #QDC, RSC
mp

# Save out
pdf(file = "oyster_fst_world_map_figure2.pdf", width = 8.5, height = 4.5
    )
mp
dev.off()
