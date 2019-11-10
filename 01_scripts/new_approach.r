# New approach
## ATTEMPT TO ADD UNKN MARKERS
#plotting_data.bck <- plotting_data

# Where is the chr info?
chr.info
# what is the avg length of a chr
mean_chr_length <- mean(chr.info$chr.lengths)

names(plotting_data)

temp_data <- NULL; temp_plotting_data.list <- list()
for(i in 1:length(names(plotting_data))){
  
  # Extract to work on the plotting data
  temp_data <- plotting_data[[i]]
  temp_data <- temp_data[order(temp_data$cum.pos), ]
  
  head(temp_data)
  
  
  ## ATTEMPT TO ADD UNKN MARKERS
  # # what is the average number of markers per chr
  # mean_markers_per_chr <- mean(table(temp_data# names(temp_plotting_data.list)
  # tail(temp_plotting_data.list[["bc"]])
  # 
  # 
  # temp_data2 <- temp_plotting_data.list[["bc"]]
  # head(temp_data2)
  # unkn.chr.left <- min(temp_data2[is.na(temp_data2$chr), "cum.pos"])
  # unkn.chr.right <- max(temp_data2[is.na(temp_data2$chr), "cum.pos"])
  # 
  # 
  # # known.chr.left <- min(temp_data2[!is.na(temp_data2$chr), "cum.pos"])
  # known.chr.right <- max(temp_data2[!is.na(temp_data2$chr), "cum.pos"])
  # $chr))
  # marker_distance <- mean_chr_length / mean_markers_per_chr
  # 
  # # Add another mean chr length distance to cum.pos
  # 
  # for(j in 1:nrow(temp_data)){
  #   
  #   # If the cumulative position is NA
  #   if(is.na(temp_data$cum.pos[j])){
  #     
  #     # Add the marker_distance to the first NA marker
  #     temp_data$cum.pos[j] <- (max(temp_data$cum.pos, na.rm = TRUE) + marker_distance)
  #     
  #   } else {
  #     
  #     #print("Already has position")
  #   
  #   }
  #   
  # }
  # 
  
  unique(temp_data$chr)
  odd.chr <- unique(temp_data$chr)[seq(from = 1, to = (length(unique(temp_data$chr))-1), by = 2)]
  even.chr <- unique(temp_data$chr)[seq(from = 2, to = (length(unique(temp_data$chr))-1), by = 2)]
  
  temp_data$marker.colour[temp_data$chr %in% odd.chr] <- "darkgrey"
  temp_data$marker.colour[temp_data$chr %in% even.chr] <- "lightgrey"
  
  head(temp_data)
  
  # Add back into a list
  temp_plotting_data.list[[names(plotting_data[i])]] <- temp_data
  
}


plotting_data <- temp_plotting_data.list 

names(plotting_data)
head(plotting_data[["bc"]])


 ### WORKING ###
# names(temp_plotting_data.list)
# tail(temp_plotting_data.list[["bc"]])
# 
# 
# temp_data2 <- temp_plotting_data.list[["bc"]]
# head(temp_data2)
# unkn.chr.left <- min(temp_data2[is.na(temp_data2$chr), "cum.pos"])
# unkn.chr.right <- max(temp_data2[is.na(temp_data2$chr), "cum.pos"])
# 
# 
# # known.chr.left <- min(temp_data2[!is.na(temp_data2$chr), "cum.pos"])
# known.chr.right <- max(temp_data2[!is.na(temp_data2$chr), "cum.pos"])
