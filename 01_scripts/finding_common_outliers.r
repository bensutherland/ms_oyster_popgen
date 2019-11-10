### This appears to be a very preliminary script

files <- list.files(path = "13_selection/", pattern = "outliers_common")

temp.df  <- NULL
for(i in 1:length(files)){
  
  temp <- read.delim(paste0("13_selection/", files[i]), sep = "\t", header = T)
  temp.df <- rbind(temp.df, temp)
  
}

head(temp.df)

table(temp.df$matching.mname) > 1
