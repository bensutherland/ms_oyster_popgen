# Identifies which outliers are found with more than one method: 

data.files <- list.files(path = "13_selection/", pattern = "_all_data_output.txt")

selection_results.list <- list() ; datatype <- NULL

for(i in 1:length(data.files)){

  # Identify datatype
  datatype <- gsub(x = data.files[i], pattern = "_.*", replacement = "") 
    
  # Read in data
  selection_results.list[[datatype]] <- read.table(file = paste0("13_selection/", data.files[i]), header = T, sep = "\t")
  
}


head(selection_results.list[["bc"]])

# data[!is.na(data$outlier_loading) & data$qval < 0.01, ]

concordance_results.list <- list() ; datatype <- NULL ; data <- NULL ; rows_of_interest <- NULL; output.FN <- NULL

# Characterize
for(j in 1:length(selection_results.list)){
  
  datatype <- names(selection_results.list)[j]
  
  data <- selection_results.list[[datatype]]

  # Identify which rows contain the data of interest
  t1 <- which(!is.na(data$outlier_loading))
  t2 <- which(data$qval < 0.01)
  
  rows_of_interest <- data[intersect(t1, t2),]
  
  concordance_results.list[[datatype]] <- rows_of_interest
  
  # output results
  output.FN <- paste0("13_selection/selection_concordance_", datatype, ".csv")
  write.csv(output.FN, x = concordance_results.list[[datatype]])

}

