# Uses output data from previous in order to find top outliers
# 1. a) What outliers have at least two pieces of evidence per comparison? (how many per comparison?)
#    b) What of these were seen in more than one comparison?

datasets <- list.files("13_selection/", pattern = "all_data_output.txt")

dataset <- NULL; dataset.FN <- NULL; all_data.list <- list() ; input.df <- NULL; all_data <- NULL

# Bring in all data and evaluate evidence
for(i in 1:length(datasets)){
  
  # Load datasets
  dataset.FN <- datasets[i]
  dataset <- gsub(x = dataset.FN, pattern = "_.*", replacement = "") # build short name
  input.df <- read.table(file = paste0("13_selection/", dataset.FN), header = TRUE, sep = "\t", stringsAsFactors = F)
  colnames(input.df)
  
  # Retain only specific columns
  keep_cols <- c("matching.mname", "Fst", "chr", "pos", "pval", "qval", "top.pc", "pcs.of.interest", "is.top.pc", "outlier_loading", "fst_percentile", "fst_percentile_cutoff")
  input.df <- input.df[, keep_cols]
  
  # Calculate evidence per test
  input.df$high.fst <- input.df$Fst > input.df$fst_percentile
  input.df$sig.pcadapt <- input.df$qval < 0.01 & input.df$is.top.pc == TRUE # Need to pull the qval cutoff earlier and include here from the table, as per fst_percentile above #TODO#
  input.df$sig.rda <- !is.na(input.df$outlier_loading)
  
  # Per locus, per contrast, how many pieces of evidence?
  input.df$top.outliers.sum <- rowSums(x = input.df[,c("high.fst", "sig.pcadapt", "sig.rda")], na.rm = T)
  
  # Make vector indicating the contrast
  input.df$contrast <- dataset
  all_data.list[[dataset]] <- input.df 
  
  # Also put all into one big df
  all_data <- rbind(all_data, input.df)
  
  }



# See what we've got
colnames(all_data)
dim(all_data)
table(all_data$contrast) # how many markers per contrast

# Keep only those values where at least X evidence showed significance (per marker, per contrast)
num.of.ev <- 2
sig_all_data <- all_data[all_data$top.outliers.sum >= num.of.ev, ]
dim(sig_all_data)
head(sig_all_data)

# Find the number of times each marker was seen to be signif (>= X pieces of evidence); i.e., in how many contrasts
num_contrasts_per_sig_marker <- as.data.frame(table(sig_all_data$matching.mname))
head(num_contrasts_per_sig_marker)
colnames(num_contrasts_per_sig_marker) <- c("mname", "freq.contrasts")
head(num_contrasts_per_sig_marker)

sig_all_data <- merge(x = sig_all_data, y = num_contrasts_per_sig_marker, by.x = "matching.mname", by.y = "mname", all.x = TRUE)

head(sig_all_data)

# Order by num contrasts, mname
# sig_all_data <- sig_all_data[with(sig_all_data, order(sig_all_data$freq.contrasts, sig_all_data$top.outliers.sum, sig_all_data$Fst, decreasing = T)), ]
sig_all_data <- sig_all_data[with(sig_all_data, order(sig_all_data$freq.contrasts, sig_all_data$matching.mname, decreasing = T)), ]
head(sig_all_data)

# 

# Output filename
sig_all_data.FN <- paste0("13_selection/", "top_outliers_w_", num.of.ev, "_ev.csv")
write.csv(x = sig_all_data, file = sig_all_data.FN, quote = F, row.names = F)

## This file actually contains all of the information needed to consider 'top' outliers

### MORE INVESTIGATION BELOW, STILL IN IMPLEMENTATION ###

#### Ranking outliers #### (HERE TODAY)
# Could do some work adding columns that says how many contrasts the marker is significant for...
# Lets view the data for some general stats
# how many markers are seen more than once with three evidences?
sig_all_data_triple_ev <- sig_all_data[sig_all_data$top.outliers.sum == 3, ]
table(sig_all_data_triple_ev$matching.mname) # if all 1, then no repeat view

sig_all_data_double_ev <- sig_all_data[sig_all_data$top.outliers.sum == 2, ]
table(sig_all_data_double_ev$matching.mname) # if all 1, then no repeat view
table(sig_all_data_double_ev$matching.mname)==2 # if all 1, then no repeat view

# Which markers from those with at least 2 pieces of evidence are seen in multiple data
mnames_double <- sig_all_data_double_ev[duplicated(x = sig_all_data_double_ev$matching.mname), "matching.mname"]

temp <- sig_all_data_double_ev[sig_all_data_double_ev$matching.mname %in% mnames_double, ]
temp <- temp[with(temp, order(temp$matching.mname)), ]

write.csv(x = temp, file = "13_selection/outliers_w_two_evidence_seen_two_contrasts.csv", quote = F, row.names = F)

### FOR ALL
# Which markers from all are seen in multiple data
mnames_more_than_one <- sig_all_data[duplicated(x = sig_all_data$matching.mname), "matching.mname"]

temp2 <- sig_all_data[sig_all_data$matching.mname %in% mnames_more_than_one, ]
temp2 <- temp2[with(temp2, order(temp2$matching.mname)), ]

write.csv(x = temp2, file = "13_selection/outliers_w_any_evidence_seen_two_contrasts.csv", quote = F, row.names = F)



