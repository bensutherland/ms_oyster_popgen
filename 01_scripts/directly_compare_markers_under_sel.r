### Directly compare data without any genome information ####
datatypes <- c("bc", "fr", "ch")
perloc_stats_list <- list()
comps <- c("bc_fr", "bc_ch", "fr_ch")

input_dir <- "11-adegenet_analysis/"
output_dir <- "13_selection/"

percentile <- 0.95

# Loop
for(i in 1:length(datatypes)){
  datatype <- datatypes[i]
  print(datatype)
  filename <- paste(input_dir, datatype, "_sel_perloc_stats.csv", sep = "")

  # Selected data
  perloc_stats <- read.csv(file = filename)
  perloc_stats_list[[datatype]]  <- perloc_stats
}

str(perloc_stats_list)

pdf(file = paste0(output_dir, "selection_expt_all_by_all.pdf"), width = 8, height = 8)
par(mfrow=c(2,2))


# Loop to plot
pop1 <- NULL; pop2 <- NULL
for(i in 1:length(comps)){
  
  # Identify the two comparisons
  pop1 <- gsub(x = comps[i], pattern = "_.*", replacement = "")
  pop2 <- gsub(x = comps[i], pattern = ".*_", replacement = "")
  
  # Report
  print(paste0("Working on ", pop1, " vs. ", pop2))
  
  # Plot
  plot(perloc_stats_list[[pop1]]["Fst"]$Fst ~ perloc_stats_list[[pop2]]["Fst"]$Fst
       , ylim = c(-0.05, 0.3)
       , xlim = c(-0.05, 0.3)
       , ylab = paste0(toupper(pop1), " (Fst)")
       , xlab = paste0(toupper(pop2), " (Fst)")
       , las = 1
       )
  
  # Find the 95th percentile of the Fst for each pop
  pop1.quantile <- quantile(x = perloc_stats_list[[pop1]]["Fst"]$Fst, na.rm = TRUE, probs = percentile)
  pop2.quantile <- quantile(x = perloc_stats_list[[pop2]]["Fst"]$Fst, na.rm = TRUE, probs = percentile)
  abline(h = pop1.quantile, lty = 3)
  abline(v = pop2.quantile, lty = 3)
  
  # How many markers are in the top-right quadrant?
  table(perloc_stats_list[[pop1]]["Fst"]$Fst > pop1.quantile)
  table(perloc_stats_list[[pop2]]["Fst"]$Fst > pop2.quantile)
  
  shared_true <- table(perloc_stats_list[[pop1]]["Fst"]$Fst > pop1.quantile & perloc_stats_list[[pop2]]["Fst"]$Fst > pop2.quantile)["TRUE"] # shared
  
  text(x = 0.2, y = 0.25, labels = paste0("markers in n.", (percentile*100), "= ", shared_true))
  
}

dev.off()





# # Find out which markers correspond to this:
# perloc_ch.gid <- as.data.frame(perloc_stats_list$ch.gid)
# perloc_fr.gid <- as.data.frame(perloc_stats_list$fr.gid)
# perloc_bc.sel.gid <- as.data.frame(perloc_stats_list$bc.sel.gid)
# 
# FR <- perloc_fr.gid[which(perloc_fr.gid$Fst > 0.1), "X"]
# BC <- perloc_bc.sel.gid[which(perloc_bc.sel.gid$Fst > 0.1), "X"]
# 
# intersect(x = FR, y = BC) # tells you which markers
# 
# 
# 
# 
