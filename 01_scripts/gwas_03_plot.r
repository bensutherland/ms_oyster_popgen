# Plot differentiation statistics across genome
## REQUIRES: already run gwas_02_load_set_contrasts_info.r

# Summary # your data is here:
names(plotting_data)
head(plotting_data[[1]])
datatypes # datatypes

##### Plot data #####
###TODO # Should not use all_data, already a variable
all_data <- NULL 
metric_plot.FN <- NULL

for(i in 1:length(datatypes)){
  
  print(paste0("Working on datatype ", datatypes[i]))
  all_data <- plotting_data[[datatypes[i]]]
  
  # Make sure in same format
  all_data$chr <- as.character(all_data$chr)
  
  # Plot
  metric_plot.FN <- paste0("13_selection/", datatypes[i], "_sel_genome_plot.pdf")
  pdf(file = metric_plot.FN, width = 10, height = 10)
  par(mfrow=c(2,1))
  plot.types <- c("Fst", "neg.log10.qval")
  
  for(l in 1:length(plot.types)){
    
    print(paste0("Plotting", datatypes[i]))
    
    metric <- plot.types[l]
    
    plot(x = all_data$cum.pos, y = all_data[, metric], xlab = "cumulative pos (bp)", ylab = metric
         , las = 1
         #, ylim = c(-0.02, 0.4)
         , xaxt = "n"
         , cex = 0.6)
    
    abline(v = chr.info$cum.lengths)
    # abline(h = outlier_thresh_2xSD, lty = 2) # need to do this dynamically
    #text(x = position.of.chr.label, y = max(all_data[, metric], na.rm = T), labels = chr.info$chr.names)
    
    axis(side = 1, at = position.of.chr.label, labels= chr.info$chr.names, las=2) 
    
  }
  dev.off()
  
  print(paste0("Metric comparison for ", datatypes[i]))
  metric_comparison.FN <- paste0("13_selection/", datatypes[i], "_metric_comparison.pdf") 
  pdf(file = metric_comparison.FN, width = 5, height = 5)
  plot(all_data$Fst, all_data$neg.log10.qval)
  dev.off()

}
