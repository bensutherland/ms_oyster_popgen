# Plot differentiation statistics across genome
## REQUIRES: already run gwas_02_load_set_contrasts_info.r

# Summary # your data is here:
names(plotting_data)
head(plotting_data[[1]])
datatypes # datatypes

# Set variable
qval_plot_thresh <- 0.05
col_pcadapt <- "red"
col_rda <- "blue"
FST_percentile <- 0.99

##### Plot data #####
all_data <- NULL 
metric_plot.FN <- NULL; fst.outlier <- NULL; FST_percentile_line <- NULL

for(i in 1:length(datatypes)){
  
  print(paste0("Working on datatype ", datatypes[i]))
  all_data <- plotting_data[[datatypes[i]]]
  
  # Make sure in same format
  all_data$chr <- as.character(all_data$chr)
  
  head(all_data)
  
  # Plot
  metric_plot.FN <- paste0("13_selection/", datatypes[i], "_sel_genome_plot.pdf")
  pdf(file = metric_plot.FN, width = 10, height = 5)

  print(paste0("Plotting", datatypes[i]))
    
  # Plot
  plot(x = all_data$cum.pos, y = all_data$Fst, type = "n", las = 1
       , ylab = paste0("Fst (contrast ", datatypes[i], ")")
       , xaxt = "n", cex = 0.6, xlab = "Genomic Location (bp)"
       #, ylim = c(-0.02, 0.4)
       )
  
  # Add 0 line
  abline(h = 0, lty = 1)
  
  points(x = all_data$cum.pos, y = all_data$Fst)
  # Add RDA info
  points(x = all_data$cum.pos[!is.na(all_data$outlier_loading)], y = all_data$Fst[!is.na(all_data$outlier_loading)], pch = 8
         , col = col_rda, cex = 1.2
         ) # plot rda outliers
  # Add pcadapt info
  points(x = all_data$cum.pos[all_data$qval < qval_plot_thresh], y = all_data$Fst[all_data$qval < qval_plot_thresh], pch = 3
         , col = col_pcadapt, cex = 1.5
        ) # plot pcadapt outliers
  
  # Add outlier line for FST
  
  # OLD METHOD (2sd of mean)
  # fst.avg <- mean(all_data$Fst, na.rm = T)
  # fst.sd <- mean(all_data$Fst, na.rm = T)
  # fst.outlier <- fst.avg + (2*fst.sd)
  # END OLD METHOD
  
  # 95th percentile method
  FST_percentile_line <- quantile(x = all_data$Fst, probs = FST_percentile, na.rm = T)
  abline(h = FST_percentile_line, lty = 3)

  # Add dividers between chr
  abline(v = chr.info$cum.lengths, lty = 2)
  
  # Add chromosome names
  axis(side = 1, at = position.of.chr.label, labels= chr.info$chr.names, las=2) 

  dev.off()
  
  # # Metric comparison (old)
  # print(paste0("Metric comparison for ", datatypes[i]))
  # metric_comparison.FN <- paste0("13_selection/", datatypes[i], "_metric_comparison.pdf") 
  # pdf(file = metric_comparison.FN, width = 5, height = 5)
  # plot(all_data$Fst, all_data$neg.log10.qval)
  # dev.off()

}
