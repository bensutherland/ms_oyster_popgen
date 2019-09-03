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
metric_plot.FN <- NULL; fst.outlier <- NULL; FST_percentile_line <- NULL ; grob.list <- list()

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
  
  rda_data <- all_data[!is.na(all_data$outlier_loading), ]
  dim(rda_data)
  
  pcadapt_data <- all_data[all_data$qval < qval_plot_thresh, ]
  pcadapt_data <- all_data[!is.na(pcadapt_data$qval), ]
  dim(pcadapt_data)
  
  FST_percentile_line <- quantile(x = all_data$Fst, probs = FST_percentile, na.rm = T)
  
  
  
  ### NEW SECTION
  #library("ggplot2")
  
  datatype.grob <- ggplot() + 
    geom_point(data = all_data, aes(x = cum.pos, y = Fst), size = 1, shape = 1) + 
    
    # Add pcadapt outliers
    geom_point(data =  pcadapt_data, aes(x = cum.pos, y = Fst), colour = "red", shape = 11, size = 2) +
  
    # Add RDA outliers
    geom_point(data = rda_data, aes(x = cum.pos, y = Fst), colour = "blue", size = 2) +  # expect some to not have positions
    
    # Add horizontal line
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = FST_percentile_line, linetype = "dashed", colour = "darkgray") +
  
    geom_vline(xintercept = chr.info$cum.lengths, linetype = "dashed", colour = "darkgray") +
    
    theme(axis.text.x=element_blank()) +  # remove xaxis labels
    
    labs(x = "Chromosome and position (bp)", y = paste0("Fst (", datatypes[i], ")")) +
    
    ylim(-0.035, max(all_data$Fst, na.rm = T) + 0.02) + 
  
    # geom_text(x = 15560284, y = 0.1, label = "test")
    
    annotate("text", x=position.of.chr.label, y=-0.03, label= chr.info$chr.names)
    
    # geom_label() # still need to do
  
    # # Add chromosome names
    # axis(side = 1, at = position.of.chr.label, labels= chr.info$chr.names, las=2)
  
    print(datatype.grob)

  ### END NEW SECTION
  
  
    dev.off()
  
     grob.list[[datatypes[i]]] <- datatype.grob
  
  
  
  # # Metric comparison (old)
  # print(paste0("Metric comparison for ", datatypes[i]))
  # metric_comparison.FN <- paste0("13_selection/", datatypes[i], "_metric_comparison.pdf") 
  # pdf(file = metric_comparison.FN, width = 5, height = 5)
  # plot(all_data$Fst, all_data$neg.log10.qval)
  # dev.off()

}

save(grob.list, file = "13_selection/grob_list.Rdata")

