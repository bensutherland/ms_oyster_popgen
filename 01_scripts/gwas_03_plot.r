# Plot differentiation statistics across genome
## REQUIRES: already run gwas_02_load_set_contrasts_info.r
# or, object containing post this script

#### Front Matter ####
# Clean space
# rm(list=ls())

# Load libraries
library("ggplot2")


load(file = "13_selection/constant_and_contrasts_loaded.Rdata")

# Summary # your data is here:
names(plotting_data)
datatypes # datatypes
head(plotting_data[["fr"]])


# Set variables
qval_plot_thresh <- 0.01
col_pcadapt <- "red"
cex_pcadapt <- 0.5
pch_pcadapt <- 16
col_rda <- "blue"
cex_rda <- 1
pch_rda <- 16
FST_percentile <- 0.99

limit.to.top.pc <- TRUE # if only want to print top pc pcadapt vals

##### Set up plot data #####
all_data <- NULL 
metric_plot.FN <- NULL; fst.outlier <- NULL; FST_percentile_line <- NULL ; grob.list <- list(); common_ID.RDA_pcadapt <- NULL; rda_and_pcadapt_data <- NULL

for(i in 1:length(datatypes)){
  
  print(paste0("Working on datatype ", datatypes[i]))
  all_data <- plotting_data[[datatypes[i]]]
  
  # Make sure in same format
  all_data$chr <- as.character(all_data$chr)
  
  head(all_data)

  
  ##### Subset all_data to make it easier to plot ####
  # Note: this selection retains all markers, not just those with alignment positions
  # Identify data with a significant RDA value (i.e. has an outlier_loading shown)
  # TODO: this value goes above 0.05 ##
  
  rda_data <- all_data[!is.na(all_data$outlier_loading), ]
  dim(rda_data)
  
  # Identify data with a significant pcadapt value (i.e. below user set qval threshold)
  ### HERE, also limit if "is.top.pc == TRUE" #####
  
  if(limit.to.top.pc==TRUE){
    
    pcadapt_data <- all_data[all_data$qval < qval_plot_thresh & all_data$is.top.pc==TRUE, ]
    pcadapt_data <- pcadapt_data[!is.na(pcadapt_data$qval), ] # get rid of NAs also
    dim(pcadapt_data)
    
  } else if(limit.to.top.pc==FALSE){
   
    pcadapt_data <- all_data[all_data$qval < qval_plot_thresh, ]
    pcadapt_data <- pcadapt_data[!is.na(pcadapt_data$qval), ] # get rid of NAs also
    dim(pcadapt_data)
     
  }
  
  # Identify common markers in RDA AND pcadapt data
  common_ID.RDA_pcadapt <- intersect(x = rda_data$rda.id, y = pcadapt_data$rda.id)
  rda_and_pcadapt_data  <- all_data[which(all_data$rda.id %in% common_ID.RDA_pcadapt),] # collect the common data
  
  # Save out the common data (RDA, pcadapt)
  write.table(file = paste0("13_selection/outliers_common_pcadapt_rda_", datatypes[i], ".txt"), x = rda_and_pcadapt_data
              , quote = F, row.names = F, col.names = T, sep = "\t")
  
  # What is the FST percentile line level?
  FST_percentile # percentile chosen
  FST_percentile_line <- quantile(x = all_data$Fst, probs = FST_percentile, na.rm = T)
  
  
  # Note: your data is in rda_data, pcadapt_data and rda_and_pcadapt_data
  
  #### Plot the results ####
  genome_plot.FN <- paste0("13_selection/", datatypes[i], "_sel_genome_plot.pdf")
  pdf(file = genome_plot.FN, width = 10, height = 5)

  print(paste0("Plotting ", datatypes[i]))
  
  ## Actual Plotting
  
  datatype.grob <- ggplot() + 
    
    # Add horizontal line
    geom_hline(yintercept = 0, linetype = "solid") +
    
    # Plot all Fst vals
    geom_point(data = all_data, aes(x = cum.pos, y = Fst), size = 1, shape = 16, colour = "darkgrey") + 
    
    # Add pcadapt outliers
    geom_point(data =  pcadapt_data, aes(x = cum.pos, y = Fst), colour = "red", shape = 16, size = 1) +
  
    # Add RDA outliers
    geom_point(data = rda_data, aes(x = cum.pos, y = Fst), colour = "blue", size = 1) +
    
    # Add common outliers
    geom_point(data = rda_and_pcadapt_data, aes(x = cum.pos, y = Fst), colour = "black", size = 2, shape = 8) +
    
    
    geom_hline(yintercept = FST_percentile_line, linetype = "dashed", colour = "darkgray") +
  
    geom_vline(xintercept = chr.info$cum.lengths, linetype = "dashed", colour = "darkgray") +
    
    theme(axis.text.x=element_blank()) +  # remove xaxis labels
    
    labs(x = "Chromosome and position (bp)", y = paste0("Fst (", datatypes[i], ")")) +
    
    ylim(-0.04, max(all_data$Fst, na.rm = T) + 0.02) + 
  
    # geom_text(x = 15560284, y = 0.1, label = "test")
    
    annotate("text", x=position.of.chr.label, y=-0.04, label= chr.info$chr.names)
    
    # geom_label() # still need to do
  
    # # Add chromosome names
    # axis(side = 1, at = position.of.chr.label, labels= chr.info$chr.names, las=2)
  
    print(datatype.grob)
  
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

