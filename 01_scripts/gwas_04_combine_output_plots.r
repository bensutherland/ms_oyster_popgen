# Combine pdfs

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


# Install packages
install.packages("cowplot")
library("cowplot")
# 
# genome_plots.files <- list.files("13_selection/", pattern = "*_genome_plot.pdf")
# selection_plots <- genome_plots.files[grep(pattern = "bc|ch|fr", perl = T, x = genome_plots.files)]
# domestication_plots <- genome_plots.files[grep(pattern = "bc|ch|fr", perl = T, x = genome_plots.files, invert = T)]
# 
# selection_plots.FN <- paste0("13_selection/", selection_plots)

load(file = "13_selection/grob_list.Rdata")

names(grob.list)


# Plot selection
pdf(file = "gwas_plot_selection.pdf", width = 5, height = 6)

plot_grid(grob.list[["bc"]], grob.list[["fr"]], grob.list[["ch"]]
          , nrow = 3)

dev.off()

# Plot domestication
pdf(file = "gwas_plot_domestication.pdf", width = 5, height = 10)
plot_grid(grob.list[["dpb"]], grob.list[["ros"]], grob.list[["gur"]], grob.list[["qdc"]], grob.list[["rsc"]]
          , nrow = 5)

dev.off()

