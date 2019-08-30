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

genome_plots.files <- list.files("13_selection/", pattern = "*_genome_plot.pdf")
selection_plots <- genome_plots.files[grep(pattern = "bc|ch|fr", perl = T, x = genome_plots.files)]
domestication_plots <- genome_plots.files[grep(pattern = "bc|ch|fr", perl = T, x = genome_plots.files, invert = T)]

selection_plots.FN <- paste0("13_selection/", selection_plots)

ggdraw() + draw_image(selection_plots.FN[1], y = 0.5)

ggdraw() + 
  draw_image(selection_plots.FN[2], 0, 0.5, 0.5, 0.5)
  draw_image(selection_plots.FN[1], 0.5, 0, 0, 0.5) +
  draw_image(selection_plots.FN[2], 0, 0.5, 0.5, 0.5)