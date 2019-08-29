# Method to detect outliers, developed by Anne-Laure Ferchaud (U Laval) and Ben Sutherland (DFO)
# Uses BayeScan methods as per http://cmpg.unibe.ch/software/BayeScan/index.html

#### Front Matter ####
# Clean space
# rm(list=ls())

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

## Info
# sessionInfo()

### Load libraries
source("/home/ben/programs/BayeScan2.1/R functions/plot_R.r") # load bayescan plotting functions

# Identify contrast files
# All BayeScan analyses should have a pop_map file
contrasts <- list.files(path = "13_selection/", pattern = "_pop_map.txt")
contrasts <- gsub(pattern = "\\_pop\\_map\\.txt", replacement = "", x = contrasts)

# Temp
contrasts <- contrasts[1:3]

#### 1. Read in results ####
data.FN <- paste0("13_selection/", contrasts[1], "/", contrasts[1], "_fst.txt")
data.df <- read.table(file = fst.FN, header = TRUE
                   #, row.names = TRUE
                   )
dim(data.df)

## e.g. ##
# # set nulls
# pcadapt.obj <- NULL; pcadapt.FN <- NULL; pcadapt.out <- NULL
# pcadapt.out.list <- list(); obj.name <- NULL; pcadapt.obj.list <- list()
# 
# for(i in 1:length(files)){
#   
#   # Select file
#   pcadapt.FN <- files[i]
#   print(paste0("Working on ", pcadapt.FN))
#   
#   # Create object name
#   obj.name <- gsub(pattern = ".*//", replacement = "", x = files[i])
#   obj.name <- gsub(pattern = "\\..*", replacement = "", x = obj.name)
#   
#   # Read in data
#   pcadapt.obj <- read.pcadapt(input = pcadapt.FN, type = "vcf")
#   
#   # Save input
#   pcadapt.obj.list[[obj.name]] <- pcadapt.obj
#   
#   # Run pcadapt
#   print("Running pcadapt with K = 20")
#   pcadapt.out <- pcadapt(input = pcadapt.obj, K = 20)
#   
#   # Save output
#   pcadapt.out.list[[obj.name]] <- pcadapt.out
#   
# }
# end e.g. ##

####  2:  work with BAYESCAN results ####
bayescan_plot <- plot_bayescan(res = data.FN, FDR = 0.01)

bayescan_plot <- plot_bayescan("bayes_input_fst.txt", FDR = 0.01)
bayescan_res <- read.table("bayes_input_fst.txt")
dim(bayescan_res)
SNPs <- read.table("snpkey.txt", header = F)
head(SNPs)
dim(SNPs)
colnames(SNPs) <- c("order", "mname")
annot_res <- cbind(SNPs$mname, bayescan_res)
head(annot_res)
colnames(annot_res) <- c("SNP", "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST")
head(annot_res)
# POST_PROB = 1 & Q_VALUE = 0 == 0.0001 
str(annot_res)
annot_res[annot_res$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 # replaces with floor (does not exist)

#Note still need to do 'BUNCH OF ROUNDING'

# #CHN_CHNF results####
# bayescan_plot_CHN_CHNF = plot_bayescan("CHN_CHNF_bayescan_fst.txt", FDR = 0.01)
# bayescan_res_CHN_CHNF= read.table("CHN_CHNF_bayescan_fst.txt")
# SNPb_CHN_CHNF=read.table("CHN_CHNF_bayescan_markers_dictionary.tsv",header=T)
# bayescan_CHN_CHNF=cbind(SNPb_CHN_CHNF$MARKERS, bayescan_res_CHN_CHNF)
# colnames(bayescan_CHN_CHNF)=c("SNP","POST_PROB","LOG10_PO","Q_VALUE","ALPHA","FST")





# attach(bayescan_CHN_CHNF)
# class(bayescan_CHN_CHNF$Q_VALUE)
# bayescan_CHN_CHNF$Q_VALUE <- as.numeric(bayescan_CHN_CHNF$Q_VALUE)
# bayescan_CHN_CHNF[bayescan_CHN_CHNF$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001

# BUNCH OF ROUNDING
# bayescan_CHN_CHNF$POST_PROB <- (round(bayescan_CHN_CHNF$POST_PROB, 4))
# bayescan_CHN_CHNF$LOG10_PO <- (round(bayescan_CHN_CHNF$LOG10_PO, 4))
# bayescan_CHN_CHNF$Q_VALUE <- (round(bayescan_CHN_CHNF$Q_VALUE, 4))
# bayescan_CHN_CHNF$ALPHA <- (round(bayescan_CHN_CHNF$ALPHA, 4))
# bayescan_CHN_CHNF$FST <- (round(bayescan_CHN_CHNF$FST, 6))

# ADD a column for selection grouping
bayescan_CHN_CHNF$SELECTION <- ifelse(bayescan_CHN_CHNF$ALPHA>=0&bayescan_CHN_CHNF$Q_VALUE<=0.05,"diversifying",ifelse(bayescan_CHN_CHNF$ALPHA>=0&bayescan_CHN_CHNF$Q_VALUE>0.05,"neutral","balancing"))
bayescan_CHN_CHNF$SELECTION<- factor(bayescan_CHN_CHNF$SELECTION)
levels(bayescan_CHN_CHNF$SELECTION)

positive_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="diversifying",]
neutral_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="neutral",]
balancing_CHN_CHNF <- bayescan_CHN_CHNF[bayescan_CHN_CHNF$SELECTION=="balancing",]
xtabs(data=bayescan_CHN_CHNF, ~SELECTION)
write.table(neutral_CHN_CHNF, "neutral_CHN_CHNF.txt", row.names=F, quote=F)
write.table(balancing_CHN_CHNF, "balancing_CHN_CHNF.txt", row.names=F, quote=F)
write.table(positive_CHN_CHNF, "positive_CHN_CHNF", row.names=F, quote=F)
detach(bayescan_CHN_CHNF)
