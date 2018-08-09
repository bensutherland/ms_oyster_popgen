# Explore genetic diversity differences among hatchery and wild stock
# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_no_Qc/09-diversity_stats/")

#install.packages("rlist")
require("rlist")

#### 1. Nucleotide Diversity ####

# Import data
datatypes <- c("wild_bc", "rosewall", "guernsey", "france_wild", "france_cultured", "china", "bc_wild_to_farm")

datatypes.list <- list(); datatype <- NULL

for(i in 1:length(datatypes)){
  
  #identify the datatype name
  datatype <- datatypes[i]
  
  # identify the filename
  filename <- paste(datatype,"_samples_out.sites.pi", sep = "")
  
  # read it in
  datatypes.list[[datatype]] <- read.table(filename, header = T)
  
}

str(datatypes.list)


wild_bc <- as.data.frame(list.extract(datatypes.list, 1))
str(wild_bc)


# See average Pi per sample type
# wild_bc <- NULL; rosewall <- NULL; guernsey <- NULL; france_wild <- NULL; france_cultured <- NULL ; china <- NULL;  bc_wild_to_farm <- NULL
data <- NULL
data <-  as.data.frame(list.extract(datatypes.list, 1))[3]
summary(data$PI)

# Extract elements from list
for(i in 2:length(datatypes.list)){
  print(datatypes[i])
  print(summary(datatypes.list[[i]]$PI))
  
  data <- cbind(data, as.data.frame(list.extract(datatypes.list, i))[3])
  
}

str(data)
colnames(data) <- datatypes

str(data)
colnames(data)

# Plot
par(mfrow=c(3,2))
for(i in 2:length(datatypes)){
  plot(data$wild_bc ~ data[,i], xlab = colnames(data)[i])
}

par(mfrow=c(1,1), mar = c(3,4,3,3))
#str(data)
boxplot(data, las = 1, ylab = "Pi")

# RUN A MODEL, NEED TO SEPARATE BY HATCHERY / WILD SAMPLES

# #### 1.1. Pi per locus across types ####
# 
# pdf(file=paste("Pi_per_locus_by_type.pdf",sep=""), height=5, width=7)
# par(mfrow=c(1,3))
# plot(x$Pi.hatchery, ylab = "Pi (hatchery)", xlab = "Locus.index", ylim = c(0,1), las = 1)
# plot(x$`Pi.wild-farm`, ylab = "Pi (wild-to-farm)", xlab = "Locus.index", ylim = c(0,1), las = 1)
# plot(x$Pi.wild, ylab = "Pi (wild)", xlab = "Locus.index", ylim = c(0,1), las = 1)
# dev.off()
# 
# 
# par(mfrow=c(1,2))
# plot(x$Pi.qc, ylab = "Pi (qc)", xlab = "Locus.index", ylim = c(0,1), las = 1)
# plot(x$Pi.bc, ylab = "Pi (bc)", xlab = "Locus.index", ylim = c(0,1), las = 1)
# 
# #### 1.2 Pi correlation ####
# # How correlated are the pi values across loci for the different types (hatchery, wild and wild-farm)?
# #install.packages('corrplot')
# library('corrplot')
# cor.test(h$PI, w$PI) # correlated (cor = 0.52)
# cor.test(h$PI, wf$PI) # correlated (cor = 0.48)
# cor.test(wf$PI, w$PI) # highly correlated (cor = 0.84)
# 
# cor.set <- cor(x[,c("Pi.hatchery","Pi.wild-farm","Pi.wild")]
#                , use = "pairwise.complete.obs" # cor bw ea pair var is comp using all compl pairs of obs on those var # note: only works with "pearson" method    
# )
# 
# 
# # Plot correlation
# pdf(file=paste("corr_of_pi_bw_groups.pdf",sep=""), height=5, width=5)
# par(mfrow=c(1,1), mar = c(3,3,4,5))
# corrplot(cor.set, method="number", addshade="all", tl.col = "black", tl.cex=0.8, type = "lower", diag = T, tl.pos = "lt", cl.pos = "r")
# corrplot(cor.set, method="circle", addshade="all", tl.col = "black", tl.cex=0.8, type = "upper", add = T, tl.pos = "n")
# dev.off()

# #### 1.3. Pi differences across groups (model) ####
# dim(h)
# dim(wf)
# dim(w)
# all <- rbind(h,wf,w)
# dim(all)
# type <- c(rep("hatchery", times = dim(h)[1])
#           , rep("wild-to-farm", times = dim(wf)[1])
#           , rep("wild", times = dim(w)[1]))
# table(type)
# 
# data <- cbind(all,type)
# head(data)
# tail(data)
# 
# par(mfrow=c(1,1), mar = c(4,4,4,4))
# mod1 <- aov(data$PI ~ data$type)
# 
# boxplot(data$PI ~ data$type, las = 1, ylab = "Pi", ylim = c(0,0.6))
# text(paste("avg=", round(mean(x$Pi.hatchery, na.rm = T), digits = 3)), x = 1, y = 0.45)
# text(paste("avg=", round(mean(x$Pi.wild, na.rm = T), digits = 3)), x = 2, y = 0.45)
# text(paste("avg=", round(mean(x$`Pi.wild-farm`, na.rm = T), digits = 3)), x = 3, y = 0.45)
# summary(mod1)
# TukeyHSD(mod1)
# text(x = 1, y = 0.55, labels = "a")
# text(x = 2, y = 0.55, labels = "b")
# text(x = 3, y = 0.55, labels = "b")
# #text(paste("pval=", round(summary(mod1)[[1]][["Pr(>F)"]][1], digits = 3)), x = 2, y = 0.42)


#### 2. Heterozygosity ####
het <- read.table("batch_1.het", header = T, row.names = 1)
head(het)
het$sample <- rownames(het)

# Add population
sub(pattern = "\\_.*", replacement = "", reads.het$sample)
poptype <- gsub(pattern = "\\_.*", replacement = "", x = het$sample )
het <- cbind(het, poptype)
head(het)

# remove deepbay
het <- het[het$poptype!="DeepBay", ]
het <- droplevels(het)
table(het$poptype)

# Boxplot by population
pdf(file = "Fis_per_pop_boxplot.pdf", width = 8, height = 5)
boxplot(het$F ~ het$poptype, las = 3, ylab = "Fis")
dev.off()

mod1 <- aov(het$F ~ het$poptype)
summary(mod1)
TukeyHSD(mod1)

pdf(file=paste("Fis_per_indiv.pdf",sep=""), height=9, width=14)
par(mar=c(6,5,3,3), mfrow=c(2,1))
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis", las = 1)
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)

#### 3. Read depth per sample
options(scipen = 1)
reads <- read.table("../04-all_samples/reads_per_sample_table.txt", col.names = c("sample","reads"))
head(reads)
reads$sample <- gsub(pattern = ".fq.gz", replacement = "", x = reads$sample)

# plot(reads$reads, xaxt = "n", xlab = "", ylab = "Read depth", las = 1)
# axis(side = 1, labels = reads$sample, at = c(1:nrow(reads)), las = 2, cex.axis = 0.5)
# 
reads.het <- merge(x = reads, y = het, by = "sample")
head(reads.het)
# 
# plot(reads.het$F ~ reads.het$reads, las = 1) # this does not explain it
# # this really shows the two types! Could plot this w/ diff 
# head(reads.het, n = 10)

plot(reads.het$F ~ reads.het$reads, las = 1)
dev.off()


### OLD Qc/BC comparison ###
# # Separate the quebec and bc data
# qc.data <- reads.het[grep(x= reads.het$sample, pattern = "_s", fixed = T), ]
# bc.data <- reads.het[-c(grep(x= reads.het$sample, pattern = "_s", fixed = T)), ]
# 
# plot(reads.het$F ~ reads.het$reads, las = 1, type = "n", ylab = "Fis", xlab = "Number reads per sample")
# points(qc.data$reads, qc.data$F, col = 2)
# points(bc.data$reads, bc.data$F, col = 1)
