# Explore genetic diversity differences among hatchery and wild stock
# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_mopup_chips/09-diversity_stats/")

#install.packages("rlist")
require("rlist")

#### 1. Nucleotide Diversity ####
# Set datatypes
datatypes <- c("bc_wild", "bc_wild_to_farm"
               , "chinaQDC", "chinaRSC", "china_wild"
               , "deepbay", "france_wild_to_farm", "france_wild", "guernsey", "japan", "rosewall")

# Import each of the datatypes
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

# Set colors for datatypes
datatypes
cols <- c("green", "darkgreen"
          , "coral3", "coral4", "coral"
          , "royalblue"
          , "purple4", "purple"
          , "paleturquoise"
          , "lightsteelblue"
          , "palevioletred")

# Collect from the datatype list
data <- NULL

# Start the data set by extracting the Pi (3rd address) value for the first population (1)
data <-  as.data.frame(list.extract(datatypes.list, 1))[3]
summary(data$PI)

# Extract the rest of the elements from list
for(i in 2:length(datatypes.list)){
  # Print summaries
  print(datatypes[i])
  print(summary(datatypes.list[[i]]$PI))
  
  # Collect in that datatype into the data dataframe
  data <- cbind(data, as.data.frame(list.extract(datatypes.list, i))[3])
}

# Rename the dataframe columns with the datatypes value
str(data)
colnames(data) <- datatypes
colnames(data)

# Plot
par(mfrow=c(3,4))
for(i in 2:length(datatypes)){
  plot(data$bc_wild ~ data[,i], xlab = colnames(data)[i])
}

# Boxplot Pi per locus averaged across the population
par(mfrow=c(2,1), mar = c(6,4,3,3))
boxplot(data, las = 3, ylab = "Per locus Pi (population mean)"
        , col = cols
        #, xaxt = "n"
        )


#### 2. Heterozygosity ####
het <- read.table("batch_1.het", header = T, row.names = 1)
head(het)
het$sample <- rownames(het)

# Add population to het object
poptype <- NULL
poptype <- sub(pattern = "\\_.*", replacement = "", het$sample)
poptype <- gsub(pattern = "Japan.*", replacement = "Japan", x = poptype ) # group all Japan into one
poptype <- gsub(pattern = "Japan.*", replacement = "Japan", x = poptype ) # group ChinaQDW and ChinaBe into one
poptype <- gsub(pattern = "ChinaQDW|ChinaBe", replacement = "ChinaWild", perl = T, x = poptype)
poptype <- gsub(pattern = "Hisnit|Pendrell|Pipestem|Serpentine", replacement = "BCWild", perl = T, x = poptype)
het <- cbind(het, poptype)
head(het)

# How many per
table(het$poptype)

# Boxplot by population
#pdf(file = "Fis_per_pop_boxplot.pdf", width = 8, height = 5)
boxplot(het$F ~ het$poptype, las = 3, ylab = "Homozygous F statistic (VCFtools)"
        , col = cols)
#dev.off()

mod1 <- aov(het$F ~ het$poptype)
summary(mod1)
result <- TukeyHSD(mod1)

write.csv(x = result$`het$poptype`, file = "tukeyHSD_poptype_het.csv")

pdf(file=paste("Fis_per_indiv.pdf",sep=""), height=9, width=14)
par(mar=c(6,5,3,3), mfrow=c(2,1))
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis", las = 1)
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)


#### 3. Read depth per sample
options(scipen = 1)
reads <- read.table("../04-all_samples/reads_per_sample_table.txt", col.names = c("sample","reads"))
head(reads)
reads$sample <- gsub(pattern = ".fq.gz", replacement = "", x = reads$sample)

reads.het <- merge(x = reads, y = het, by = "sample")
head(reads.het)

pdf(file="het_by_reads.pdf", height=5, width=7)
par(mfrow=c(1,1))
plot(reads.het$F ~ reads.het$reads, las = 1, ylab = "F", xlab = "Reads per Sample")
dev.off()


### OLD CODE: Qc/BC comparison ###
# # Separate the quebec and bc data
# qc.data <- reads.het[grep(x= reads.het$sample, pattern = "_s", fixed = T), ]
# bc.data <- reads.het[-c(grep(x= reads.het$sample, pattern = "_s", fixed = T)), ]
# 
# plot(reads.het$F ~ reads.het$reads, las = 1, type = "n", ylab = "Fis", xlab = "Number reads per sample")
# points(qc.data$reads, qc.data$F, col = 2)
# points(bc.data$reads, bc.data$F, col = 1)
