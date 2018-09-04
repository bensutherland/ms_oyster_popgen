# Explore genetic diversity differences among hatchery and wild stock
# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/09-diversity_stats/")

#install.packages("rlist")
require("rlist")

#### 1. Nucleotide Diversity ####
# Set datatypes
datatypes <- c("bc_wild", "bc_wild_to_farm"
               , "chinaQDC", "chinaRSC", "china_wild", "china_wild_to_farm"
               , "deepbay", "france_wild", "france_wild_to_farm", "guernsey", "japan", "rosewall")

datatypes.short <- c("BCW", "BCWF"
               , "QDC", "RSC", "CHN", "CHNF"
               , "DPB", "FRA", "FRAF", "GUR", "JPN", "ROS")


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
my.cols <- read.csv("../01-info_files/my_cols.csv")

datatypes
cols <- as.character(my.cols[c(9,10,12,14,1,2,3,4,5,6,8,13), 2]) # select colors into the same order that they are plotted
# in names(data)

# cols <- c("mediumorchid1"
#           , "purple1"
#           , "darkseagreen1", "darkseagreen"
#           , "darkseagreen4"
#           , "darkseagreen3"
#           , "dodgerblue2"
#           , "gold1", "gold3"
#           , "darkslategray2"
#           , "turquoise4"
#           , "red")

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
pdf(file = "Pi_per_locus_all_comp_to_bc_wild.pdf", width = 15, height = 7)
par(mfrow=c(3,4), mar=c(5,5,3,3))
for(i in 2:length(datatypes)){
  plot(data$bc_wild ~ data[,i], xlab = colnames(data)[i])
}
dev.off()

# Boxplot Pi per locus averaged across the population
pdf(file = "Pi_per_pop_mean.pdf", width = 8, height = 5)
par(mfrow=c(1,1), mar = c(8,4,3,3))
boxplot(data, las = 3, ylab = "Per locus Pi (population mean)"
        , col = cols
        , xaxt = "n"
        , ylim = c(0,0.6)
        )
axis(side = 1, at = c(1:length(datatypes.short)), labels = datatypes.short, las = 2)
dev.off()

# put data into a longform dataframe that can be used for stats
str(data)
data.long <- NULL
data.long.pop <- NULL

for(i in 1:length(names(data))){
  data.long <- c(data.long, data[,i])
  data.long.pop <- c(data.long.pop, rep(colnames(data)[i], times = length(data[,i])))
}

length(data.long)
length(data.long.pop)

data.long.df <- as.data.frame(cbind(data.long.pop, data.long))
dim(data.long.df)
colnames(data.long.df) <- c("pop", "per.locus.pi")
data.long.df$per.locus.pi <- as.numeric(as.character(data.long.df$per.locus.pi))
str(data.long.df)

head(data.long.df)

#

# test w/ boxplot
boxplot(data.long.df$per.locus.pi ~ data.long.df$pop, ylim = c(0,0.6), col = cols, las = 2)

mod.pi <- aov(data.long.df$per.locus.pi ~ data.long.df$pop) 
summary(mod.pi)
TukeyHSD(mod.pi)

levels(data.long.df$pop)

mod.kr.pi <- kruskal.test(data.long.df$per.locus.pi ~ data.long.df$pop)

pairwise.wilcox.test(data.long.df$per.locus.pi, data.long.df$pop, p.adjust.method = "BH")

# Just compare between deepbay and bc wild
data.long.bcw.dpb <- data.long.df[data.long.df$pop=="bc_wild"|data.long.df$pop=="deepbay",]
wilcox.test(data.long.bcw.dpb$per.locus.pi ~ data.long.bcw.dpb$pop)
summary(aov(data.long.bcw.dpb$per.locus.pi ~ data.long.bcw.dpb$pop))



#### 2. Heterozygosity ####
het <- read.table("batch_1.het", header = T, row.names = 1)
head(het)
het$sample <- rownames(het)
head(het)

# Add population to het object
poptype <- NULL
poptype <- sub(pattern = "\\_.*", replacement = "", het$sample)
# poptype <- gsub(pattern = "HIS|PEN|PIP|SER", replacement = "BCWild", perl = T, x = poptype) # this could be an option..
het <- cbind(het, poptype)
head(het)

# How many per
table(het$poptype)

# Conduct anova and pairwise tukeys to investigate differences among poptypes
mod.fis <- aov(het$F ~ het$poptype)
summary(mod.fis)
result <- TukeyHSD(mod.fis)

write.csv(x = result$`het$poptype`, file = "tukeyHSD_poptype_het.csv")

# Based on these results, make significance value indicators for the plot (if doesn't share letter, sig diff p< 0.05)
sig.indicators <- c("a", "a", "b", "bc", "bc"
                  , "bc", "ab", "ab", "ac", "ac"
                  , "ab", "ab", "bc", "ab", "ac")
length(sig.indicators)

# Boxplot by population
pdf(file = "Fis_per_pop_boxplot.pdf", width = 8, height = 5)
boxplot(het$F ~ het$poptype, las = 1, ylab = "Fis"
        , col = as.character(my.cols$my.cols)
        , ylim = c(-0.1,0.3)
        , xaxt = "n")
axis(side = 1, labels = levels(het$poptype), at = c(1:length(levels(het$poptype))), las = 2)
text(x = seq(1:length(levels(het$poptype))), y = 0.25
     , labels = sig.indicators
     )

dev.off()

# Get some summaries
library("dplyr")
summarized <- het %>%
  select(F, poptype) %>%
  group_by(poptype) %>%
  summarise(F.avg = mean(F), F.med = median(F))

summarized
write.csv(x = summarized, file = "F_avg_and_median.csv", quote = F)

### 3. Plot Fis by read depth per sample
pdf(file=paste("Fis_and_read_depth_per_indiv.pdf",sep=""), height=9, width=14)
par(mar=c(6,5,3,3), mfrow=c(2,1))
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis", las = 1)
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)


#### 3. Read depth per sample
options(scipen = 1)
reads <- read.table("../04-all_samples/reads_per_sample_table.txt", col.names = c("sample","reads"))
head(reads)
reads$sample <- gsub(pattern = ".fq.gz", replacement = "", x = reads$sample)

# combine the reads and het dataframes
reads.het <- merge(x = reads, y = het, by = "sample")
head(reads.het)

#pdf(file="het_by_reads.pdf", height=5, width=7)
#par(mfrow=c(1,1))
plot(reads.het$F ~ reads.het$reads, las = 1, ylab = "F", xlab = "Reads per Sample")
dev.off()

