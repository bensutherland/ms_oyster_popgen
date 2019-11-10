# Explore genetic diversity differences among selected groups
# Clear space
# rm(list=ls())

#### 0. Set up #####
# Install packages
# install.packages("rlist") # note: will require XML (may need to install XML via terminal)
require("rlist")

# Set working directory for stark or xavier
# if on a different system, prompt for working directory
if(Sys.info()["nodename"] == "stark"){ 
  print("On Stark, ready to go")
  setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/09-diversity_stats/") # stark
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, ready to go")
  setwd("/hdd/01_moore_oyster_project/stacks_workflow/09-diversity_stats/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}


#### 1. Nucleotide Diversity ####

### 1.1 Set up nucleotide diversity data ####
## Set datatypes (and order)
datatypes <- c("bc_wild", "bc_wild_to_farm"
               , "deepbay", "rosewall"
               , "france_wild", "france_wild_to_farm"
               , "guernsey"
               , "chinaQDC", "chinaRSC", "china_wild", "china_wild_to_farm"
               )

datatypes.short <- c("BCW", "BCWF"
                      , "DPB", "ROS" 
                      , "FRA", "FRAF" 
                      , "GUR"
                      , "QDC", "RSC", "CHN", "CHNF"
                      )


## Import each of the datatypes
datatypes.list <- list(); datatype <- NULL
for(i in 1:length(datatypes)){
  
  #identify the datatype name
  datatype <- datatypes[i]
  
  # identify the filename
  filename <- paste(datatype,"_samples_out.sites.pi", sep = "")
  
  # read it in
  datatypes.list[[datatype]] <- read.table(filename, header = T)
  
}

names(datatypes.list) # these are the data types
head(datatypes.list[["bc_wild"]]) # this is what the data looks like

## Set colors for datatypes
my.cols <- read.csv(file = "../../ms_oyster_popgen/00_archive/my_cols.csv")
rownames(my.cols) <- my.cols$my.pops

datatypes.short
cols <- as.character(my.cols[c("PEN", "PENF", "DPB", "ROS", "FRA", "FRAF", "GUR", "QDC", "RSC", "CHN", "CHNF"), 2]) # select colors into the same order that they are plotted


## Extract data from list into dataframe
data <- NULL

# Start the data set by extracting the Pi (3rd address) value for the first population (1)
# This is required because cannot fill into an empty list
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
head(data)
colnames(data) <- datatypes.short
head(data)


#### 1.2 Plot nucleotide diversity ####
### 1.2.1 Plot pi per locus all pops vs BC wild (PEN) ####

pdf(file = "pi_per_locus_all_vs_BCW.pdf", width = 15, height = 7)
par(mfrow=c(3,4), mar=c(5,5,3,3))

# Set nulls
adj.rsq <- NULL; mod <- NULL

for(i in 2:length(datatypes)){
  
  # Get the rsq val per comparison
  mod <- lm(data$BCW ~ data[, i])
  adj.rsq <- summary(mod)$adj.r.squared
  adj.rsq <- round(x = adj.rsq, digits = 2)
  
  plot(data$BCW ~ data[,i]
       , xlab = colnames(data)[i]
       , ylab = "bc_wild"
       , ylim = c(0,0.55)
       )
  
  text(x = 0.05, y = 0.52
      , labels = paste0("adj.R2 = ", adj.rsq))
  # legend(x = 0.03, y = 0.45, legend = paste0("adj.R2 = ", adj.rsq), cex = 0.8)
}
dev.off()

### 1.2.2 Plot pi per population each pop individually ####

### WORKING ####
# Need to put data into long form
colnames(data)
NROW(data)
long.form.grouping <- rep(x = colnames(data), each = (NROW(data)))
head(long.form.grouping)
tail(long.form.grouping)

long.form.vals <- as.vector(as.matrix(data))
long.form.data <- cbind(long.form.grouping, long.form.vals)
long.form.data.df <- as.data.frame(long.form.data)
head(long.form.data.df)
colnames(long.form.data.df) <- c("pop", "nuc.div")
long.form.data.df$nuc.div <- as.numeric(as.character(long.form.data.df$nuc.div))
head(long.form.data.df)
head(data) # confirm its ok

# Test plotting, note that the order will be different than below's method
boxplot(long.form.data.df$nuc.div ~ long.form.data.df$pop) # looks good though
mod1 <- aov(long.form.data.df$nuc.div ~ long.form.data.df$pop)
summary(mod1)
TukeyHSD(x = mod1)

# Working to find a solution for non-normal data
kruskal.test(long.form.data.df$nuc.div ~ long.form.data.df$pop)

# install.packages("pgirmess")
# library("pgirmess")
# kruskalmc(long.form.data.df$nuc.div ~ long.form.data.df$pop)


#### PLOTTING ####

pdf(file = "pi_per_pop_mean.pdf", width = 8, height = 5)
par(mfrow=c(1,1), mar = c(5,5,3,3))
boxplot(data, las = 1, ylab = "Nucleotide Diversity (pi)"
        , col = cols
        , xaxt = "n"
        , ylim = c(0,0.6)
        )
axis(side = 1, at = c(1:length(datatypes.short)), labels = datatypes.short, las = 2)
dev.off()


##### Plot distributions #### WORKING
percent.zero <- NULL; nuc_div.FN <- NULL

colnames(data)

pdf(file = "nuc_div_hist.pdf", width =7, height = 9)
par(mfrow = c(4,3), mar = c(4,2,2,2))
for(i in 1:length(colnames(data))){
  
  # Find your percent zero
  percent.zero <- (table(data[,i]==0)["TRUE"] / sum(table(data[,i]))) * 100
  percent.zero <- round(percent.zero, digits = 1)
  percent.zero <- as.numeric(percent.zero)
  
  # Plot
  hist(data[,i], breaks = 100, ylim = c(0, 10000), main = ""
       , xlab = paste0(colnames(data)[i], " nucleotide diversity")
       )
  text(x = 0.4, y = 1000, labels = paste0("% 0 = ", percent.zero))
  text(x = 0.4, y = 3000, labels = paste0("avg "
                                          , round(mean(x = data[,i]), digits = 4)
                                          )
       )
  text(x = 0.4, y = 5000, labels = paste0("med "
                                          , round(median(x = data[,i]), digits = 4)
  )
  )
  
  
  
}
dev.off()

# Testing
par(mfrow = c(2,2), mar = c(3,3,3,3))

boxplot(data$FRA)
hist(data$FRA, breaks = 50)

hist(data$DPB, breaks = 50)
boxplot(data$DPB)



#DPB
percent.zero <- (table(data$DPB==0)["TRUE"] / sum(table(data$DPB))) * 100
percent.zero <- round(percent.zero, digits = 1)
percent.zero<- as.numeric(percent.zero)

hist(data$DPB, breaks = 100, ylim = c(0,10000))
text(x = 0.4, y = 1000, labels = paste0("% 0 = ", percent.zero))



### 1.3 Nucleotide diversity stats ####
summary(data)


# USING LONG FORM DF (NOT READY)
# # Note: Not yet implemented as not sure on the best approach. May just use descriptive statistics
# # head(data)
# 
# # Convert data into long form
# data.long <- NULL; data.long.pop <- NULL; data.long.df <- NULL
# 
# for(i in 1:length(names(data))){
#   # Extract the values for the population, and make a second column with the pop name
#   data.long <- c(data.long, data[,i])
#   data.long.pop <- c(data.long.pop, rep(colnames(data)[i], times = length(data[,i])))
#   data.long.df <- as.data.frame(cbind(data.long.pop, data.long)) # combine
# 
#   # Name columns
#   colnames(data.long.df) <- c("pop", "per.locus.pi")
# 
#   # Change the nucleotide diversity to numeric
#   data.long.df$per.locus.pi <- as.numeric(as.character(data.long.df$per.locus.pi))
#   #data.long.df$pop <- as.character(data.long.df$pop)
# }
# 
# head(data.long.df)
# tail(data.long.df)
# str(data.long.df)
# 
# # Confirm the data still looks as expected via boxplot (need to set an order of levels)
# pop.o <- ordered(data.long.df$pop, levels = c("BCW", "BCWF", "DPB", "ROS", "FRA", "FRAF", "GUR", "QDC", "RSC", "CHN", "CHNF"))
# str(pop.o)
# boxplot(per.locus.pi ~ pop.o, data = data.long.df, ylim = c(0,0.6), col = cols, las = 2)
# 
# ## The following also works, but it will not be in the correct order for comparison
# # boxplot(data.long.df$per.locus.pi ~ data.long.df$pop, ylim = c(0,0.6), las = 2
# #         )
# 
# # Run model
# mod.pi <- aov(data.long.df$per.locus.pi ~ data.long.df$pop)
# summary(mod.pi) # significant overall
# TukeyHSD(mod.pi) # no pairwise significant differences
# 
# # What are the unique population names?
# unique(data.long.df$pop)
# 
# # Other methods (rank sum)
# mod.kr.pi <- kruskal.test(data.long.df$per.locus.pi ~ data.long.df$pop)
# mod.kr.pi # highly significant
# 
# # Pairwise comparisons using Wilcoxon rank sum test
# pairwise.wilcox.test(data.long.df$per.locus.pi, data.long.df$pop, p.adjust.method = "BH")
# 
# # Just compare between deepbay and bc wild
# data.long.bcw.dpb <- data.long.df[data.long.df$pop=="BCW"|data.long.df$pop=="DPB",]
# wilcox.test(data.long.bcw.dpb$per.locus.pi ~ data.long.bcw.dpb$pop) # highly signfiicant
# summary(aov(data.long.bcw.dpb$per.locus.pi ~ data.long.bcw.dpb$pop)) # not significant at all (p = 0.913); why such large differences?


#### 2. Heterozygosity ####
### 2.1 Set up heterozygosity data ####
het <- read.table("populations.snps.het", header = T, row.names = 1)
het$sample <- rownames(het)
head(het)

# Add population to het object
het$poptype <- sub(pattern = "\\_.*", replacement = "", het$sample)
het$poptype <- gsub(x = het$poptype, pattern = "HIS|PEN|PIP|SER", replacement = "BCW", perl = T) # merge all BC wild pops
het <- het[grep(pattern = "JPN", x = het$poptype, invert = T), ] # drop japan due to low sample size
het$poptype <- as.factor(het$poptype) # make factor for plotting

head(het)
table(het$poptype)

### 2.2 Plot heterozygosity data ####
# Conduct anova and pairwise tukeys to investigate differences among poptypes
mod.fis <- aov(het$F ~ het$poptype)
summary(mod.fis)
result <- TukeyHSD(mod.fis)
het_result <- result[[1]]

write.csv(x = het_result, file = "tukeyHSD_poptype_het.csv")

# Based on these results, make significance value indicators for the plot (if doesn't share letter, sig diff p< 0.05)
sig.indicators <- c("ac", "abc", "b", "ab", "abc", "ab", "ab", "abc", "ac", "c", "ac")
length(sig.indicators)

# Boxplot by population
pdf(file = "Fis_per_pop_boxplot.pdf", width = 8, height = 5)
pop.o.het <- ordered(het$poptype, levels = c("BCW", "BCWF", "DPB", "ROS", "FRA", "FRAF", "GUR", "QDC", "RSC", "CHN", "CHNF")) # set order
boxplot(F ~ pop.o.het, data = het 
        , col = cols, las = 2
        , ylim = c(-0.1, 0.3)
        , ylab = "Fis")
text(x = seq(1:length(levels(het$poptype))), y = 0.28
     , labels = sig.indicators
      )
dev.off()

#### 2.3 Summarize het data ###
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

# First plot het by sample index
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis", las = 1)
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)

# Second plot het by reads
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


#### 3. Plot both in same image
# this is just copied from above for combining

pdf(file = "genetic_diversity_all.pdf", width = 6, height = 8)
par(mfrow=c(2,1), mar = c(4,4,3,3))

# Plot 1 (pi)
boxplot(data, las = 1, ylab = "Nucleotide Diversity (pi)"
        , col = cols
        , xaxt = "n"
        , ylim = c(0,0.6)
)
axis(side = 1, at = c(1:length(datatypes.short)), labels = datatypes.short, las = 2)

# Plot 2 (het)
pop.o.het <- ordered(het$poptype, levels = c("BCW", "BCWF", "DPB", "ROS", "FRA", "FRAF", "GUR", "QDC", "RSC", "CHN", "CHNF")) # set order
boxplot(F ~ pop.o.het, data = het 
        , col = cols, las = 2
        , ylim = c(-0.1, 0.3)
        , ylab = "Inbreeding Coeff. (Fis)"
        )
text(x = seq(1:length(levels(het$poptype))), y = 0.28
     , labels = sig.indicators
)
dev.off()



