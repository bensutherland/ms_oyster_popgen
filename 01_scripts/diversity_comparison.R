# Explore genetic diversity differences among hatchery and wild stock
# Clear space
# rm(list=ls())

# Set working directory
setwd("~/Documents/01_moore_oyster_project/stacks_workflow/09-diversity_stats/")

#### 1. Nucleotide Diversity ####

# Import data
h <- read.table("hatchery_out.sites.pi", header=T)
wf <- read.table("wild_to_farm_out.sites.pi", header = T)
w <- read.table("wild_out.sites.pi", header=T)
head(h)
head(wf)
head(w)

# Explore data
summary(h$PI)
summary(wf$PI)
summary(w$PI)

# Merge
x <- merge(h, wf, by = c("CHROM","POS"))
x <- merge(x, w, by = c("CHROM","POS"))
colnames(x)[3] <- "Pi.hatchery"
colnames(x)[4] <- "Pi.wild-farm"
colnames(x)[5] <- "Pi.wild"

head(x)

#### 1.1. Pi per locus across types ####

pdf(file=paste("Pi_per_locus_by_type.pdf",sep=""), height=5, width=7)
par(mfrow=c(1,3))
plot(x$Pi.hatchery, ylab = "Pi (hatchery)", xlab = "Locus.index", ylim = c(0,1), las = 1)
plot(x$`Pi.wild-farm`, ylab = "Pi (wild-to-farm)", xlab = "Locus.index", ylim = c(0,1), las = 1)
plot(x$Pi.wild, ylab = "Pi (wild)", xlab = "Locus.index", ylim = c(0,1), las = 1)
dev.off()

#### 1.2 Pi correlation ####
# How correlated are the pi values across loci for the different types (hatchery, wild and wild-farm)?
#install.packages('corrplot')
library('corrplot')
cor.test(h$PI, w$PI) # correlated (cor = 0.52)
cor.test(h$PI, wf$PI) # correlated (cor = 0.48)
cor.test(wf$PI, w$PI) # highly correlated (cor = 0.84)

cor.set <- cor(x[,c("Pi.hatchery","Pi.wild-farm","Pi.wild")]
               , use = "pairwise.complete.obs" # cor bw ea pair var is comp using all compl pairs of obs on those var # note: only works with "pearson" method    
)


# Plot correlation
pdf(file=paste("corr_of_pi_bw_groups.pdf",sep=""), height=5, width=5)
par(mfrow=c(1,1), mar = c(3,3,4,5))
corrplot(cor.set, method="number", addshade="all", tl.col = "black", tl.cex=0.8, type = "lower", diag = T, tl.pos = "lt", cl.pos = "r")
corrplot(cor.set, method="circle", addshade="all", tl.col = "black", tl.cex=0.8, type = "upper", add = T, tl.pos = "n")
dev.off()

#### 1.3. Pi differences across groups (model) ####
dim(h)
dim(wf)
dim(w)
all <- rbind(h,wf,w)
dim(all)
type <- c(rep("hatchery", times = dim(h)[1])
          , rep("wild-to-farm", times = dim(wf)[1])
          , rep("wild", times = dim(w)[1]))
table(type)

data <- cbind(all,type)
head(data)
tail(data)

par(mfrow=c(1,1), mar = c(4,4,4,4))
mod1 <- aov(data$PI ~ data$type)

boxplot(data$PI ~ data$type, las = 1, ylab = "Pi", ylim = c(0,0.6))
text(paste("avg=", round(mean(x$Pi.hatchery, na.rm = T), digits = 3)), x = 1, y = 0.45)
text(paste("avg=", round(mean(x$Pi.wild, na.rm = T), digits = 3)), x = 2, y = 0.45)
text(paste("avg=", round(mean(x$`Pi.wild-farm`, na.rm = T), digits = 3)), x = 3, y = 0.45)
summary(mod1)
TukeyHSD(mod1)
text(x = 1, y = 0.55, labels = "a")
text(x = 2, y = 0.55, labels = "b")
text(x = 3, y = 0.55, labels = "b")
#text(paste("pval=", round(summary(mod1)[[1]][["Pr(>F)"]][1], digits = 3)), x = 2, y = 0.42)


#### 2. Heterozygosity ####
het <- read.table("batch_1.het", header = T, row.names = 1)
head(het)

pdf(file=paste("Fis_per_indiv.pdf",sep=""), height=5, width=14)
par(mar=c(6,5,3,3), mfrow=c(1,1))
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis")
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)
dev.off()

