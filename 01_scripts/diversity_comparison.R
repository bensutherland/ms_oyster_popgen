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
mean(h$PI)
summary(wf$PI)
mean(wf$PI)
summary(w$PI)
mean(w$PI)

# Merge
x <- merge(h, wf, by = c("CHROM","POS"))
x <- merge(x, w, by = c("CHROM","POS"))
colnames(x)[3] <- "PI.h"
colnames(x)[4] <- "PI.wf"
colnames(x)[5] <- "PI.w"

head(x)
par(mfrow=c(3,1))
plot(x$PI.h, ylab = "Pi (hatchery)", ylim = c(0,1))
plot(x$PI.wf, ylab = "Pi (wild-to-farm)", ylim = c(0,1))
plot(x$PI.w, ylab = "Pi (wild)", ylim = c(0,1))

# How correlated are these values?
cor.test(h$PI, w$PI) # highly correlated
cor.test(wf$PI, w$PI) # very highly correlated

# Model
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

par(mfrow=c(1,1))
mod1 <- aov(data$PI ~ data$type)

boxplot(data$PI ~ data$type, las = 1, ylab = "Pi", ylim = c(0,0.6))
text(paste("avg=", round(mean(x$PI.h), digits = 3)), x = 1, y = 0.45)
text(paste("avg=", round(mean(x$PI.w), digits = 3)), x = 2, y = 0.45)
text(paste("avg=", round(mean(x$PI.wf), digits = 3)), x = 3, y = 0.45)
summary(mod1)
TukeyHSD(mod1)
text(x = 1, y = 0.55, labels = "a")
text(x = 2, y = 0.55, labels = "b")
text(x = 3, y = 0.55, labels = "ab")
#text(paste("pval=", round(summary(mod1)[[1]][["Pr(>F)"]][1], digits = 3)), x = 2, y = 0.42)


#### 2. Heterozygosity ####
het <- read.table("batch_1.het", header = T, row.names = 1)
head(het)
par(mar=c(6,5,3,3), mfrow=c(1,1))
plot(het$F, xaxt ="n", xlab = "", ylab = "Fis")
axis(side = 1, labels = rownames(het), at = c(1:nrow(het)), las = 2, cex.axis = 0.5)


