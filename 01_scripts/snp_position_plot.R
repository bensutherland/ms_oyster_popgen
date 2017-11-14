setwd("~/Documents/01_moore_oyster_project/stacks_workflow_2017-11-09")

snp.pos <- read.table(file = "snp_positions_only_positions.txt")
snp.pos.num <- as.numeric(snp.pos$V1)
hist(snp.pos.num, main = "", las = 1, breaks=100
     , xlab = "Position in Tag (bp)"
     , ylab = "Number SNPs")

text(x = (max(snp.pos.num)-0.1*max(snp.pos.num))
     , y = 1700, labels = paste("total =", length(snp.pos.num)))
