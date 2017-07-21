# To perform some basic stats on the reads per sample table

# Read in data
data <- read.table("04-all_samples/reads_per_sample_table.txt", col.names = c("sample", "reads"))
# Clean it up
data$sample <- gsub(pattern = ".fq.gz", replacement = "", x = data$sample)

# Read in other data
data2 <- read.table("04-all_samples/mappings_per_sample_table.txt"
                    , col.names = c("sample", "mappings"))

# Clean it up
data2$sample <- gsub(pattern = ".bam", replacement = "", x = data2$sample)

# Merge
collective <- merge(x = data, y = data2, by = "sample")


# Plot

jpeg("reads_and_mappings_per_sample.jpg")
par(mfrow=c(1,1), mar= c(2,5,0.5,1) + 0.2, mgp = c(2,0.75,0))

par(mfrow=c(1,1), mar= c(2,5,0.5,1) + 0.2, mgp = c(2,0.75,0))
barplot(collective$reads, las = 1, names.arg = collective$sample, ylab = "Num reads or mappings", xaxt = "n", col = "white")

x<- barplot(collective$mappings, las = 1, add = TRUE, col = "grey", xaxt = "n")

labs <- paste(collective$sample)
text(cex = 0.8, x=x-0.25, y=-1.25, labs, xpd=T, srt=45)
legend(x = "topright", legend = c("reads","mappings"), fill = c("white","black"))
dev.off()

# Produce general summary statistic results to screen
print("summary"); summary(data$reads); print("mean"); mean(data$reads)

# Produce general summary statistic results to screen
print("summary"); summary(data2$mappings); print("mean"); mean(data2$mappings)

