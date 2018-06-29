# Plotting script for a barplot showing the number of reads (white)
# by the number of mappings
# B. Sutherland 2017-07-26

reads.data <- read.table("04-all_samples/reads_per_sample_table.txt"
                   , col.names = c("sample", "reads"))
reads.data

reads.data$sample <- gsub(pattern = ".fq.gz", replacement = "", x = reads.data$sample)
reads.data

mappings.data <- read.table("04-all_samples/mappings_per_sample_table.txt"
                    , col.names = c("sample", "mappings"))
mappings.data

mappings.data$sample <- gsub(pattern = ".bam", replacement = "", x = mappings.data$sample)
mappings.data

collective <- merge(x = reads.data, y = mappings.data, by = "sample")

# The following will save the image to the directory which the script is launched from
pdf("reads_and_mappings_current.pdf", width = 16, height = 10)

par(mfrow=c(1,1), mar= c(4,4.5,0.5,1) + 0.2, mgp = c(3,0.75,0))
barplot(collective$reads, las = 1, names.arg = collective$sample, ylab = "Num reads or mappings", xaxt = "n", col = "white")
x<- barplot(collective$mappings, las = 1, add = TRUE, col = "grey", xaxt = "n")

labs <- paste(collective$sample)
text(cex = 0.7, 
     x=x-0.25, y=1.25, labs, xpd=T, 
     srt=90
     , pos=1 # 1 is below
     , offset=2.5 # gives distance to offset from pos
)
legend(x = "topright", legend = c("reads","mappings"), fill = c("white","grey"))
abline(h = 1500000, lty = 2)
abline(h = 1000000, lty = 2)


# The following turns off the saving out of the image
dev.off()
# save out as 16 x 10 in Portrait


# Produce general summary statistic results to screen
print("summary"); summary(reads.data$reads); print("mean"); mean(reads.data$reads)

# Produce general summary statistic results to screen
print("summary"); summary(mappings.data$mappings); print("mean"); mean(mappings.data$mappings)

# Find average percent mapping
print("percent mapping rate")
print(mean(mappings.data$mappings/reads.data$reads * 100))


