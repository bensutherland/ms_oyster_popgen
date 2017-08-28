# Plotting script for a barplot showing the number of reads (white)
# by the number of unique scaffolds being mapped
# B. Sutherland 2017-07-28

reads.data <- read.table("04-all_samples/reads_per_sample_table.txt"
                   , col.names = c("sample", "reads"))
reads.data

reads.data$sample <- gsub(pattern = ".fq.gz", replacement = "", x = reads.data$sample)
reads.data

scaff.data <- read.table("04-all_samples/scaff_per_sample_table.txt"
                    , col.names = c("sample", "scaff"))
scaff.data

scaff.data$sample <- gsub(pattern = ".bam", replacement = "", x = scaff.data$sample)
scaff.data

collective <- merge(x = reads.data, y = scaff.data, by = "sample")

# Bring in coverage data from pstacks
loci.cov <- read.csv(file = "output_pstacks_results.csv")

# Merge with other data
all.data <- merge(x=collective, y = loci.cov, by = "sample")

# Plot
# The following will save the image to the directory which the script is launched from
pdf("num_loci_per_sample_and_avg_coverage.pdf", width = 7, height = 4)

par(mfrow=c(1,2), mar= c(4,4.5,0.5,1) + 0.2, mgp = c(3,0.75,0))
plot(all.data$reads, all.data$avg.cov
     , ylab = "mean cov per locus", xlab = "num reads per sample"
     , las = 1)
plot(all.data$reads, all.data$num.loci
     , ylab = "total loci per sample", xlab = "num reads per sample"
     , las = 1)


# The following turns off the saving out of the image
dev.off()
