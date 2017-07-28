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

# The following will save the image to the directory which the script is launched from
pdf("reads_and_scaff_current.pdf", width = 16, height = 10)

par(mfrow=c(1,1), mar= c(4,4.5,0.5,1) + 0.2, mgp = c(3,0.75,0))
barplot(collective$reads/100, las = 1, names.arg = collective$sample
        , ylab = "Num reads/100 & uniq scaffolds", xaxt = "n", col = "white")
x<- barplot(collective$scaff, las = 1, add = TRUE, col = "grey", xaxt = "n")

labs <- paste(collective$sample)
text(cex = 0.7, 
     x=x-0.25, y=1.25, labs, xpd=T, 
     srt=90
     , pos=1 # 1 is below
     , offset=2.5 # gives distance to offset from pos
)
legend(x = "topright", legend = c("reads","scaff"), fill = c("white","grey"))

# The following turns off the saving out of the image
dev.off()
# save out as 16 x 10 in Portrait


# Produce general summary statistic results to screen
print("summary"); summary(scaff.data$scaff); print("mean"); mean(scaff.data$scaff)

# # Find average percent mapping
# print("percent mapping rate")
# print(mean(scaff.data$scaff/reads.data$reads * 100))


