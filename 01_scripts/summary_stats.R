data <- read.table("reads_per_sample_table.txt") 
print("summary"); summary(data$V2); print("mean"); mean(data$V2)
