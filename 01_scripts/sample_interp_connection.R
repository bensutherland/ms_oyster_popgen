# Script to obtain data about specific samples
interp <- read.csv(file = "~/Documents/01_moore_oyster_project/sample_database_exported_2018-05-30.csv", header = T, sep = ",")
dim(interp)
head(interp)
colnames(interp)

interp$Sample.Number <- as.numeric(interp$Sample.Number)


indiv.used
merge(x = indiv.used, y = interp, by.x = "indiv.used", by.y = "Sample.Number")
