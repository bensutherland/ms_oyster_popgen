### Directly compare data without any genome information ####
datatypes <- c("bc", "fr", "ch")
perloc_stats_list <- list()

input_dir <- "11-adegenet_analysis/"


# Loop
for(i in 1:length(datatypes)){
  datatype <- datatypes[i]
  print(datatype)
  filename <- paste(input_dir, datatype, "_sel_perloc_stats.csv", sep = "")

  # Selected data
  perloc_stats <- read.csv(file = filename)
  perloc_stats_list[[datatype]]  <- perloc_stats
}

str(perloc_stats_list)

pdf(file = "11-adegenet_analysis/selection_expt_all_by_all.pdf", width = 8, height = 8)
par(mfrow=c(2,2))

# France vs China
plot(perloc_stats_list[["ch"]]["Fst"]$Fst ~ perloc_stats_list[["fr"]]["Fst"]$Fst
     , ylim = c(-0.05, 0.3)
     , xlim = c(-0.05, 0.3)
     , ylab = "China (Fst)", xlab = "France (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)

# BC vs France
plot(perloc_stats_list[["bc"]]["Fst"]$Fst ~ perloc_stats_list[["fr"]]["Fst"]$Fst
     , ylim = c(-0.05, 0.3)
     , xlim = c(-0.05, 0.3)
     , ylab = "BC (Fst)", xlab = "France (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)

# BC vs China
plot(perloc_stats_list[["bc"]]["Fst"]$Fst ~ perloc_stats_list[["ch"]]["Fst"]$Fst
     , ylim = c(-0.05, 0.3)
     , xlim = c(-0.05, 0.3)
     , ylab = "BC (Fst)", xlab = "China (Fst)")
abline(h = 0.1, lty = 3)
abline(v = 0.1, lty = 3)
dev.off()


# # Find out which markers correspond to this:
# perloc_ch.gid <- as.data.frame(perloc_stats_list$ch.gid)
# perloc_fr.gid <- as.data.frame(perloc_stats_list$fr.gid)
# perloc_bc.sel.gid <- as.data.frame(perloc_stats_list$bc.sel.gid)
# 
# FR <- perloc_fr.gid[which(perloc_fr.gid$Fst > 0.1), "X"]
# BC <- perloc_bc.sel.gid[which(perloc_bc.sel.gid$Fst > 0.1), "X"]
# 
# intersect(x = FR, y = BC) # tells you which markers
# 
# 
# 
# 
