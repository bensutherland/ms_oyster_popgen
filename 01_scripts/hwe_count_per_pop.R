# Identify loci consistently out of HWE
# Requires that there are at least two populations found to be violating HWE for each SNP (within locus)
# If there are multiple SNP in locus, all will be considered individually, and if any
## SNP is found to be out of HWE for multiple pops, then the entire locus is removed

#### Front Matter ####
# Clean space
# rm(list=ls())

# User set variables
max_pops_out_hwe <- 2
pval_cutoff <- 0.01

# Set working directory
setwd("/mnt/data/01_moore_oyster_project/stacks_workflow/")

# Read in sumstats (contains HWE val)
populations_sumstats <- read.table(file = "07-filtered_vcfs/populations.sumstats.tsv", header = F, sep = "\t"
          )
head(populations_sumstats)
nrow(populations_sumstats) # Number of pop x SNP
ncol(populations_sumstats)

# Read in header
header <- as.matrix(read.table(file = "07-filtered_vcfs/populations.sumstats_header.txt", sep = "\t", stringsAsFactors = F))
header <- gsub(pattern = " ", replacement = "_", x = header) # remove spaces


# Annotate sumstats
colnames(populations_sumstats) <- header

head(populations_sumstats)

df <- populations_sumstats

#df <- head(x = df, n = 40000)

locus.ids <- unique(df$Locus_ID)
length(locus.ids)

# Count the instances of HWE equilibrium violation per SNP
sec <- NULL; LOI <- NULL; black.list <- NULL

for(i in 1:length(unique(df$Locus_ID))){
  
  LOI <- locus.ids[i]
  
  # print(locus.ids[i])
  sec <- df[df$Locus_ID==LOI, ]
  sec <- sec[sec$`HWE_P-value` < pval_cutoff, ] # identify how many instances of HWE diseq.
  
  # If the same SNP is out of HWE for at least X populations, put this locus into black list
  if(any(table(sec$Col) > max_pops_out_hwe) == TRUE){
    
     # Add locus to black.list
     black.list <- c(black.list, unique(sec$Locus_ID))
    
  }
}

length(black.list)

print(paste0("Of ", length(locus.ids), " RAD loci, ", length(black.list), " had more than "
             , max_pops_out_hwe, " pops out of HWE at a pval cutoff of ", pval_cutoff))

output.FN <- paste0("01-info_files/blacklist_hwe_less_than_", pval_cutoff, "_more_than_", max_pops_out_hwe, "_pops.txt")

write.table(x = black.list, file = output.FN, sep = "\t", quote = F, row.names = F, col.names = F)
