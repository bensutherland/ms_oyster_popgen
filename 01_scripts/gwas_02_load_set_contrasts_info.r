# Plot differentiation statistics across genome
## REQUIRES: already run gwas_01_load_constant_info.r

#### 2. Read in contrast files  ####
## Identify file names for contrasts
# Identify sel_perloc_stats, as these are all the contrasts to be used
perloc_stats.files <- list.files(path = "11-adegenet_analysis/", pattern = "_sel_perloc_stats.csv")
datatypes <- gsub(pattern = "_sel_perloc_stats.csv", replacement = "", x = perloc_stats.files) # make short form

# Identify the pcadapt qval files (requires checking)
pcadapt_qvals.files <- list.files(path = "13_selection/", pattern = "qvals_")
rosetta <- as.data.frame(pcadapt_qvals.files, stringsAsFactors = F)
rosetta$datatypes <- c("ch", "dpb", "fr", "gur", "bc", "qdc", "ros", "rsc")
rosetta # confirm that these are corresponding

## Load in per locus contrast metrics (qval and FST)
#perloc.list <- list() ; qval.list <- list() ; 
plotting_data <- list()

for(i in 1:length(datatypes)){
  
  # Identify datatype and report
  datatype <- datatypes[i]
  print(paste0("Now working on ", datatype))
  
  # Identify perloc.FN
  perloc_fst.FN <- paste0("11-adegenet_analysis/", datatype, "_sel_perloc_stats.csv")
  
  # Read in perloc FST data
  perloc_fst <- read.csv(file = perloc_fst.FN)
  head(perloc_fst)
  colnames(perloc_fst)[which(colnames(perloc_fst)=="X")] <- "marker" # rename column with mname
  
  # Format the perloc marker names
  perloc_fst$marker <- gsub(pattern = "X", replacement = "", x = perloc_fst$marker) # drop the X from marker name
  perloc_fst$marker <- substr(x = perloc_fst$marker, start = 1, stop = nchar(perloc_fst$marker)-2) # drop underscore and allele
  head(perloc_fst)
  
  # Order by marker names
  perloc_fst <- perloc_fst[with(perloc_fst, order(marker)), ]
  head(perloc_fst)
  
  ### OLD METHOD 
  # # Put the perloc FST into a list
  # perloc.list[[datatype]] <- perloc_fst
  # head(perloc.list[[datatype]])
  ### END OLD METHOD
  
  # Merge perloc stats with the map info bring the stacks.chr, to match with essential alignment info
  perloc_fst_w_stacks_chr <- merge(x = perloc_fst, y = map, by.x = "marker", by.y = "locus_pos")
  head(perloc_fst_w_stacks_chr)
  dim(perloc_fst_w_stacks_chr) # How many records?
  dim(alignments) # How many records in alignments?
  head(perloc_fst_w_stacks_chr)
  head(alignments)

  ## Combine alignment info with perloc FST 
  # Create matching name with elements from perloc FST
  perloc_fst_w_stacks_chr <- separate(data = perloc_fst_w_stacks_chr, col = "marker", into = c("marker", "marker.locus"), sep = "_")
  head(perloc_fst_w_stacks_chr)
  perloc_fst_w_stacks_chr$matching.mname <- paste0("CLocus_", perloc_fst_w_stacks_chr$marker, "-", perloc_fst_w_stacks_chr$stacks.chr)
  head(perloc_fst_w_stacks_chr)
  
  head(alignments)
  # Bring in the chromosome info
  perloc_fst_w_full_chr_pos <- merge(x = perloc_fst_w_stacks_chr, y = alignments, by.x = "matching.mname", by.y = "mname", all.x = T)
  dim(perloc_fst_w_full_chr_pos)
  head(perloc_fst_w_full_chr_pos)
  
  ## General Stats
  # How many markers per chromosome (Should be constant for all contrasts..)
  table(perloc_fst_w_full_chr_pos$target.chr)
  sum(table(perloc_fst_w_full_chr_pos$target.chr))
  table(is.na(perloc_fst_w_full_chr_pos$target.chr))
  # TRUE shows the markers that are not mapped to the genome (i.e. 2223 of )
  sum(table(is.na(perloc_fst_w_full_chr_pos$target.chr)))
  
  
  # Identify qval.FN (pcadapt); note: this uses the rosetta translation of datatypes to pcadapt qval files
  qval.FN <- paste0("13_selection/", rosetta[rosetta$datatypes==datatype, "pcadapt_qvals.files"])
  qval.df <- read.csv(file = qval.FN)
  head(qval.df)
  qval.df$mname <- gsub(pattern = ">", replacement = "", x = qval.df$mname) # Formatting to match other datasets
  head(qval.df)
  # Calculate -log10(qvalue) for plotting
  qval.df$neg.log10.qval <- -log10(qval.df$qval)
  head(qval.df)
  
  ### OLD METHOD
  # # Put the qval into a list
  # qval.list[[datatype]] <- qval.df
  ### END OLD METHOD

  plotting_data.df <- merge(x = perloc_fst_w_full_chr_pos, y = qval.df, by.x = "matching.mname", by.y = "mname", all.x = TRUE)
  head(plotting_data.df)
  
  # Save into a list
  plotting_data[[datatype]] <- plotting_data.df
}

# Summary # your data is here:
names(plotting_data) # qval per datatype
head(plotting_data[[1]])
datatypes # datatypes
