# Plot differentiation statistics across genome
## REQUIRES: 
# already run gwas_01_load_constant_info.r

# Set variables
FST_percentile <- 0.99

#### 1.0 Read in rda analysis output ####
DOMESTICATION.rda_loadings <- read.csv("13_selection/rda_analysis/RDA_loading_vals_DOMESTICATION.csv")
head(DOMESTICATION.rda_loadings)

FARM.rda_loadings <- read.csv("13_selection/rda_analysis/RDA_loading_vals_FARM.csv")
head(FARM.rda_loadings)

# put into list
rda_results.list <- NULL
as.list(rda_results.list)
rda_results.list[["DOMESTICATION"]] <- DOMESTICATION.rda_loadings
rda_results.list[["FARM"]] <- FARM.rda_loadings


## Bring in correspondence file to link the rda to the genome plot
correspondence.df <- read.csv("13_selection/chrom_pos_id.csv", header = F)
colnames(correspondence.df) <- c("CHROM", "POS", "ID")
head(correspondence.df)

# Create identifier same as that used by genome plot
correspondence.df$plot.id <- paste0("CLocus_", correspondence.df$ID, "-", correspondence.df$CHROM)
head(correspondence.df)

# Create identifier same as that used by rda_loadings
correspondence.df$rda.id <- paste0(correspondence.df$CHROM, "_", correspondence.df$POS)
head(correspondence.df)


#### 2. Read in contrast files  ####
## Identify file names for contrasts
# Identify sel_perloc_stats, as these are all the contrasts to be used
perloc_stats.files <- list.files(path = "11-adegenet_analysis/", pattern = "_sel_perloc_stats.csv")
datatypes <- gsub(pattern = "_sel_perloc_stats.csv", replacement = "", x = perloc_stats.files) # make short form

# Identify the pcadapt qval files (requires checking)

### NOTE: HARD CODED ###
#pcadapt_qvals.files <- list.files(path = "13_selection/", pattern = "qvals_")
pcadapt_qvals.files <- list.files(path = "13_selection/", pattern = "outliers_details")
rosetta <- as.data.frame(pcadapt_qvals.files, stringsAsFactors = F)
rosetta$datatypes <- c("ch", "dpb", "fr", "gur", "bc", "qdc", "ros", "rsc")
rosetta # confirm that these are corresponding
# Also provide domestication/farm identifier for RDA results
rosetta$rda.type <- c("FARM", "DOMESTICATION", "FARM", "DOMESTICATION", "FARM", "DOMESTICATION", "DOMESTICATION", "DOMESTICATION")

## Load in per locus contrast metrics (FST, qval, RDA loading)
#perloc.list <- list() ; qval.list <- list() ; 
plotting_data <- list() ; rda_datatype <- NULL

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
  
  # Merge perloc stats w/ map info to bring in the stacks.chr to match w/ alignments
  perloc_fst_w_stacks_chr <- merge(x = perloc_fst, y = map, by.x = "marker", by.y = "locus_pos", sort = F)
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
  
  # Bring in the chromosome info
  head(alignments)
  perloc_fst_w_full_chr_pos <- merge(x = perloc_fst_w_stacks_chr, y = alignments, by.x = "matching.mname", by.y = "mname", all.x = T, sort = F)
  dim(perloc_fst_w_full_chr_pos)
  head(perloc_fst_w_full_chr_pos)
  
  ## General Stats
  # How many markers per chromosome (Should be constant for all contrasts..)
  table(perloc_fst_w_full_chr_pos$chr)
  sum(table(perloc_fst_w_full_chr_pos$chr))
  table(is.na(perloc_fst_w_full_chr_pos$chr)) # TRUE shows the markers that are not mapped to the genome (i.e. 2223 of )
  sum(table(is.na(perloc_fst_w_full_chr_pos$chr)))
  
  
  # Identify qval.FN (pcadapt); note: this uses the rosetta translation of datatypes to pcadapt qval files
  qval.FN <- paste0("13_selection/", rosetta[rosetta$datatypes==datatype, "pcadapt_qvals.files"])
  qval.df <- read.csv(file = qval.FN)
  head(qval.df)
  qval.df$mname <- gsub(pattern = ">", replacement = "", x = qval.df$mname) # Formatting to match other datasets
  head(qval.df)
  
  # Calculate -log10(qvalue) for plotting
  qval.df$neg.log10.qval <- -log10(qval.df$qval)
  head(qval.df)
  
  plotting_data.df <- merge(x = perloc_fst_w_full_chr_pos, y = qval.df, by.x = "matching.mname", by.y = "mname", all.x = TRUE, sort = F)
  head(plotting_data.df)
  
  
  # Bring in RDA info
  rda_datatype <- rosetta[which(rosetta$datatypes==datatype), "rda.type"]
  rda_results.df <- rda_results.list[[rda_datatype]]
  rda_results.df$rda.id <- as.character(rda_results.df$rda.id)
  str(rda_results.df)
  
  # Connect to plotting data
  # First bring in a file that has a connector between the plotting data and rda
  head(plotting_data.df)
  head(correspondence.df)
  plotting_data.df <- merge(x = plotting_data.df, y = correspondence.df, by.x = "matching.mname", by.y = "plot.id", all.x = T, sort = F)
  
  # Then bring in the rda
  head(plotting_data.df)
  head(rda_results.df)
  
  plotting_data.df <- merge(x = plotting_data.df, y = rda_results.df, by = "rda.id", all.x = T, sort = F)
  
  # Add in the FST 99th percentile to the output file
  # 95th percentile method
  FST_percentile_line <- quantile(x = plotting_data.df$Fst, probs = FST_percentile, na.rm = T)
  plotting_data.df$fst_percentile <- rep(x = FST_percentile_line, times = nrow(plotting_data.df))
  plotting_data.df$fst_percentile_cutoff <- rep(x = FST_percentile, times = nrow(plotting_data.df))
  head(plotting_data.df)

  # Save into a list to use for plotting
  plotting_data[[datatype]] <- plotting_data.df
  
  # Save into a text file for potential future use
  output.FN <- paste0("13_selection/", datatype, "_all_data_output.txt")
  write.table(x = plotting_data.df, file = output.FN, sep = "\t", row.names = F, quote = F, col.names = T)
  
}

# Summary # your data is here:
names(plotting_data) # qval per datatype
head(plotting_data[[1]])
datatypes # datatypes

save.image("13_selection/constant_and_contrasts_loaded.Rdata")
