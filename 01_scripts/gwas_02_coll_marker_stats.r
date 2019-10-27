# Bring in all marker statistics together
# inputs (see set variables below)

# Set variables
FST_percentile <- 0.99

dom.FN  <- "13_selection/rda_analysis/RDA_loading_vals_DOMESTICATION.csv"
farm.FN <- "13_selection/rda_analysis/RDA_loading_vals_FARM.csv"

#### 1.0 Read in rda analysis output ####
DOMESTICATION.rda_loadings <- read.csv(dom.FN)
head(DOMESTICATION.rda_loadings)

FARM.rda_loadings <- read.csv(farm.FN)
head(FARM.rda_loadings)

# Put RDA into a subsettable list
rda_results.list <- list()
rda_results.list[["DOMESTICATION"]] <- DOMESTICATION.rda_loadings
rda_results.list[["FARM"]] <- FARM.rda_loadings


## Bring in correspondence file to link the rda to the genome plot (from the populations.snps.vcf)
correspondence.df <- read.csv("13_selection/chrom_pos_id.csv", header = F)
colnames(correspondence.df) <- c("CHROM", "POS", "ID")
head(correspondence.df)

# Create identifier same as that used by genome plot
correspondence.df$plot.id <- paste0("CLocus_", correspondence.df$ID, "-", correspondence.df$CHROM)
head(correspondence.df)

# Create identifier same as that used by rda_loadings
correspondence.df$rda.id <- paste0(correspondence.df$CHROM, "_", correspondence.df$POS)
head(correspondence.df)


## Read in map file (to bring in FST info)
map <- read.table(file = "11-adegenet_analysis/populations.plink.map")
map <- map[,c(1,2)] # Only need first two cols (stacks chr, locus_pos)
colnames(map) <- c("stacks.chr", "locus_pos")
head(map)
map$locus_pos <- as.character(map$locus_pos) # make mname character
map <- map[with(map, order(locus_pos)), ] # Order by locus_pos for data checking later
head(map)

#### 2. Read in contrast files  ####
# Identify the perloc FST files
perloc_stats.files <- list.files(path = "11-adegenet_analysis/", pattern = "_sel_perloc_stats.csv")
datatypes <- gsub(pattern = "_sel_perloc_stats.csv", replacement = "", x = perloc_stats.files) # make short form

# Identify the pcadapt qval files
pcadapt_qvals.files <- list.files(path = "13_selection/", pattern = "outliers_details")

# MANUAL, create the linking file to connect FST and perloc pop names
rosetta <- as.data.frame(pcadapt_qvals.files, stringsAsFactors = F)
rosetta$datatypes <- c("ch", "dpb", "fr", "gur", "bc", "qdc", "ros", "rsc")
rosetta # MANUAL, confirm that these are corresponding

# Also provide domestication/farm identifier for RDA results
rosetta$rda.type <- c("FARM", "DOMESTICATION", "FARM", "DOMESTICATION", "FARM", "DOMESTICATION", "DOMESTICATION", "DOMESTICATION")

## Load in per locus contrast metrics (FST, qval, RDA loading)
plotting_data <- list() ; rda_datatype <- NULL; datatype <- NULL

for(i in 1:length(datatypes)){
  
  # Identify datatype and report
  datatype <- datatypes[i]
  print(paste0("Now working on ", datatype))
  
  # FST filename
  perloc_fst.FN <- paste0("11-adegenet_analysis/", datatype, "_sel_perloc_stats.csv")
  
  # FST read-in
  perloc_fst <- read.csv(file = perloc_fst.FN)
  head(perloc_fst)
  colnames(perloc_fst)[which(colnames(perloc_fst)=="X")] <- "marker" # rename column with mname
  
  # FST format mnames
  perloc_fst$marker <- gsub(pattern = "X", replacement = "", x = perloc_fst$marker) # drop the X from marker name
  perloc_fst$marker <- substr(x = perloc_fst$marker, start = 1, stop = nchar(perloc_fst$marker)-2) # drop underscore and allele
  head(perloc_fst)
  
  # FST order
  perloc_fst <- perloc_fst[with(perloc_fst, order(marker)), ]
  head(perloc_fst)
  
  # FST merge w/ map info to be able to match w/ alignments
  perloc_fst_w_stacks_chr <- merge(x = perloc_fst, y = map, by.x = "marker", by.y = "locus_pos", sort = F)
  head(perloc_fst_w_stacks_chr)
  dim(perloc_fst_w_stacks_chr) # How many records?
  dim(alignments) # How many records in alignments?
  head(perloc_fst_w_stacks_chr)
  head(alignments)

  # FST prepare matching name for alignments
  perloc_fst_w_stacks_chr <- separate(data = perloc_fst_w_stacks_chr, col = "marker", into = c("marker", "marker.locus"), sep = "_")
  head(perloc_fst_w_stacks_chr)
  perloc_fst_w_stacks_chr$matching.mname <- paste0("CLocus_", perloc_fst_w_stacks_chr$marker, "-", perloc_fst_w_stacks_chr$stacks.chr)
  head(perloc_fst_w_stacks_chr)
  
  
  # CHROM read-in
  head(alignments)
  
  # CHROM merge w/ FST
  perloc_fst_w_full_chr_pos <- merge(x = perloc_fst_w_stacks_chr, y = alignments, by.x = "matching.mname", by.y = "mname", all.x = T, sort = F)
  dim(perloc_fst_w_full_chr_pos)
  head(perloc_fst_w_full_chr_pos)
  
  ## CALC STATS on MARKER per CHROM
  # How many markers per chromosome (Should be constant for all contrasts..)
  table(perloc_fst_w_full_chr_pos$chr)
  sum(table(perloc_fst_w_full_chr_pos$chr))
  table(is.na(perloc_fst_w_full_chr_pos$chr)) # TRUE shows the markers that are not mapped to the genome (i.e. 2223 of )
  sum(table(is.na(perloc_fst_w_full_chr_pos$chr)))
  
  
  # PCADAPT read-in
  # Identify qval.FN (pcadapt); note: this uses the rosetta translation of datatypes to pcadapt qval files
  qval.FN <- paste0("13_selection/", rosetta[rosetta$datatypes==datatype, "pcadapt_qvals.files"])
  qval.df <- read.csv(file = qval.FN)
  head(qval.df)
  qval.df$mname <- gsub(pattern = ">", replacement = "", x = qval.df$mname) # Formatting to match other datasets
  head(qval.df)
  
  # PCADAPT formatting
  qval.df$neg.log10.qval <- -log10(qval.df$qval) # Calculate -log10(qvalue) for plotting
  head(qval.df)
  
  # PCADAPT merge w/ CHROM & FST
  plotting_data.df <- merge(x = perloc_fst_w_full_chr_pos, y = qval.df, by.x = "matching.mname", by.y = "mname", all.x = TRUE, sort = F)
  head(plotting_data.df)
  
  # RDA read-in
  rda_datatype <- rosetta[which(rosetta$datatypes==datatype), "rda.type"]
  rda_results.df <- rda_results.list[[rda_datatype]]
  rda_results.df$rda.id <- as.character(rda_results.df$rda.id)
  str(rda_results.df)
  
  # RDA merge w/ PCADAPT & CHROM & FST, through correspondence.df
  head(plotting_data.df) 
  head(correspondence.df) # file to connect plot data and rda
  plotting_data.df <- merge(x = plotting_data.df, y = correspondence.df, by.x = "matching.mname", by.y = "plot.id", all.x = T, sort = F)
  head(plotting_data.df)
  head(rda_results.df)
  plotting_data.df <- merge(x = plotting_data.df, y = rda_results.df, by = "rda.id", all.x = T, sort = F)
  
  # REPORT the FST 99th percentile used
  # 95th percentile method
  FST_percentile_line <- quantile(x = plotting_data.df$Fst, probs = FST_percentile, na.rm = T)
  plotting_data.df$fst_percentile <- rep(x = FST_percentile_line, times = nrow(plotting_data.df))
  plotting_data.df$fst_percentile_cutoff <- rep(x = FST_percentile, times = nrow(plotting_data.df))
  head(plotting_data.df)

  # SAVE all into a list for plotting
  plotting_data[[datatype]] <- plotting_data.df
  
  # SAVE all into text file
  output.FN <- paste0("13_selection/", datatype, "_all_data_output.txt")
  write.table(x = plotting_data.df, file = output.FN, sep = "\t", row.names = F, quote = F, col.names = T)
  
}

# Summary # your data is here:
names(plotting_data) # qval per datatype
head(plotting_data[[1]])
datatypes # datatypes

save.image("13_selection/constant_and_contrasts_loaded.Rdata")
