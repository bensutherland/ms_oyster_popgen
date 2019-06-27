# Plot Fst across genome (oyster)

### 1. read in constant files ####
## Read in alignments
alignments <- read.csv(file = "12_genome_plot/aligned_marker_info.csv", header = F)
colnames(alignments) <- c("mname", "target.chr", "target.loc", "mapq")
head(alignments)

# Keep only top mapq
head(alignments)
length(unique(alignments$mname))
length(alignments$mname) #only a handful of multiple alignments
# TODO# Keep top mapping locus (using mapq, keeping only one row per alignment query)

# Read in map
map <- read.table(file = "11-adegenet_analysis/populations.plink.map")
colnames(map) <- c("mname.chr", "mname.plink", "some.col", "some.col2")
head(map)


#### 2. bring in expt-specific files ####
# Read in necessary, multiple-type data files
datatypes <- c("bc", "ch", "dpb", "fr", "gur", "qdc", "ros", "rsc")
datatype <- "bc"

in_sel_info.FN <- paste0("11-adegenet_analysis/", datatype, "_sel_perloc_stats.csv")
out.FN <- NULL

sel_info <- read.csv(file = in_sel_info.FN)
head(sel_info)
colnames(sel_info)[which(colnames(sel_info)=="X")] <- "marker" # rename columns
head(sel_info)

#### 3. combine data ####
# Combine sel_info and map by "marker" (sel_info) and "mname.chr" and "mname.plink" (map)
# match up the marker name
sel_info$marker <- gsub(pattern = "X", replacement = "", x = sel_info$marker)
head(sel_info)
# remove the last two characters
sel_info$marker <- substr(x = sel_info$marker, start = 1, stop = nchar(sel_info$marker)-2)
head(sel_info)

str(map)
map$mname.plink <- as.character(map$mname.plink)
str(sel_info)
str(map)

# Put into the same order to confirm IDs
map <- map[with(map, order(mname.plink)), ]
head(map)

sel_info <- sel_info[with(sel_info, order(marker)), ]
head(sel_info)

# Merge
data <- merge(x = sel_info, y = map, by.x = "marker", by.y = "mname.plink")
head(data)
dim(data)

# Bring in the alignments, but keep all the data, even those without alignments
dim(data)
dim(alignments)
head(data)
head(alignments)

# merge using mname (in alignments) and marker w/ mname.chr (in data)
require("tidyr")
data <- separate(data = data, col = "marker", into = c("marker", "marker.locus"), sep = "_")
data$matching.mname <- paste0("CLocus_", data$marker, "-", data$mname.chr)

head(data)
head(alignments)
all_data <- merge(x = data, y = alignments, by.x = "matching.mname", by.y = "mname", all.x = T)
dim(all_data)
head(all_data)

# How many markers per chromosome
table(all_data$target.chr)

sum(table(all_data$target.chr))
table(is.na(all_data$target.chr))
# TRUE shows the markers that are not mapped to the genome (i.e. 2223 of )
sum(table(is.na(all_data$target.chr)))

#### 2. Find outliers via quantile or 2xSD of mean ####
head(all_data)

plot(all_data$Fst)
avg_val <- mean(all_data$Fst, na.rm = T)
sd_val  <- sd(all_data$Fst, na.rm = T)
outlier_thresh_2xSD <- avg_val + (2 * sd_val)
outlier_thresh_2xSD

plot(all_data$Fst)
abline(h = outlier_thresh_2xSD)

outlier_thresh_95th_percentile <- quantile(x = all_data$Fst, probs = 0.95, na.rm = T)

outlier_thresh_2xSD
outlier_thresh_95th_percentile

# how many markers are outliers? 
table(all_data$Fst > outlier_thresh_2xSD) # TRUE = outlier

all_data$outlier <- all_data$Fst > outlier_thresh_2xSD
head(all_data)

# Looks like there are some markers that do not have Fst vals
table(is.na(all_data$Fst))
table(is.na(sel_info$Fst)) # these did not come in with Fst vals, they are probably monomorphic


# save out the results (if want) but all is present in all_data

#data is in all_data


##### PLOT ON CHROMOSOMES #####






