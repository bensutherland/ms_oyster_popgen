#### TOPLOADING STUFF ####
# Identify the top loading marker names
toploading.markers <- dimnames(dapc2$var.contr)[[1]][which(dapc2$var.contr > 0.0005)]
toploading.markers <- gsub(pattern = "\\..*", replacement = "", toploading.markers) # match the genlight name format (remove .1 or .2)
toploading.markers

# # subset the gid file w/ toploading markers
# locNames(wild.bc.gid)
# 
# wild.bc.gid.toploading <- wild.bc.gid[, loc=toploading.markers]
# wild.bc.gid.toploading

# Subset the original my.data genlight file with the toploading markers and convert back to a genlight file
my.data.mat <- as.matrix(my.data)
my.data.mat[1:5,1:5]
# retain only top loading markers
my.data.mat.toploading <- my.data.mat[,toploading.markers]
# retain only wild samples (removing DeepBay only currently)
samples.to.keep <- rownames(my.data.mat)[grep(pattern = "DeepBay", x = rownames(my.data.mat), invert = T)]
my.data.mat.toploading.wild <- my.data.mat.toploading[samples.to.keep,]

# Convert back to genlight format
my.data.toploading.gl  <- as.genlight(my.data.mat.toploading.wild)

# Determine distance between individuals
# First use seploc to create a list of ten blocks of loci, each to be computed in parallel
x <- seploc(my.data.toploading.gl, n.block=10, parallel=F)

## Use dist in a lapply loop to compute pairwise distances between individuals for each block
lD <- lapply(x, function(e) dist(as.matrix(e)))

# Obtain general distance matrix by summing these
D <- Reduce("+", lD)

# Plot a neighbor-joining tree using ape's nj function on the distance matrix
plot(nj(D), type="fan", cex=0.7)
# this fails if there are too few markers
# plot(njs(D), type="fan", cex=0.7)
