# To predict the number of cut sites using an in silico digest

# install packages
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("SimRAD") # does this work?
require("SimRAD")

# for GC percent
install.packages("seqinr")
require(seqinr)

# Import the reference genome
# GenBank
rfsq <- ref.DNAseq(FASTA.file = "/Users/wayne/Documents/z-genome-Cgig/GCA_000297895.1_oyster_v9_genomic.fa"
           , subselect.contigs = T
           #, prop.contigs = 0.1
           , prop.contigs = 1
           )

# PAG
rfsq <- ref.DNAseq(FASTA.file = "/Users/wayne/Documents/miller/Moore_oyster/00_resources/Cgigas_genome_PAG/Assembly_1/C_gigas_assembly_1_all_chr.fa"
                   , subselect.contigs = T
                   #, prop.contigs = 0.1
                   , prop.contigs = 1
)


# salmonid control
rfsq <- ref.DNAseq(FASTA.file = "/Users/wayne/Documents/z-genome-Ssalar/UVIC/Ssa_ASM_3.6.fasta"
                   , subselect.contigs = T
                   , prop.contigs = 0.5
                   #, prop.contigs = 1
)


# What is the size of the ref genome?
width(rfsq) 
# 10% contigs = 56596888
# 100% contigs = 557717710

# What is the GC content:
GC(s2c(rfsq)) 
# 10% = 0.33
# 100% = 0.33

##### ddRAD Set your cutters: ####
# https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities
# ddRAD w/ NsiI and MspI
# NsiI = ATGCA/T
cs_5p1 <- "ATGCA"
cs_3p1 <- "T"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"

# ddRAD w/ PstI and MspI
# PstI = CTGCA/G
cs_5p1 <- "CTGCA"
cs_3p1 <- "G"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"

# ddRAD w/ SbfI and MspI
# SbfI = CCTGCA/GG
cs_5p1 <- "CCTGCA"
cs_3p1 <- "GG"

# MspI = C/CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"


# The following enzymes do not have compatible barcodes to the standard (PstI)
# # ddRAD w/ EcoRI and HinfI
# # EcoRI = G/AATTC
# cs_5p1 <- "G"
# cs_3p1 <- "AATTC"
# 
# # HinfI = G/ANTC
# cs_5p2 <- "G"
# cs_3p2 <-  "ANTC"

# # ddRAD w/ EcoRI and MspI
# # EcoRI = G/AATTC
# cs_5p1 <- "G"
# cs_3p1 <- "AATTC"
# 
# # MspI = C/CGG
# cs_5p2 <- "C"
# cs_3p2 <- "CGG"


#### Simulate digest ####
rfsq.dig <- insilico.digest(rfsq, cs_5p1, cs_3p1, cs_5p2, cs_3p2
                              , verbose = T
                              )

# From digest, select only those reads digested by the two enzymes
rfsq.sel <- adapt.select(rfsq.dig, type="AB+BA"
                         , cs_5p1, cs_3p1, cs_5p2, cs_3p2)

# From those selected w/ correct cut sites, select based on size
# 100-250 bp (standard at IBIS, not including barcodes etc)
wid.rfsq.sel <- size.select(rfsq.sel,  min.size = 100, max.size = 250
                            , graph=TRUE, verbose=TRUE)


# # 200-270 bp
# wid.rfsq.sel <- size.select(rfsq.sel,  min.size = 200, max.size = 270
#                             , graph=TRUE, verbose=TRUE)
# # 210-260 bp
# nar.rfsq.sel <- size.select(rfsq.sel,  min.size = 210, max.size = 260
#                             , graph=TRUE, verbose=TRUE)




#### RAD ####
# w/ SbfI = CCTGCA/GG
cs_5p1 <- "CCTGCA"
cs_3p1 <- "GG"

# w/ EcoRI = G/AATTC
cs_5p1 <- "G"
cs_3p1 <- "AATTC"

#### RAD Digest ####
rfsq.dig.RAD <- insilico.digest(rfsq, cs_5p1, cs_3p1
                                , verbose = T)



