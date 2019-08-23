
setwd("/mnt/data/01_moore_oyster_project/example_bayescan_files_2019-08-09/")

require(devtools)

devtools::install_github("thierrygosselin/radiator")
library("radiator")

genomic_converter(data = "DBP_PEN.recode.vcf"
                , strata = "DBP_PEN_strata.txt"
                , output = c('bayescan')
                , filename = 'DPB_TEST_BSCN'
                )

# genomic_converter(data = "13_selection/CHN_CHNF.recode.vcf"
#                   , strata = "13_selection/CHN_CHNF_strata.txt"
#                   , output = c('bayescan')
#                   , filename = 'CHN_TEST_bayescan'
# )
