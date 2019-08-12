# Plot differentiation statistics across genome (continued)
# note: requires that 01_scripts/gwas_02_load_set_contrasts_info.r already was run

# Summary # your data is here:
names(plotting_data) # qval per datatype
head(plotting_data[[1]])
datatypes # datatypes

#### Find outliers via quantile or 2xSD of mean ####
##NOT CURRENTLY BEING USED ##
all_data <- plotting_data[[1]]
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
table(is.na(perloc_fst$Fst)) # these did not come in with Fst vals, they are probably monomorphic
# save out the results (if want) but all is present in all_data
#data is in all_data


##### PLOT ON CHROMOSOMES #####
# Go to 01_scripts/sam_to_gwas.R
