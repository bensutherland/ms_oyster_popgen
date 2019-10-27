# Calculate cumulative positions and chromosome midpoints

#### 1.1 Find cumulative positions ####
alignments_w_cum_pos <- alignments
head(alignments_w_cum_pos)

# Sort by chromosome, then position
alignments_w_cum_pos <- alignments_w_cum_pos[order(alignments_w_cum_pos$chr, alignments_w_cum_pos$pos), ]
head(alignments_w_cum_pos)

#### EXTRACT - HERE DID NOT MAP ###
### subset rows that did not map
#unplaced.alignments_w_cum_pos <- alignments_w_cum_pos[which(is.na(alignments_w_cum_pos$chr)) , ]
#dim(unplaced.alignments_w_cum_pos) 
#dim(alignments_w_cum_pos)
## Keep only rows that mapped
#alignments_w_cum_pos <- alignments_w_cum_pos[which(!is.na(alignments_w_cum_pos$chr)) , ]
#dim(alignments_w_cum_pos)
#head(alignments_w_cum_pos)
#tail(alignments_w_cum_pos)
#### END EXTRACT ####

# Add a new column to alignments_w_cum_pos that is a continuous value of total position on the genome
cum.pos <- NULL; cum.pos.string <- NULL; chr <- NULL; identical.row <- NULL; adding.factor <- NULL

for(j in 1:nrow(alignments_w_cum_pos)){
  
  # Find the chromosome for this marker
  chr <- as.character(alignments_w_cum_pos$chr[j])
  print(chr)
  
  # Determine details of this chromosome in the chr.info file
  identical.row <- which(chr.info$chr.names==chr)
  
  # Identify the cumulative position of the preceeding chromosome (unless at first chr)
  if(identical.row > 1){
    adding.factor <- chr.info[identical.row - 1, "cum.lengths"]
  } else {
    adding.factor <- 0
  }
  
  # Add this cumulative position to the current marker
  cum.pos[j] <- alignments_w_cum_pos$pos[j] + adding.factor
}

# Add the cumulative position vector to the alignments_w_cum_pos
alignments_w_cum_pos <- cbind(alignments_w_cum_pos, cum.pos)
head(alignments_w_cum_pos)
tail(alignments_w_cum_pos)


#### 1.2 Find chr midpoints for plotting ####
position.of.chr.label <- NULL

for(k in 1:nrow(chr.info)){
  position.of.chr.label[k] <- (chr.info$chr.lengths[k] / 2) + chr.info$cum.lengths[k - 1]
}

position.of.chr.label[is.na(position.of.chr.label)] <- (chr.info$chr.lengths[1] / 2) # replace the first one due to issue of -1 from 1 position
position.of.chr.label

alignments <- alignments_w_cum_pos

# 
# Output is:
head(alignments)
head(map)
head(chr.info)
