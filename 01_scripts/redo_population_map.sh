#!/bin/bash
# Identify what samples are truly in your study, and rebuild the population map using only these files

ls 04-all_samples/*.fq.gz | awk -F"/" '{ print $2 }' - | sed 's/\.fq\.gz//g' - > 00-info_files/all_retained_samples.txt

grep -f 01-info_files/all_retained_samples.txt 01-info_files/population_map.txt > 01-info_files/population_map_retained.txt

# Then point the populations module to this new population map file
