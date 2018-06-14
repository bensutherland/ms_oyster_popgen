#!/bin/bash
# Run this script from within the stacks_workflow main repo
# This will depend that all of the following individual scripts have been edited to the desired settings

# pstacks
./00-scripts/stacks_1b_pstacks.sh

# cstacks
./00-scripts/stacks_2_cstacks.sh

# sstacks
./00-scripts/stacks_3_sstacks.sh

# populations (round 1)
./00-scripts/stacks_4_populations.sh

# rxstacks
./00-scripts/stacks_5b_rxstacks.sh

# cstacks (2)
./00-scripts/stacks_6_cstacks_rx.sh

# sstacks (2)
./00-scripts/stacks_7_sstacks_rx.sh

# populations (round 2)
./00-scripts/stacks_8_populations_rx.sh

echo "Done!"
