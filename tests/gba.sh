#!/bin/bash

# Array of different graph files
graphs=("1-Matrix_ML_NET_Female_adj.csv" "1-Matrix_ML_NET_Male_adj.csv" "1-Matrix_ML_NET_Male_Female_adj.csv")
# graphs=("1-Matrix_ML_NET_Female_adj.csv")

# Outer loop for each graph file
for graph in "${graphs[@]}"
do
    # Inner loop for the different trajectory files
    for i in 0 1 2 3 4 5 8 9 10
    do
        python ../Python/main_gba.py --s_nodes "../data/Male_female_and_trajectories/trajectory_${i}.csv" --adj_path "../data/Male_female_and_trajectories/${graph}" &
    done
    wait  # Wait for all background processes in the inner loop to finish
done

echo DONE!

