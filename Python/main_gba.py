import csv
import argparse as ap
import os

import numpy as np
import pandas as pd

import gba_analysis

def main(args=None):

    parser = ap.ArgumentParser(
        description="A script that performs Guilt by association.")
    parser.add_argument(
        "--adj_path", help="Path to adjacency matrix CSV with columns, no index", type=str, required=True
    )
    parser.add_argument(
        "--s_nodes", help="Path to S nodes CSV", type=str, required=True
    )
    parser.add_argument(
        "-m",
        "--metric",
        help="radial or boundary",
        type=str,
        default="boundary",
    )
    parser.add_argument(
        "-o",
        "--outpath",
        help="Output path",
        type=str,
        default="./results",
    )

    args = parser.parse_args()
    metric = args.metric

    adj_df = pd.read_csv(args.adj_path, skip_blank_lines=True)
    adj_df.index = adj_df.columns
    s_nodes = pd.read_csv(args.s_nodes).astype(str).squeeze().values

    # DEBUG: TODO delete this!
    # has_non_zero = (adj_df.loc['6905'] != 0).any()

    # Construct information hop matrix used for information distance-based neighbourhood extension
    IHM = gba_analysis.construct_information_hop_matrix(adj_df, s_nodes)

    # Check the conductance of S
    conductance = gba_analysis.calculate_conductance(s_nodes, adj_df)
    print(f"Conductance of S is {conductance}.")

    distances_from_S, S_diameter, distances_inside_S = gba_analysis.calculate_distances_from_S(IHM, s_nodes)
    if (metric == "radial"):
        # Calculate rankings by radial distance
        distance:gba_analysis.RadialDistanceMetric = gba_analysis.RadialDistanceMetric(distances_from_S, use_min_distance=True)            
    elif (metric == "boundary"):
        # Calculate rankings by boundary distance
        distance:gba_analysis.BoundaryDistanceMetric = gba_analysis.BoundaryDistanceMetric(adj_df, distances_from_S, use_min_distance=True)

    rankings_df = gba_analysis.get_rankings_by_gba(distances_from_S, distance_metric=distance)
    os.makedirs(args.outpath, exist_ok=True)
    
    outfname = os.path.join(args.outpath, os.path.basename(args.s_nodes).replace('.csv', f'_{np.round(conductance, 3)}_gba_ranks.csv'))
    rankings_df.to_csv(outfname, index=False)

if __name__ == "__main__":
    main()

    #python Python/main_gba.py --s_nodes data/Male_female_and_trajectories/trajectory_0.csv --adj_path data/Male_female_and_trajectories/1-Matrix_ML_NET_Male_Female_adj.csv