#!/usr/bin/env python3
"""
Pipeline Script for CTD + GBA Analysis
--------------------------------------

This script automates the workflow:
1. Run CTD.py with adjacency matrix + s_module input.
2. Extract F_most_connected_nodes → save as one-column CSV.
3. Run main_gba.py with adj_matrix + generated F_most_connected_nodes CSV.
4. Process resulting *_gba_ranks.csv file:
   - Take top (len(S_nodes) - len(F_nodes)) Node_id entries.
   - Report intersection size with original S_nodes.

Usage:
    python run_pipeline.py --adj_matrix path/to/adj.csv --s_module path/to/s_module.csv
"""

import argparse
import json
import os
import subprocess
import pandas as pd
import glob

# ═══════════════════════════════════════════════════════════════
# 1. Parse arguments
# ═══════════════════════════════════════════════════════════════
def parse_args():
    parser = argparse.ArgumentParser(description="CTD + GBA automated pipeline")
    parser.add_argument("--adj_matrix", required=True, help="Adjacency matrix CSV file")
    parser.add_argument("--s_module", required=True, help="Single-column CSV file with nodes")
    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════
# 2. Run CTD.py and generate JSON output
# ═══════════════════════════════════════════════════════════════
def run_ctd(adj_matrix, s_module):
    """
    Calls CTD.py tool with the provided adjacency matrix and s_module.
    Returns parsed JSON as dictionary.
    """
    json_output_file = s_module.replace('.csv', '.json')#"ctd_output.json"

    subprocess.run([
        "python", "CTD.py",
        "--adj_matrix", adj_matrix,
        "--s_module", s_module,
        "--output_name", json_output_file,
        "-v", "1",
    ], check=True)

    with open(json_output_file, "r") as f:
        ctd_json = json.load(f)

    return ctd_json


# ═══════════════════════════════════════════════════════════════
# 3. Save F_most_connected_nodes → one-column CSV
# ═══════════════════════════════════════════════════════════════
def save_f_nodes_csv(f_nodes, s_module):
    """
    Converts F_most_connected_nodes into a CSV file with header 'x'.
    File is saved with prefix 'f_' before s_module filename.
    """
    # Construct filename
    s_dir = os.path.dirname(s_module)
    s_base = os.path.basename(s_module)
    f_nodes_file = os.path.join(s_dir, "f_" + s_base)

    # Save CSV
    pd.DataFrame(f_nodes, columns=["x"]).to_csv(f_nodes_file, index=False)

    return f_nodes_file


# ═══════════════════════════════════════════════════════════════
# 4. Run main_gba.py with adj_matrix + f_nodes CSV
# ═══════════════════════════════════════════════════════════════
def run_gba(adj_matrix, f_nodes_file):
    """
    Runs main_gba.py with required arguments.
    Returns path to the generated *_gba_ranks.csv file.
    """
    subprocess.run([
        "python", "Python/main_gba.py",
        "--adj_path", adj_matrix,
        "--s_nodes", f_nodes_file
    ], check=True)

    # Locate the newest *_gba_ranks.csv in same directory as f_nodes_file
    out_dir = 'results'  # os.path.dirname(f_nodes_file)
    candidates = glob.glob(os.path.join(out_dir, "*_gba_ranks.csv"))
    if not candidates:
        raise FileNotFoundError("No *_gba_ranks.csv found in output directory")
    else:
        print(f'Created GBA ranks: {candidates[0]}')

    # Pick the most recently modified one
    gba_file = max(candidates, key=os.path.getmtime)
    return gba_file


# ═══════════════════════════════════════════════════════════════
# 5. Post-process GBA results
# ═══════════════════════════════════════════════════════════════
def analyze_results(ctd_json, gba_file):
    """
    Computes top GBA nodes and intersection with S_nodes.
    """
    # Extract nodes from CTD.json
    s_nodes = ctd_json["S_perturbed_nodes"]
    f_nodes = ctd_json["F_most_connected_nodes"]

    # Read GBA results
    gba_nodes = pd.read_csv(gba_file)

    # Select top nodes
    n_top = max(0, len(s_nodes) - len(f_nodes))
    top_gba_nodes = gba_nodes["Node_id"].iloc[:n_top].values

    # Compute intersection
    intersection = set(top_gba_nodes).intersection(s_nodes)

    # Report
    print("══════════════════════════════════")
    print(f"Total S_nodes: {len(s_nodes)}")
    print(f"Total F_nodes: {len(f_nodes)}")
    print(f"Top GBA nodes selected: {len(top_gba_nodes)}")
    print(f"Intersection size: {len(intersection)}")
    print("══════════════════════════════════")


# ═══════════════════════════════════════════════════════════════
# 6. Main driver
# ═══════════════════════════════════════════════════════════════
def main():
    args = parse_args()

    # Step 1: Run CTD
    ctd_json = run_ctd(args.adj_matrix, args.s_module)

    # Step 2: Save F_nodes CSV
    f_nodes_file = save_f_nodes_csv(ctd_json["F_most_connected_nodes"], args.s_module)

    # Step 3: Run GBA
    gba_file = run_gba(args.adj_matrix, f_nodes_file)

    # Step 4: Analyze results
    analyze_results(ctd_json, gba_file)


if __name__ == "__main__":
    main()
# Usage: 
# python run_pipeline.py 
#     --adj_matrix data/example_argininemia/adj.csv 
#     --s_module data/example_argininemia/s_nodes.csv
