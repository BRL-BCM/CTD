#!/usr/bin/env python3
import os
import json
import subprocess
import tempfile
import argparse
import math

def run_ctd(adj_matrix, s_module):
    """
    Calls CTD.py tool with the provided adjacency matrix and s_module.
    Returns parsed JSON as dictionary.
    """
    json_output_file = s_module.replace('.csv', '.json')

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

def read_nodes_from_csv(path):
    """Reads one-column CSV of nodes and returns as list of strings."""
    nodes = []
    with open(path, "r") as f:
        for line in f:
            node = line.strip()
            if node and node.lower() != "x":  # ignore header
                nodes.append(node)
    return nodes

def write_nodes_to_csv(nodes, path):
    """Writes list of nodes into one-column CSV with header 'x'."""
    with open(path, "w") as f:
        f.write("x\n")
        for n in nodes:
            f.write(f"{n}\n")

def calculate_penalized_score(p_union, f_union, f_a, f_b):
    """
    Calculate penalized p-value score for union analysis.

    Formula:
    Score_union = -log(p_union) * (|F_union| / max(|F_A|, |F_B|))

    This score combines:
    - Statistical significance: -log(p_union) - higher when p-value is smaller
    - Improvement ratio: |F_union| / max(|F_A|, |F_B|) - how much union improves over best individual

    Parameters:
    -----------
    p_union : float
        P-value for the union analysis
    f_union : int
        Number of F nodes found in union
    f_a : int
        Number of F nodes found in set A
    f_b : int
        Number of F nodes found in set B

    Returns:
    --------
    float : Penalized score, or None if calculation not possible

    Notes:
    ------
    - Returns None if p_union is 0, NA, or invalid
    - Returns None if max(|F_A|, |F_B|) is 0 (to avoid division by zero)
    - Higher scores indicate better combinations (more significant + more improvement)
    """
    try:
        # Validate inputs
        if p_union is None or p_union == "NA" or p_union <= 0:
            return None

        f_max = max(f_a, f_b)
        if f_max == 0:
            return None

        # Calculate score
        # -log(p) increases as p decreases (more significant)
        # |F_union| / max(|F_A|, |F_B|) increases with improvement
        score = -math.log(p_union) * (f_union / f_max)

        return score

    except (ValueError, TypeError, ZeroDivisionError):
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Run CTD on multiple s_modules and their union."
    )
    parser.add_argument(
        "--adj_matrix", 
        required=True, 
        help="Path to adjacency matrix CSV file"
    )
    parser.add_argument(
        "--s_module", 
        required=True, 
        nargs="+", 
        help="Paths to one or more s_module CSV files"
    )
    args = parser.parse_args()

    adj_matrix = args.adj_matrix
    s_modules = args.s_module

    results = {}
    all_nodes = set()

    # Run CTD on each S
    for s_module in s_modules:
        print(f"Running CTD for {s_module}...")
        ctd_result = run_ctd(adj_matrix, s_module)
        results[s_module] = ctd_result

        # Collect nodes for union
        nodes = read_nodes_from_csv(s_module)
        all_nodes.update(nodes)

    # Write union to temporary file
    union_file = tempfile.mktemp(suffix="_union.csv")
    write_nodes_to_csv(sorted(all_nodes), union_file)

    print("Running CTD for union of all S modules...")
    union_result = run_ctd(adj_matrix, union_file)
    results["UNION"] = union_result

    # Calculate penalized score for union
    penalized_score = None
    if len(s_modules) >= 2:
        # Get F node counts from first two modules (assuming we're comparing two sets)
        f_counts = []
        for s_module in s_modules[:2]:
            F_nodes = results[s_module].get("F_most_connected_nodes", [])
            f_counts.append(len(F_nodes))

        # Get union metrics
        p_union = union_result.get("p_value", "NA")
        F_union_nodes = union_result.get("F_most_connected_nodes", [])
        f_union = len(F_union_nodes)

        # Calculate penalized score
        if len(f_counts) >= 2:
            penalized_score = calculate_penalized_score(
                p_union,
                f_union,
                f_counts[0],
                f_counts[1]
            )

    # Reporting
    print("\n=== CTD Results Summary ===")
    for key, res in results.items():
        p_val = res.get("p_value", "NA")
        F_nodes = res.get("F_most_connected_nodes", [])
        print(f"\nSubset: {key}")
        print(f"  p-value: {p_val}")
        print(f"  |F|: {len(F_nodes)}")
        print(f"  F nodes (first 10): {F_nodes[:10]}{'...' if len(F_nodes) > 10 else ''}")

    # Print penalized score
    if penalized_score is not None:
        print(f"\n=== Penalized P-Value Score ===")
        print(f"  Score_union = -log(p_union) * (|F_union| / max(|F_A|, |F_B|))")
        print(f"  Score: {penalized_score:.4f}")
        print(f"\n  Interpretation:")
        print(f"    - Higher score = Better combination (more significant + more improvement)")
        print(f"    - Combines statistical significance with improvement ratio")
    else:
        print(f"\n=== Penalized P-Value Score ===")
        print(f"  Score: NA (insufficient data or invalid p-value)")

if __name__ == "__main__":
    main()
