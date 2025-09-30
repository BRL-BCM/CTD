#!/usr/bin/env python3
import os
import json
import subprocess
import tempfile
import argparse

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

    # Reporting
    print("\n=== CTD Results Summary ===")
    for key, res in results.items():
        p_val = res.get("p_value", "NA")
        F_nodes = res.get("F_most_connected_nodes", [])
        print(f"\nSubset: {key}")
        print(f"  p-value: {p_val}")
        print(f"  |F|: {len(F_nodes)}")
        print(f"  F nodes (first 10): {F_nodes[:10]}{'...' if len(F_nodes) > 10 else ''}")

if __name__ == "__main__":
    main()
