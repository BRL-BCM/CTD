import pandas as pd
import numpy as np
from collections import defaultdict

# Load the available genes
print("Loading available genes from subtype_spec_genes_available.csv...")
available_genes_df = pd.read_csv('../../data/scRNA_CHD/subtype_spec_genes_available.csv')
available_genes = set(available_genes_df['Gene'].unique())
print(f"Number of available genes: {len(available_genes)}")

# Load the adjacency matrix
print("\nLoading adjacency matrix...")
adj_matrix = pd.read_csv('../../data/scRNA_CHD/adj_TOF.csv')
print(f"Original matrix shape: {adj_matrix.shape}")

# The columns are node names, rows correspond to the same nodes in order
node_names = adj_matrix.columns.tolist()
adj_matrix.index = node_names

# Filter to only include available genes
genes_in_matrix = [gene for gene in node_names if gene in available_genes]
print(f"Number of available genes in matrix: {len(genes_in_matrix)}")

# Subset the adjacency matrix to only available genes
adj_matrix = adj_matrix.loc[genes_in_matrix, genes_in_matrix]
print(f"Filtered matrix shape: {adj_matrix.shape}")

# Calculate node degrees (number of connections for each node)
print("\nCalculating node degrees...")
# For symmetric adjacency matrix, sum each row (or column)
node_degrees = adj_matrix.sum(axis=1)
node_degrees = node_degrees.sort_values(ascending=False)

print(f"\nTop 10 nodes by degree:")
print(node_degrees.head(10))

# Get the top nodes
top_nodes = node_degrees.head(60).index.tolist()  # Get top 60 to work with

# Function to find highly connected subgraph
def find_connected_group(adj_matrix, candidate_nodes, max_size=30):
    """
    Find a highly connected group by starting with the highest degree node
    and adding nodes that have the most connections to the current group
    """
    if len(candidate_nodes) == 0:
        return []

    # Start with the highest degree node
    group = [candidate_nodes[0]]
    remaining = set(candidate_nodes[1:])

    while len(group) < max_size and remaining:
        # For each remaining node, count connections to current group
        best_node = None
        best_connections = 0

        for node in remaining:
            connections = 0
            for group_node in group:
                # Check both directions in adjacency matrix
                if node in adj_matrix.index and group_node in adj_matrix.columns:
                    connections += adj_matrix.loc[node, group_node]
                if group_node in adj_matrix.index and node in adj_matrix.columns:
                    connections += adj_matrix.loc[group_node, node]

            if connections > best_connections:
                best_connections = connections
                best_node = node

        if best_node is None or best_connections == 0:
            break

        group.append(best_node)
        remaining.remove(best_node)

    return group

# Find first group from top nodes
print("\nFinding first highly connected group...")
group1 = find_connected_group(adj_matrix, top_nodes, max_size=30)
print(f"Group 1 size: {len(group1)}")

# Find second group from remaining top nodes
remaining_nodes = [n for n in top_nodes if n not in group1]
print("\nFinding second highly connected group...")
group2 = find_connected_group(adj_matrix, remaining_nodes, max_size=30)
print(f"Group 2 size: {len(group2)}")

# Calculate internal connectivity for each group
def calculate_internal_connections(adj_matrix, group):
    """Calculate total connections within a group"""
    connections = 0
    for i, node1 in enumerate(group):
        for node2 in group[i+1:]:
            if node1 in adj_matrix.index and node2 in adj_matrix.columns:
                connections += adj_matrix.loc[node1, node2]
            if node2 in adj_matrix.index and node1 in adj_matrix.columns:
                connections += adj_matrix.loc[node2, node1]
    return connections

print(f"\nGroup 1 internal connections: {calculate_internal_connections(adj_matrix, group1)}")
print(f"Group 2 internal connections: {calculate_internal_connections(adj_matrix, group2)}")

# Save groups to CSV files
print("\nSaving groups to CSV files...")
group1_df = pd.DataFrame({'x': group1})
group2_df = pd.DataFrame({'x': group2})

group1_df.to_csv('../../data/scRNA_CHD/group1_TOF.csv', index=False)
group2_df.to_csv('../../data/scRNA_CHD/group2_TOF.csv', index=False)

print("\nDone! Groups saved to:")
print("- ../../data/scRNA_CHD/group1_TOF.csv")
print("- ../../data/scRNA_CHD/group2_TOF.csv")
