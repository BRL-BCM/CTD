import json
import random
import pandas as pd

import csv
import numpy as np
import warnings
import scipy.sparse.csgraph
import abc
import math

import scipy.sparse as sp
from scipy.sparse.linalg import eigs, eigsh

def split_80_20(array:list[str], generator_seed:int=42)->tuple[list[str],list[str]]:
    # Calculate the split index (80% of the list)
    split_index = int(len(array) * 0.8)

    # Shuffle the list randomly
    random.seed(generator_seed)
    random.shuffle(array)

    # Split the list into two lists
    anchors = array[:split_index]  # First 80% elements
    targets = array[split_index:]  # Remaining 20% elements

    return (anchors, targets)

def find_centers_of_S(distances_inside_S : pd.DataFrame, verbose=False)->tuple[list[str],float]:
    # Calculate the eccentricity for each node
    eccentricities : pd.Series = np.amax(distances_inside_S, axis=1)
    sorted_eccentricities = eccentricities.sort_values()

    # Find the index of the node with the minimum eccentricity
    min_eccentricity = np.amin(eccentricities) 
    center_node_indices = list(np.where(eccentricities == min_eccentricity)[0])
    center_node_labels = distances_inside_S.columns[center_node_indices]

    if(verbose):
        print(f"Sorted node eccentricities are \n {sorted_eccentricities}")
        print(f"The center of the graph consists of nodes {list(center_node_labels)}.")
        print(f"Graph radius is {min_eccentricity}.")
    return center_node_labels, min_eccentricity

# Returns nodes in S that have a neighbour not in S
def find_inner_boundary_of_S(adjacency_matrix: pd.DataFrame, S_list: list[str])->list[str]:

    degrees = adjacency_matrix.sum(axis = 1)
    adj_S = adjacency_matrix.loc[S_list, S_list]   
    degrees_in_S = adj_S.sum(axis = 1) 

    # Set for storing nodes that are in the inner boundary
    inner_boundary_nodes = set()    

    # Iterate over nodes in S
    for node in S_list:
        if (degrees[node] != degrees_in_S[node]):
            inner_boundary_nodes.add(node)
    return list(inner_boundary_nodes)

class DistanceMetric(abc.ABC):
    def __init__(self, distances_from_S:pd.DataFrame, use_min_distance:bool=True):
        self._distances_from_S = distances_from_S
        self._use_min_distance = use_min_distance
        self._distance_dict = {}
        self._key_points = None

    def set_key_points(self, key_points):
        self._key_points = key_points

    def get_key_points(self)->list[str]:
        return self._key_points

    @abc.abstractmethod
    def get_key_points(self):
        pass

    @abc.abstractmethod
    def calculate_distance(self, Y_label):
        pass 

    @abc.abstractmethod
    def print_info(self):
        pass

class RadialDistanceMetric(DistanceMetric):
    def __init__(self, distances_from_S: pd.DataFrame, use_min_distance:bool=True):
        super().__init__(distances_from_S, use_min_distance)

    def get_key_points(self)->list[str]:
        if self._key_points is None:
            S = list(self._distances_from_S.index)
            distances_inside_S : pd.DataFrame = self._distances_from_S[S]
            self._key_points, self._radius = find_centers_of_S(distances_inside_S)
        return self._key_points

    def calculate_distance(self, Y_label):
        if (Y_label in self._distance_dict.keys()):
            return self._distance_dict[Y_label]
        else:
            if (self._use_min_distance):
                distance = min(list(self._distances_from_S.loc[self.get_key_points(), Y_label]))
            else:
                distance = max(list(self._distances_from_S.loc[self.get_key_points(), Y_label]))
            self._distance_dict[Y_label] = distance
            return distance

    def print_info(self):
        print("Radial distance metric")
        if (self._use_min_distance):
            print(f"Using minimal distance from the centers {self._key_points}.")
        print(f"Radius of S is {self._radius}")

class BoundaryDistanceMetric(DistanceMetric):
    def __init__(self, adjacency_matrix: pd.DataFrame, distances_from_S: pd.DataFrame, use_min_distance:bool=True):
        super().__init__(distances_from_S, use_min_distance)
        self._adjacency_matrix = adjacency_matrix

    # Key points for the boundary distance metric are the inner graph boundary
    def get_key_points(self)->list[str]:
        if self._key_points is None:
            S = list(self._distances_from_S.index)
            # A node belongs to the inner boundary of S if it has a neighbour outside of S 
            self._key_points = find_inner_boundary_of_S(adjacency_matrix=self._adjacency_matrix, S_list=S)
        return self._key_points

    def calculate_distance(self, Y_label):
        if (Y_label in self._distance_dict.keys()):
            return self._distance_dict[Y_label]
        else:
            if (self._use_min_distance):
                distance = min(list(self._distances_from_S.loc[self.get_key_points(), Y_label]))
            else:
                distance = max(list(self._distances_from_S.loc[self.get_key_points(), Y_label]))
            self._distance_dict[Y_label] = distance
            return distance

    def print_info(self):
        print("Boundary distance metric")
        if (self._use_min_distance):
            print(f"Using minimal distance from the boundary formed from nodes {self._key_points}.")
        

def construct_information_hop_matrix(adj_mat : pd.DataFrame, S : list[str], verbose=False):

    # For information distance calculation, only the absolute values of the edge weights are used
    adj_mat = adj_mat.abs()

    # print("(Non-negative) adjacency matrix")
    # print(adj_mat)

    outdegrees = adj_mat.sum(axis=1)
    indegrees = adj_mat.sum(axis=0)
    
    # Print nodes with outdegree equal to 0
    outdeg_zero_nodes = list(np.where(outdegrees == 0)[0])
    outdeg_zero_node_labels = adj_mat.columns[outdeg_zero_nodes]
    if verbose:
        print(f"Nodes with outdegree 0 are {outdeg_zero_node_labels}")

    # Print nodes with indegree equal to 0
    indeg_zero_nodes = list(np.where(indegrees == 0)[0])
    indeg_zero_node_labels = adj_mat.columns[indeg_zero_nodes]
    if verbose:
        print(f"Nodes with indegree 0 are {indeg_zero_node_labels}")

    S_zero_outdeg_labels = [node for node in S if node in outdeg_zero_node_labels]
    S_zero_indeg_labels = [node for node in S if node in indeg_zero_node_labels]
    if (len(S_zero_outdeg_labels) > 0):
        warnings.warn(f"There are {len(S_zero_outdeg_labels)} zero outdegree nodes in S and they are {S_zero_outdeg_labels}.")
    if (len(S_zero_indeg_labels) > 0):
        warnings.warn(f"There are {len(S_zero_indeg_labels)} zero indegree nodes in S and they are {S_zero_indeg_labels}.")
    S_zero_total_deg_labels = list(set(S_zero_outdeg_labels).intersection(S_zero_indeg_labels))
    if (len(S_zero_total_deg_labels) > 0):
        warnings.warn(f"There are {len(S_zero_total_deg_labels)} nodes in S with a total degree of 0 and they are {S_zero_total_deg_labels}.")

    total_deg_zero_node_labels = set(indeg_zero_node_labels).intersection(set(outdeg_zero_node_labels))
    if verbose:
        print(f"Dropping {len(total_deg_zero_node_labels)} nodes with total degree equal to 0 from the dataframe.")
    adj_mat = adj_mat.drop(total_deg_zero_node_labels, axis=0)
    adj_mat = adj_mat.drop(total_deg_zero_node_labels, axis=1)
    
    #outdegrees = outdegrees[outdegrees != 0]
    outdegrees = adj_mat.sum(axis=1)

    # Normalize each element of the adjacency matrix by the degree of the source node
    with np.errstate(divide='ignore', invalid='ignore'):
        normalized_adj_mat = adj_mat.div(outdegrees, axis=0)
    
    if verbose:
        print("Normalized adjacency matrix")
        print(normalized_adj_mat)

    # Create the information hop matrix
    IHM = pd.DataFrame()

    # Apply the transformation -log(normalized_adj_mat[i,j])
    # Use np.where to handle division by zero and logarithm of zero gracefully
    with np.errstate(divide='ignore', invalid='ignore'):
        IHM = np.abs(-np.log(normalized_adj_mat))

    if verbose:
        print("Information hop distance matrix")
        print(IHM)

    return IHM

def calculate_distances_from_S(adj_mat:pd.DataFrame, S:list[str], verbose=False)->tuple[pd.DataFrame, float, pd.DataFrame]:
    S_indices = [adj_mat.index.get_loc(label) for label in S]

    # scipy dijkstra implementation treats edge weights of 0 as non-existent edges, therefore,
    # we need to replace all 0s with 1e-7 (smaller values are treated as 0), and all np.inf values with 0
    replacement_for_zero = 1e-7
    adjusted_adjacency_matrix = adj_mat.replace(0, replacement_for_zero)
    adjusted_adjacency_matrix.replace(np.inf, 0, inplace = True) 

    if (verbose):
        print("Adjusted adjacency matrix:")
        print(adjusted_adjacency_matrix)

    dst_mat = scipy.sparse.csgraph.dijkstra(csgraph=adjusted_adjacency_matrix, directed=True, 
                                            indices=S_indices, return_predecessors=False, unweighted=False, 
                                            limit=np.inf, min_only=False)
    distances_from_S = pd.DataFrame(dst_mat, index=S, columns=adjusted_adjacency_matrix.columns)

    if (verbose):
        print("Distances from S to all nodes:")
        print(distances_from_S)

    distances_inside_S : pd.DataFrame = distances_from_S[S]
    if (verbose):
        print("Distances in S:")
        print(distances_inside_S)

    S_diameter = distances_inside_S.max().max()
    if (verbose):
        print(f"Diameter of S is {S_diameter}")
    if (S_diameter == np.inf):
        warnings.warn("The diameter of S is reported to be infinite!")

    return distances_from_S, S_diameter, distances_inside_S

# def criterion_all_distances_from_S(distances_from_S: pd.DataFrame, Y_label: str, T: float):
#     return distances_from_S[Y_label].max() < T  
    
# def criterion_any_distance_from_S(distances_from_S: pd.DataFrame, Y_label: str, T: float):
#     return distances_from_S[Y_label].min() < T  

# #TODO: change to work with DistanceMetric
# def criterion_radial_distance(distances_from_S: pd.DataFrame, Y_label:str, T: float)->bool:


#     S = list(distances_from_S.index)
#     distances_inside_S : pd.DataFrame = distances_from_S[S]

#     # Only calculate the centers and radius of S the first time, later just reuse it.
#     if (not hasattr(criterion_radial_distance, "centers")):
#         centers, radius = find_centers_of_S(distances_inside_S)
#         criterion_radial_distance.centers = centers
#         criterion_radial_distance.radius = radius
    
#     if (criterion_radial_distance.radius > radial_distance(distances_from_S, criterion_radial_distance.centers, Y_label) ):
#         return True
#     else:
#         return False


# for radial distance, pass RadialDistanceMetric
# TODO: for boundary distance, pass BoundaryDistanceMetric 
# returns a dataframe with node rankings by GBA principle
def get_rankings_by_gba(distances_from_S:pd.DataFrame, distance_metric:DistanceMetric, use_K=False, K=50)->pd.DataFrame:
    #TODO adjust K to depend on the size of S
    ranking_distances = {}
    S = list(distances_from_S.index)
    for label in list(distances_from_S.columns):
        ranking_distance:float = distance_metric.calculate_distance(label)
        ranking_distances[label] = ranking_distance
    ranking_distances = {node: distance for node, distance in ranking_distances.items() if node not in S}
    rankings_df = pd.DataFrame(ranking_distances, index=['Importance']).T
    rankings_df.sort_values(by='Importance', inplace=True, ascending=True)
    rankings_df.reset_index(inplace=True)
    rankings_df.columns = ['Node_id', 'Importance']
    return rankings_df

# def extend_S(adj_mat, S, inclusion_criterion = criterion_radial_distance, verbose=False):
#     distances_from_s, S_diameter, distances_inside_S = calculate_distances_from_S(adj_mat, S)
    
#     # Set the threshold to the information diameter of S
#     # S_radius = S_diameter / 2
#     #extended_S = [label for label in adj_mat.columns if criterion(distances_from_s, label, T)]
    
#     #TODO change T based on used criterion
#     T = S_diameter / 2
#     extended_S = [label for label in adj_mat.columns if inclusion_criterion(distances_from_s, label, T)]
#     if verbose:
#         print(f"Extended S consists of these {len(extended_S)} nodes: {extended_S}")
#     new_S_nodes = [node for node in extended_S if node not in S]
#     if verbose:
#         print(f"Found {len(new_S_nodes)} nodes in extended S that are not in S: {new_S_nodes}")
#     return extended_S

# TODO: reverse direction and do Dijkstra
def calculate_information_distance_from_set(S, Y):
    pass

def calculate_volume(A:list[str], adj_mat:pd.DataFrame)->float:
    submatrix:pd.DataFrame = adj_mat.loc[A, :]
    volume: float = adj_mat[A].to_numpy().sum()
    return volume

def calculate_cut(S:list[str], adj_mat:pd.DataFrame)->float:
    # Ensure S is a list or set of nodes
    S = list(S)
    # Complement of S (all nodes not in S)
    all_nodes = set(adj_mat.index)
    S_complement = list(all_nodes - set(S))
    
    # Extract the submatrix corresponding to edges between S and S_complement
    cut_matrix = adj_mat.loc[S, S_complement]
    
    # Sum the weights of the edges in the cut
    cut_value = cut_matrix.to_numpy().sum()
    
    return cut_value

# Calculate the conductance of a node set S, given by the total weight of edges starting in S and ending out of S divided 
# by the minimum of volumes of S and V\S, where the volume of node set A is given by the sum of weights of all edges starting in A
def calculate_conductance(S:list[str], adj_mat: pd.DataFrame)->float:
    # Calculate the volume of S
    vol_S:float = calculate_volume(S, adj_mat)
    
    # Calculate the volume of the complement of S
    all_nodes = set(adj_mat.index)
    S_complement = list(all_nodes - set(S))
    vol_S_complement:float = calculate_volume(S_complement, adj_mat)
    
    # Calculate the cut between S and S_complement
    cut_value = calculate_cut(S, adj_mat)

    # Compute the conductance
    conductance = cut_value / min(vol_S, vol_S_complement)
    return conductance


def find_fiedler_eigenvalue(matrix):
    """
    Finds the second smallest eigenvalue (Fiedler eigenvalue) of a given sparse matrix.
    Assumes the input matrix is symmetric and all fields are positive.
    
    Parameters:
        matrix (scipy.sparse.csr_matrix): The sparse matrix for which to find the second smallest eigenvalue.

    Returns:
        float: The second smallest eigenvalue.
        numpy.ndarray: The corresponding eigenvector (Fiedler vector).
    """

    
    # Ensure the matrix is sparse and symmetric
    assert sp.issparse(matrix), "Input matrix must be sparse."
    
    # Find the two smallest eigenvalues using the eigs method
    # Since we are interested in the second smallest eigenvalue, we use k=2
    eigenvalues, eigenvectors = eigs(matrix, k=2, which='SR')
    
    # The second smallest eigenvalue is at index 1
    fiedler_eigenvalue = eigenvalues[1]
    fiedler_vector = eigenvectors[:, 1]
    
    return fiedler_eigenvalue, fiedler_vector

def calculate_Cheeger_bounds(matrix):
    fiedler_eigenvalue, _ = find_fiedler_eigenvalue(matrix)
    lower_bound = fiedler_eigenvalue / 2
    upper_bound = math.sqrt(2*fiedler_eigenvalue)
    return lower_bound, upper_bound

# Example usage
if __name__ == "__main__":
    # Create a random sparse matrix for demonstration
    size = 100  # Size of the matrix (100x100)
    density = 0.1  # Approximate density of non-zero elements (10%)
    
    # Generate a random sparse matrix
    random_matrix = sp.random(size, size, density=density)

    random_matrix.data = np.abs(random_matrix.data)
    
    # Symmetrize the matrix (just for demonstration purposes)
    symmetric_matrix = (random_matrix + random_matrix.T) / 2
    
    # Call the function to get the second smallest eigenvalue and vector
    eigenvalue, eigenvector = find_fiedler_eigenvalue(symmetric_matrix)
    
    print(f"Second smallest eigenvalue: {eigenvalue}")
    assert eigenvalue>=0, "Fiedler eigenvalue of G can not be less than zero!"

    lower_bound, upper_bound = calculate_Cheeger_bounds(symmetric_matrix)

    print(f"Cheeger bounds for the isoperimetric number are {lower_bound}, {upper_bound}")