import pandas as pd
import csv
import numpy as np
import subprocess
import os
from dataclasses import make_dataclass
import time
import sys

# Get the root directory
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the Python directory to sys.path
python_dir = os.path.join(root_dir, 'Python')
sys.path.append(python_dir)

import gba_analysis, rankings

def write_list_to_csv(lst:list[str], path:str, header_column:str='x'):
    with open(path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([header_column])
        # Write each string as a new row in the CSV file
        for string in lst:
            writer.writerow([string])  # Each string is written as a single-element list

def read_S(s_path:str)->list[str]:
    with open(s_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)
        S_nodes = [row for row in csv_reader]

    S_nodes = np.array(S_nodes)
    S_nodes = S_nodes.flatten()
    return S_nodes

class NextID_GBATest :
    # If split_precalculated, the passed anchors and targets are used.
    # Otherwise, the split is calculated.
    # metric_description: radial_min, radial_max, boundary_min, boundary_max
    def __init__(self, disease_name:str=None, metric_description:str="radial_min", adj_matrix_path:str=None, s_path:str=None, 
                 split_precalculated:bool=False, anchors:list[str]=None, targets:list[str]=None, split_seed:int=42) -> None:
        self._disease_name = disease_name
        self._metric, self._criterion = tuple(metric_description.split("_"))
        print(f"Metric used: {self._metric} with {self._criterion} distance as criterion.")
        self._adj_matrix_path = adj_matrix_path

        self._adj_df: pd.DataFrame = pd.read_csv(adj_matrix_path, skip_blank_lines=True)
        self._adj_df.index = self._adj_df.columns

        self._S_nodes = read_S(s_path=s_path)

        if (not split_precalculated):
            print("Generating an 80-20 split of S to get anchors and targets.")
            self._anchors, self._targets = gba_analysis.split_80_20(self._S_nodes, split_seed)
        else:
            print("Using precalculated anchors and targets.")
            self._anchors = anchors
            self._targets = targets

        self._outfname = "nextid_rankings_"+disease_name+"_"+metric_description+".csv"
        self._out_path = f"results/{disease_name}/"+self._outfname

    # In order to reuse the test with the same adj and S, we need to re-split S into anchors and targets
    def reset_split(self, seed=42):
        
        self._anchors, self._targets = gba_analysis.split_80_20(self._S_nodes, seed)

    def run(self):
        # Construct information hop matrix used for information distance-based neighbourhood extension
        IHM = gba_analysis.construct_information_hop_matrix(self._adj_df, self._anchors)

        #gba_analysis.extend_S(IHM, self._anchors, gba_analysis.criterion_radial_distance)

        # distance = None
        distances_from_S, S_diameter, distances_inside_S = gba_analysis.calculate_distances_from_S(IHM, self._anchors)
        if (self._metric == "radial"):
            # Calculate rankings by radial distance
            distance:gba_analysis.RadialDistanceMetric = gba_analysis.RadialDistanceMetric(distances_from_S, use_min_distance=(self._criterion=="min"))            
        elif (self._metric == "boundary"):
            # Calculate rankings by boundary distance
            distance:gba_analysis.BoundaryDistanceMetric = gba_analysis.BoundaryDistanceMetric(self._adj_df, distances_from_S, use_min_distance=(self._criterion=="min"))

        self._rankings_df:pd.DataFrame = gba_analysis.get_rankings_by_gba(distances_from_S, distance_metric=distance)
        self._ndcg_score = rankings.ndcg(self._rankings_df, self._targets)

        # Write results
        results_directory_path = f"results/{self._disease_name}/seed_{seed}"
        if not os.path.exists(results_directory_path):
            os.makedirs(results_directory_path)
        out_path = results_directory_path+"/"+self._outfname
        self._rankings_df.to_csv(out_path, index=False)
        print(f"Finishing NextID test for {self._disease_name}, seed {seed}, metric {self._metric}_{self._criterion}.")


Disease = make_dataclass("Disease", [("name", str), ("adj_path", str), ("s_path", str)])

arthritis = Disease("arthritis", "data/arthritis/Arthritis_large_adj.csv", "data/arthritis/Arthritis_all_genes_large.csv")
asthma = Disease("asthma","data/asthma/Asthma_large_adj.csv", "data/asthma/Asthma_all_genes_large.csv")
chron_pulmo = Disease("chronic_obstructive_pulmonary_disease","data/chronic_obstructive_pulmonary_disease/Chronic_Obstructive_Pulmonary_Disease_large_adj.csv", "data/chronic_obstructive_pulmonary_disease/Chronic_Obstructive_Airway_Disease_all_genes_large.csv")
dilated_cardiomyopath = Disease("dilated_cardiomyopathy", "data/dilated_cardiomyopathy/Dilated_Cardiomyopathy_large_adj.csv", "data/dilated_cardiomyopathy/Cardiomyopathy_Dilated_all_genes_large.csv")
breast_carcinoma = Disease("invasive_breast_carcinoma", "data/invasive_breast_carcinoma/Invasive_Breast_Carcinoma_large_adj.csv", "data/invasive_breast_carcinoma/Breast_Carcinoma_all_genes_large.csv" )
lung_adenocarcinoma = Disease("lung_adenocarcinoma", "data/lung_adenocarcinoma/Lung_Adenocarcinoma_large_adj.csv", "data/lung_adenocarcinoma/Carcinoma_of_lung_all_genes_large.csv")
psoriasis = Disease("psoriasis", "data/psoriasis/Psoriasis_large_adj.csv", "data/psoriasis/Psoriasis_all_genes_large.csv")
ulcerative_colitis = Disease("ulcerative_colitis", "data/ulcerative_colitis/Ulcerative_Colitis_large_adj.csv", "data/ulcerative_colitis/Ulcerative_Colitis_all_genes_large.csv")

# disease_list = [breast_carcinoma]
disease_list = [arthritis, asthma, chron_pulmo, dilated_cardiomyopath, breast_carcinoma,
               lung_adenocarcinoma, psoriasis, ulcerative_colitis]
diseases_df = pd.DataFrame(disease_list)

split_seeds = [42, 45, 55, 100, 420] #TODO generate and use more splits?
distance_metrics_used = ["radial_min", "boundary_min"]
#distance_metrics_used = ["boundary_min"]
results_list: list[dict] = []

for index, row in diseases_df.iterrows(): # loop over diseases
    for seed in split_seeds:
        start_time = time.time()

        disease_name = row['name'] #"arthritis"
        path_disease_adj = row['adj_path'] #"data/arthritis/Arthritis_large_adj.csv"
        path_disease_S = row['s_path'] #"data/arthritis/Arthritis_all_genes_large.csv"
        test_name = f"{disease_name}_{seed}"

        print("_________________________________________________________________________")
        print(f"Test suit {index}: {disease_name}, seed {seed}")

        S = read_S(path_disease_S)
        anchors, targets = gba_analysis.split_80_20(S, generator_seed=seed)

        results_directory_path = f"results/{disease_name}/seed_{seed}"
        if not os.path.exists(results_directory_path):
            os.makedirs(results_directory_path)

        anchors_path = results_directory_path+f"/anchors_{seed}.csv"
        write_list_to_csv(anchors, anchors_path)

        targets_path = results_directory_path+f"/targets_{seed}.csv"
        write_list_to_csv(targets, targets_path)

        # Uncomment if you want to also run CTD_GBA1 on the data

        # print(f"Running CTD tests for {disease_name}")
        # ctd_start_time = time.time()
        # ctd_args = ["--adj_matrix", path_disease_adj, "--s_module", anchors_path, "--output_name", f"{test_name}.json"]
        # command = ['python', 'CTD.py'] + ctd_args
        # ctd_subprocess = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # stdout, stderr = ctd_subprocess.communicate()

        # dtype_dict = {'Node_id': 'str', 'Importance': 'float'}
        # ctd_gba_rankings_df:pd.DataFrame = pd.read_csv(f"{test_name}_gba.csv", dtype=dtype_dict)
        # ndcg_ctd_gba = rankings.ndcg(ctd_gba_rankings_df, targets)

        # ctd_end_time = time.time()
        # ctd_time = ctd_end_time - ctd_start_time
        # print(f"CTD test completed and took {ctd_time} seconds.\n")

        print(f"Running NextID_GBA tests for {disease_name}.")
        for metric in distance_metrics_used:
            
            nextid_start_time = time.time()
            test = NextID_GBATest(disease_name=disease_name, metric_description=metric, adj_matrix_path=path_disease_adj, s_path=path_disease_S, 
                                        split_precalculated=True, anchors=anchors, targets=targets)
            test.run()
            nextid_end_time = time.time()
            nextid_time = nextid_end_time - nextid_start_time
            

            # Append to results dataframe
            new_test_result = {'test_name':test_name+"_" + metric, 
                      # 'ndgc_ctd':ndcg_ctd_gba,
                      # 'time_ctd':ctd_time, 
                      'ndgc':test._ndcg_score, 
                      'time':nextid_time, 
                      # 'time_total':test_execution_time, 
                      'split_name':f"generated_{seed}"}
            results_list.append(new_test_result)
            print(f"NextID_GBA test for {disease_name} using {metric} completed and took {nextid_time} seconds.\n")

        end_time = time.time()
        test_execution_time = end_time - start_time
               
        print(f"Test suit {index} took {test_execution_time} seconds to complete.")

results_dtype_dict = {'test_name':'str', 
                      'ndgc':'float', 
                      'time':'float',  
                      'split_name':'str'}
results_columns = ['test_name', 'ndgc', 'time', 'split_name']
# results_columns = ['test_name', 'ndgc_ctd', 'time_ctd', 'ndgc_nextid', 'time_nextid', 'time_total', 'split_name']
benchmark_results_df:pd.DataFrame = pd.DataFrame(data=results_list, columns=results_columns)
benchmark_results_df.astype(results_dtype_dict)
benchmark_results_df.to_csv("results/benchmark_results.csv", index=False, header=True)