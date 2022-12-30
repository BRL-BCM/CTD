import argparse
import time
import subprocess
import zipfile
import json
import pandas as pd
import random
import os

cpu_count = os.cpu_count()
p = argparse.ArgumentParser(description="Connect The Dots - R vs Python")

p.add_argument("--out_path", help="Path to the directory for storing outputs.", default='.')
p.add_argument("--adj_matrix", help="CSV with adjacency matrix. If none is given, Kegg graph will be used.", default='')
p.add_argument("--s_module",
               help="Comma-separated list or path to CSV of graph G nodes to consider when searching for the most "
                    "connected subgraph. If 'random_N' is specified, N random nodes will be selected.")
p.add_argument("--num_processes", help="Number of worker processes to use for parallelization. Default is to use the "
                                       "number returned by os.cpu_count()", default=cpu_count, type=int)

if __name__ == '__main__':

    argv = p.parse_args()
    messages = []

    out_path = argv.out_path
    np = argv.num_processes
    arguments_dict = argv.__dict__
    arguments_dict.pop('out_path', None)
    arguments_dict.pop('num_processes', None)

    try:
        os.mkdir(f'{out_path}')
    except FileExistsError:
        pass

    try:
        os.mkdir(f'{out_path}/python_results')
    except FileExistsError:
        pass

    try:
        os.mkdir(f'{out_path}/R_results')
    except FileExistsError:
        pass

    if 'adj_matrix' in arguments_dict:
        if not arguments_dict['adj_matrix']:
            with zipfile.ZipFile('kegg_s_module.zip', 'r') as zip_ref:
                zip_ref.extractall(out_path)
            arguments_dict['adj_matrix'] = f'{out_path}/kegg_large_graph.csv'

    if 's_module' in arguments_dict:
        if arguments_dict['s_module']:
            if arguments_dict['s_module'].startswith('random_'):
                nodes = pd.read_csv(arguments_dict['adj_matrix'], nrows=2).columns
                N = int(arguments_dict['s_module'].split('_')[-1])
                arguments_dict['s_module'] = ','.join(random.sample(list(nodes), N))
            if 'csv' not in arguments_dict['s_module']:
                arguments_dict['s_module'] = '"' + arguments_dict['s_module'] + '"'

    cmdline = ' '.join([f'--{arg} {arguments_dict[arg]}' for arg in arguments_dict if arguments_dict[arg]])
    messages.append(f'Command line arguments: {cmdline}\n')
    out_name_R = f'{out_path}/R_results/out.json'
    out_name_py = f'{out_path}/python_results/out.json'
    messages.append('\n##############################################################\n\n')

    print('------------- Running R version of CTD -------------')
    start_R = time.time()
    subprocess.call(f"Rscript CTD.r {cmdline} --output_name {out_name_R}", shell=True)
    end_R = time.time()
    time_R = end_R - start_R

    print(f'------------- Running Python version of CTD, num_processes = {np} -------------')
    start_py = time.time()
    subprocess.call(f"python CTD.py {cmdline} --output_name {out_name_py} --num_processes {np}", shell=True)
    end_py = time.time()
    time_py = end_py - start_py

    messages.append(f'Execution time - R: {round(time_R, 3)}\n')
    messages.append(f'Execution time - Python, num_processes = {np}: {round(time_py, 3)}\n')

    with open(out_name_R.replace('.json', '_ranks.json'), 'r') as f:
        ranks_r = json.load(f)

    with open(out_name_R, 'r') as f:
        d_r = json.load(f)

    with open(out_name_py.replace('.json', '_ranks.json'), 'r') as f:
        ranks_py = json.load(f)

    with open(out_name_py, 'r') as f:
        d_py = json.load(f)

    messages.append(f'Ranks identical: {ranks_py == ranks_r}\n')

    out_dicts_identical = True

    for k in d_r:
        if type(d_r[k]) == float:
            diff = abs(d_r[k] - d_py[k])
            out_dicts_identical = diff < 1e-10
        elif type(d_r[k]) != list and type(d_py[k]) == list:
            out_dicts_identical = [d_r[k]] == d_py[k]
        else:
            out_dicts_identical = d_r[k] == d_py[k]

        if not out_dicts_identical:
            break

    messages.append(f'Out dicts identical: {out_dicts_identical}\n')

    with open(f'{out_path}/comparison.txt', 'w') as f:
        f.writelines(messages)

    print(f'\nDone. Comparison results written to {out_path}/comparison.txt file.')

