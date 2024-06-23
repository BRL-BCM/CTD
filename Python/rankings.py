import pandas as pd
import math

def ndcg(rankings:pd.DataFrame, marked:list[str], relevance_scores:dict=None)->float:
    if (relevance_scores is None):
        relevance_scores = {node:1 for node in marked} 
    indices_of_marked: list[int] = rankings.index[rankings['Node_id'].isin(marked)].tolist()
    print(f"Positions of marked nodes are {indices_of_marked}")
    dcg_scores = [1/math.log2(2+rank) for rank in indices_of_marked]#TODO use relevance scores!
    dcg = sum(dcg_scores)
    ideal_dcg = sum([1/math.log2(2+rank) for rank in range(0, len(marked))])#TODO use relevance scores!

    ndcg_score = dcg / ideal_dcg
    print(f"The ndcg score is {ndcg_score}.")
    return ndcg_score




    


    