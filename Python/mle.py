import pandas as pd
import numpy as np
import math


def get_encoding_length(bs, G, pvals=None, pt_id=None, pen_p_val=False):

    """
    Minimum encoding length

    This function calculates the minimum encoding length associated with a subset of variables given a background
    knowledge graph.

    Parameters
    ----------
    bs : A list of bitstrings associated with a given patient's perturbed variables.
    G : A dictionary of probabilities with keys being the node names of the background graph.
    pvals : The matrix that gives the perturbation strength significance for all variables (columns)
    for each patient (rows).
    pt_id : The row name in pvals corresponding to the patient you specifically want encoding information for.

    Returns
    -------
    res : A DataFrame that for every bitstring provided in bs input parameter, a row is returned with the
    following data: the patientID; the bitstring evaluated where T denotes a hit and 0 denotes a miss; the subsetSize,
    or the number of hits in the bitstring; the individual p-values associated with the variable's perturbations,
    delimited by '/'; Shannon's entropy, IS.null; the minimum encoding length IS.alt; and IS.null-IS.alt,the d.score.

    Examples
    --------

    """

    res = []

    for k in range(len(bs)):
        row = {}
        opt_bs = bs[k]
        mets_k = [d[0] for d in opt_bs if d[1] == 1]  # node rankings that found at least k nodes
        found = sum([d[1] for d in opt_bs])  # number of ranks that found at least k nodes
        not_found = k - found + 1
        if pen_p_val:
            e = math.comb(len(G), (not_found + 1)) + len(opt_bs) - 1  # encoding length
        else:
            e = (not_found + 1) * np.log2(len(G)) + len(opt_bs) - 1  # encoding length
        opt_bs_tmp = ''.join(['T' if d[1] == 1 else '0' for d in opt_bs])  # T denotes a hit and 0 denotes a miss

        if pt_id and pvals is not None:
            row['patientID'] = pt_id
            row['varPvalue'] = '/'.join([str(pval) for pval in pvals.loc[pt_id, mets_k]]) 

        row['optimalBS'] = opt_bs_tmp  # optimal bitstring
        row['subsetSize'] = k+1  # assume k=0 corresponds to subset of size 1
        row['opt.T'] = found
        if pen_p_val:
            row['IS.null'] = np.float32(math.comb(len(G), (k + 1)))
        else:
            row['IS.null'] = np.float32(np.log2(len(G)) * (k + 1))  # Shannon's entropy
        row['IS.alt'] = np.float32(e)
        row['d.score'] = round(row['IS.null'] - e, 3)

        res.append(row)

    return pd.DataFrame(res)


def get_pt_bs_by_k(S, ranks, num_misses=None):

    """
    Generate patient-specific bitstrings

    This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with a network walker which tries to
    find all nodes in a given subset, S, in a given network, G.

    Parameters
    ----------
    S : List of node names describing the node subset to be encoded
    ranks : A dictionary of node ranks calculated over all possible nodes, starting with each node in subset of interest.
    num_misses : The number of misses tolerated by the network walker before path truncation occurs.

    Returns
    -------
    pt_by_k : A dictionary of bitstrings.

    Examples
    --------

    """

    if not num_misses:
        num_misses = np.log2(len(ranks))

    ranks2 = {k: v for k, v in ranks.items() if k in S}

    pt_bit_string = {}

    for node in S:
        miss = 0
        thresh = len(ranks2[node])
        pt_bit_string[node] = []

        for ii in range(len(ranks2[node])):
            ind_t = int(ranks2[node][ii] in S)
            if not ind_t:
                miss += 1
                if miss >= num_misses:
                    thresh = ii
                    break
            else:
                miss = 0

            pt_bit_string[node].append((ranks2[node][ii], ind_t))

        pt_bit_string[node] = pt_bit_string[node][0:thresh]
        ind = [i for i in range(len(pt_bit_string[node])) if pt_bit_string[node][i][1] == 1]
        pt_bit_string[node] = pt_bit_string[node][0:ind[-1] + 1]

    # For each k, find the ranks that found at least k in S. Which node
    # rankings, from different start nodes, found the first k soonest?
    pt_by_k = []
    for k in range(1, len(S) + 1):
        pt_by_k_tmp = pt_bit_string
        best_ind = {}

        # Which found at least k nodes
        for node in S:
            best_ind[node] = sum([d[1] for d in pt_bit_string[node]])

        found_k = [node for node, s in best_ind.items() if s >= k]

        if found_k:
            pt_by_k_tmp = {node: v for node, v in pt_by_k_tmp.items() if node in found_k}
            # Which found the most nodes soonest
            best_ind = {}
            for node in pt_by_k_tmp:
                best_ind[node] = sum(
                    [i + 1 for i in range(len(pt_by_k_tmp[node])) if pt_by_k_tmp[node][i][1] == 1][0:k + 1])

            lst = pt_by_k_tmp[min(best_ind, key=best_ind.get)]
            max_ind = max([i for i in range(len(lst)) if lst[i][1] == 1][0:k + 1])
            pt_by_k.append(lst[0:max_ind + 1])

        else:
            pt_by_k.append(pt_by_k_tmp[max(best_ind, key=best_ind.get)])

    return pt_by_k
