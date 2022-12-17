import logging
import graph


def single_node_get_node_ranks(n, G, p1, threshold_diff, adj_mat, S=None, num_misses=None, verbose=False, out_dir='',
                               use_labels=False, coords=None):

    """
    Generate single-node node rankings ("fixed" walk)
    This function calculates the node rankings starting from a given perturbed variable in a subset of variables in
    the network.

    Parameters
    ----------
    n : The index (out of a vector of node names) of the node ranking you want to calculate.
    G : A list of probabilities with list names being the node names of the network.
    p1 : The probability that is preferentially distributed between network nodes by the probability diffusion algorithm
    based solely on network connectivity. The remaining probability (i.e., "p0") is uniformly distributed between
    network nodes, regardless of connectivity.
    threshold_diff : When the probability diffusion algorithm exchanges this amount or less between nodes, the algorithm
     returns up the call stack.
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    S : A character vector of node names in the subset you want the network walker to find.
    num_misses : The number of "misses" the network walker will tolerate before switching to fixed length codes for
    remaining nodes to be found.
    verbose : If True, print statements will execute as progress is made. Default is False.
    out_dir : If specified, a image sequence will generate in the output directory specified.
    use_labels : If True, node names will display next to their respective modes in the network. If False, node names
    will not display. Only relevant if out_dir is specified.
    coords : The x and y coordinates for each node in the network, to remain static between images.

    Returns
    -------
    curr_ns - A character vector of node names in the order they were drawn by the probability diffusion algorithm.

    Examples
    --------

    """

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    p0 = 1 - p1
    if not S and (num_misses or out_dir != ''):
        print("You must also supply S if out_dir or num_misses is supplied.")
        return 0

    if verbose:
        logging.debug("Node ranking {} of {}.".format(n, len(G)))

    stop_iterating = False
    start_node = list(G.keys())[n]
    curr_gph = G
    count_misses = 0
    curr_ns = [start_node] # current node set

    if out_dir:
        graph.net_walk_snap_shot(adj_mat=adj_mat, G=G, visited_nodes=curr_ns, S=S, output_dir=out_dir, coords=coords,
                                 img_num=len(curr_ns), use_labels=use_labels)

    while not stop_iterating:
        curr_gph = {node:0 for node in curr_gph.keys()}  # clear probabilities
        base_p = p0 / (len(curr_gph) - len(curr_ns))
        # set unvisited nodes to base_p
        curr_gph.update({node:base_p for node in curr_gph.keys() if node not in curr_ns})
        curr_gph = graph.diffuse_p1(start_node=start_node, G=curr_gph, visited_nodes=curr_ns, p1=p1, adj_mat=adj_mat,
                                    threshold_diff=threshold_diff, verbose=verbose)

        # Sanity check - p1_event should add up to roughly 1
        p1_event = sum(curr_gph.values())
        if abs(p1_event - 1) > threshold_diff:
            extra_prob_to_diffuse = 1 - p1_event
            curr_gph.update({node: 0 for node in curr_gph.keys() if node in curr_ns})
            ind = [node for node in curr_gph.keys() if node not in curr_ns]
            curr_gph.update({node: val + extra_prob_to_diffuse / len(ind) for node, val in curr_gph.items() if node in ind})

        # Set start_node to the node with the max probability in the new curr_gph
        max_prob = max(curr_gph, key=curr_gph.get)

        # Break ties: When there are ties, choose the first of the winners.
        start_node = max_prob

        if out_dir:
            graph.net_walk_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, visited_nodes=curr_ns, S=S, coords=coords
                                     , img_num=len(curr_ns), use_labels=use_labels)

        if S:  # draw until all members of S are found
            if start_node in S:
                count_misses = 0
            else:
                count_misses += 1

            curr_ns.append(start_node)

            if count_misses > num_misses or all(node in curr_ns for node in S):
                stop_iterating = True
        else:
            curr_ns.append(start_node)
            if len(curr_ns) >= len(G):
                stop_iterating = True

        if out_dir:
            graph.net_walk_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, visited_nodes=curr_ns, S=S, coords=coords
                                     , img_num=len(curr_ns), use_labels=use_labels)

    return curr_ns
