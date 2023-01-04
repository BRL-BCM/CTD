import logging
from datetime import datetime

import igraph
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def net_walk_snap_shot(adj_mat, G, output_dir, visited_nodes, S, coords, img_num=1, use_labels=True):

    """
    Capture the current location of a network walker.

    A network walker steps towards the node that inherited the highest probability from the last node that it
    stepped into.

    Parameters
    ----------
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    G : A dictionary of probabilities, with keys being the node names in the network.
    output_dir : The local directory at which you want still PNG images to be saved.
    visited_nodes : A list of node names, storing the history of previous draws in the node ranking.
    S : A list of node names in the subset you want the network walker to find.
    coords : The x and y coordinates for each node in the network, to remain static between images.
    img_num : The image number for this snapshot. If images are being generated in a sequence, this serves as an
    iterator for file naming.
    use_labels : If True, node names will display next to their respective nodes in the network. If False, node names
    will not display.

    Returns
    -------

    Examples
    --------

    """

    ig = igraph.Graph.Weighted_Adjacency(adj_mat, mode='max')

    ig.vs['color'] = ['white' if ig.vs['name'][i] not in S else 'green' for i in range(len(G))]

    if use_labels:
        ig.vs['label'] = ig.vs['name']
    else:
        ig.vs['label'] = ['' for _ in range(len(ig.vs['name']))]

    ig.vs['cluster'] = [1 if node in visited_nodes else 0 for node in ig.vs['name']]
    communities = igraph.VertexClustering.FromAttribute(ig, 'cluster')
    palette = igraph.PrecalculatedPalette(['dark red', 'white'])

    out_fig_name = f'{output_dir}/netWalkMovie{S.index(visited_nodes[0])}_{img_num}.png'

    visual_style = {
        'vertex_color': ig.vs['color'],
        'vertex_label': ig.vs['name'],
        'vertex_size': [0.3 + round(50 * v, 0) for v in G.values()],
        'edge_width': 3 * np.abs(ig.es['weight']),
        'edge_arrow_size': 0.01,
        'mark_groups': True,
        'palette': palette,
        'layout': coords

    }

    fig, ax = plt.subplots(figsize=(20, 20))

    igraph.plot(communities, target=ax, **visual_style)

    ax.title.set_text('{{{}}}'.format(''.join(['1' if node in S else '0' for node in visited_nodes])))
    jackpot_patch = mpatches.Patch(color='#02fe00', label='Jackpot Nodes')
    drawn_patch = mpatches.Patch(color='#c48080', label='Drawn Nodes')
    ax.legend(handles=[jackpot_patch, drawn_patch], loc='upper right')

    plt.savefig(out_fig_name, dpi=100)


def diffusion_snap_shot(adj_mat, G, output_dir, p1, start_node, visited_nodes, coords, recursion_level=1):
    """
    Capture the current state of probability diffusion

    Recursively diffuse probability from a starting node based on the connectivity in a network, G, where the
    probability represents the likelihood that a variable will be influenced by a perturbation in the starting node.

    Parameters
    ----------
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    G : A dictionary of probabilities, with keys being the node names in the network.
    output_dir : The local directory at which you want still PNG images to be saved.
    p1 : The probability being dispersed from the starting node, startNode, which is preferentially distributed between
    network nodes by the probability diffusion algorithm based solely on network connectivity.
    start_node : The first variable drawn in the node ranking, from which p1 gets dispersed.
    visited_nodes : A character vector of node names, storing the history of previous draws in the node ranking.
    coords : The x and y coordinates for each node in the network, to remain static between images.
    recursion_level : The current depth in the call stack caused by a recursive algorithm.

    Returns
    -------

    Examples
    --------

    """

    ig = igraph.Graph.Weighted_Adjacency(adj_mat, mode='max')

    G = {node: val for node, val in G.items() if node in ig.vs['name']}

    ig.vs['color'] = ['red' if ig.vs['name'][i] in visited_nodes else 'blue' for i in range(len(G))]

    vals = {node: 0 for node in ig.vs['name']}
    vals.update({node: val for node, val in G.items()})

    ig.vs['label'] = [f'{node}:{val:.2f}' for node, val in vals.items()]

    curr_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S.%f")[:-3]

    ig.vs['cluster'] = [1 if node == start_node else 0 for node in ig.vs['name']]
    communities = igraph.VertexClustering.FromAttribute(ig, 'cluster')
    palette = igraph.PrecalculatedPalette(['black', 'white'])

    out_fig_name = f'{output_dir}/diffusionP1Movie-{curr_time}.png'

    visual_style = {
        'vertex_color': ig.vs['color'],
        'vertex_label': ig.vs['label'],
        'vertex_label_dist': [3 for _ in G.values()],
        'edge_width': 5 * np.abs(ig.es['weight']),
        'mark_groups': True,
        'palette': palette,
        'layout': coords
    }

    fig, ax = plt.subplots(figsize=(20, 20))

    igraph.plot(communities, target=ax, **visual_style)

    ax.title.set_text(f'Diffuse {p1:.2f} from {start_node} at recursion level {recursion_level}')
    visited_patch = mpatches.Patch(color='red', label='Visited')
    unvisited_patch = mpatches.Patch(color='blue', label='Unvisited')
    ax.legend(handles=[visited_patch, unvisited_patch])

    plt.savefig(out_fig_name, dpi=100)


def connect_to_ext(adj_mat, start_node, visited_nodes):

    """
    Connect a node to its unvisited "extended" neighbors

    Parameters
    ----------
    adj_mat : The adjacency matrix that encodes the edge weights for the network.
    start_node : The node most recently visited by the network walker, from which p1 gets dispersed.
    visited_nodes : The history of previous draws in the node ranking sequence.

    Returns
    -------
    adj_mat_after : The adjacency matrix where the start_node is now connected to its unvisited "extended" neighbors.
    An extended neighbor is the neighbor of a neighbor.

    """
    ind = np.nonzero(adj_mat[start_node].values)
    start_node_nbors = adj_mat.index[ind]

    start_node_unvisited_nbors = [node for node in start_node_nbors if node not in visited_nodes]
    vN = [node for node in visited_nodes if node != start_node]  # visited nodes excluding start node

    ext_connections = {}  # or None

    if vN and not start_node_unvisited_nbors:
        keep = adj_mat.index.drop(vN)
        adj_mat_after = adj_mat.loc[keep, keep]

        conn_yes = adj_mat.index[ind]
        conn_no = adj_mat.index.drop(conn_yes)
        conn_no = conn_no.drop([node for node in conn_no if node in vN or node == start_node])

        for n1 in conn_yes:
            tmp = adj_mat.loc[n1, conn_no]
            ind = np.nonzero(tmp.values)
            for n2 in tmp.index[ind]:
                w = adj_mat.loc[n1, n2]
                ext_connections[n2] = w

        if ext_connections:
            adj_mat_after.loc[start_node, list(ext_connections.keys())] = list(ext_connections.values())
            adj_mat_after.loc[list(ext_connections.keys()), start_node] = list(ext_connections.values())

    else:
        adj_mat_after = adj_mat

    return adj_mat_after


def diffuse_p1(p1, start_node, G, visited_nodes, threshold_diff, adj_mat, verbose=None, out_dir='', r_level=1,
               coords=None):

    """
    Diffuse Probability P1 from a starting node

    Recursively diffuse probability from a starting node based on the connectivity of the network, representing the
    likelihood that a variable is most influenced by a perturbation in the starting node.


    Parameters
    ----------
    p1 : The probability being dispersed from the starting node, start_node, which is preferentially distributed between
     network nodes by the probability diffusion algorithm based solely on network connectivity.
    start_node : "Start node", or the node most recently visited by the network walker, from which p1 gets dispersed.
    G : A dictionary of probabilities, with keys being the node names in the network.
    visited_nodes : "Visited nodes", or the history of previous draws in the node ranking sequence.
    threshold_diff : When the probability diffusion algorithm exchanges this amount (threshold_diff) or less between
    nodes, the algorithm returns up the call stack.
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    verbose : If debugging or tracking a diffusion event, verbose=True will activate print statements. Default is False.
    out_dir : If specified, a image sequence will generate in the output directory specified.
    r_level : "Recursion level", or the current depth in the call stack caused by a recursive algorithm. Only relevant
    if out_dir is specified.
    coords : The x and y coordinates for each node in the network, to remain static between images. Only relevant if
    out_dir is specified.

    Returns
    -------
    G : A dictionary of returned probabilities after the diffusion of probability has truncated, with the keys being the
     node names in the network.

    Examples
    --------

    """

    n_tabs = '\t' * (r_level - 1)  # for verbose statements

    if verbose == 2:
        logging.debug(f'{n_tabs}prob. to diffuse:{p1} start_node: {start_node}, visitedNodes: {visited_nodes}')

    if out_dir:
        diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                            visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

    adj_mat2 = connect_to_ext(adj_mat=adj_mat, start_node=start_node, visited_nodes=visited_nodes)

    ind = np.nonzero(adj_mat2[start_node].values)
    start_node_nbors = adj_mat2.index[ind]
    start_node_unvisited_nbors = [node for node in start_node_nbors if node not in visited_nodes]

    if start_node_unvisited_nbors:
        w_edges_sum = sum(adj_mat2.loc[start_node_unvisited_nbors, start_node])

        for z in range(len(start_node_unvisited_nbors)):
            i_prob = p1 * abs(adj_mat2.loc[start_node_unvisited_nbors[z], start_node]) / w_edges_sum  # inherited prob
            G[start_node_unvisited_nbors[z]] = G[start_node_unvisited_nbors[z]] + i_prob

            if verbose == 2:
                logging.debug(f'{n_tabs}child#{z} {start_node_unvisited_nbors[z]} got {i_prob}')

            if out_dir:
                diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                    visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

            unv_ind = np.nonzero(adj_mat2[start_node_unvisited_nbors[z]].values)
            unv = adj_mat2.index[unv_ind]
            n_nbors = [x for x in unv if x in G]

            if n_nbors and i_prob / 2 > threshold_diff and len(visited_nodes) + 1 < len(G):
                G[start_node_unvisited_nbors[z]] = G[start_node_unvisited_nbors[z]] - i_prob / 2
                if verbose == 2:
                    logging.debug(f'{n_tabs}took {i_prob / 2} from child#{z}:{start_node_unvisited_nbors[z]} to send')

                if out_dir:
                    diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                        visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

                G = diffuse_p1(p1=i_prob / 2, start_node=start_node_unvisited_nbors[z], G=G,
                               visited_nodes=visited_nodes + [start_node_unvisited_nbors[z]],
                               threshold_diff=threshold_diff, adj_mat=adj_mat, verbose=verbose, out_dir=out_dir,
                               coords=coords, r_level=r_level + 1)

        if out_dir:
            diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)
    else:

        if verbose == 2:
            logging.debug(f'{start_node} is singleton or stranded by visited n.bors')
            logging.debug('Diffuse p1 uniformly amongst all unvisited nodes.')

        for node, val in G.items():
            if node not in visited_nodes:
                G[node] = val + p1 / (len(G) - len(visited_nodes))

        if out_dir:
            diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

    return G


def single_node_get_node_ranks(n, G, p1, threshold_diff, adj_mat, S=None, num_misses=None, verbose=None, out_dir='',
                               use_labels=False, coords=None):

    """
    Generate single-node node rankings ("fixed" walk)

    This function calculates the node rankings starting from a given perturbed variable in a subset of variables in
    the network.

    Parameters
    ----------
    n : The name of node ranking you want to calculate.
    G : A dictionary of probabilities with keys being the node names of the network.
    p1 : The probability that is preferentially distributed between network nodes by the probability diffusion algorithm
    based solely on network connectivity. The remaining probability (i.e., "p0") is uniformly distributed between
    network nodes, regardless of connectivity.
    threshold_diff : When the probability diffusion algorithm exchanges this amount or less between nodes, the algorithm
     returns up the call stack.
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    S : A list of node names in the subset you want the network walker to find.
    num_misses : The number of "misses" the network walker will tolerate before switching to fixed length codes for
    remaining nodes to be found.
    verbose : If True, print statements will execute as progress is made. Default is False.
    out_dir : If specified, a image sequence will generate in the output directory specified.
    use_labels : If True, node names will display next to their respective modes in the network. If False, node names
    will not display. Only relevant if out_dir is specified.
    coords : The x and y coordinates for each node in the network, to remain static between images.

    Returns
    -------
    curr_ns - A list of node names in the order they were drawn by the probability diffusion algorithm.

    Examples
    --------

    """

    p0 = 1 - p1
    if not S and (num_misses or out_dir != ''):
        print("You must also supply S if out_dir or num_misses is supplied.")
        return 0

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Node ranking {} of {}.".format(list(G).index(n) + 1, len(G)))

    stop_iterating = False
    start_node = n
    curr_gph = G
    count_misses = 0
    curr_ns = [start_node]  # current node set

    if out_dir:
        net_walk_snap_shot(adj_mat=adj_mat, G=G, visited_nodes=curr_ns, S=S, output_dir=out_dir, coords=coords,
                           img_num=len(curr_ns), use_labels=use_labels)

    while not stop_iterating:
        curr_gph = {node: 0 for node in curr_gph.keys()}  # clear probabilities
        base_p = p0 / (len(curr_gph) - len(curr_ns))
        # set unvisited nodes to base_p
        curr_gph.update({node: base_p for node in curr_gph.keys() if node not in curr_ns})

        curr_gph = diffuse_p1(start_node=start_node, G=curr_gph, visited_nodes=curr_ns, p1=p1, adj_mat=adj_mat,
                              threshold_diff=threshold_diff, verbose=verbose)

        # Sanity check - p1_event should add up to roughly 1
        p1_event = sum(curr_gph.values())
        if abs(p1_event - 1) > threshold_diff:
            extra_prob_to_diffuse = 1 - p1_event
            curr_gph.update({node: 0 for node in curr_gph.keys() if node in curr_ns})
            ind = [node for node in curr_gph.keys() if node not in curr_ns]
            curr_gph.update({node: val + extra_prob_to_diffuse / len(ind) for node, val in curr_gph.items()
                             if node in ind})

        # Set start_node to the node with the max probability in the new curr_gph
        max_prob = max(curr_gph, key=curr_gph.get)

        # Break ties: When there are ties, choose the first of the winners.
        start_node = max_prob

        if out_dir:
            net_walk_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, visited_nodes=curr_ns, S=S, coords=coords,
                               img_num=len(curr_ns), use_labels=use_labels)

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
            net_walk_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, visited_nodes=curr_ns, S=S, coords=coords,
                               img_num=len(curr_ns), use_labels=use_labels)

    return curr_ns, n
