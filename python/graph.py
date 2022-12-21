import igraph
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
import logging

logging.getLogger().setLevel(logging.DEBUG)


def net_walk_snap_shot(adj_mat, G, output_dir, visited_nodes, S, coords, img_num=1, use_labels=True):

    """
    Capture the current location of a network walker.
    A network walker steps towards the node that inherited the highest probability from the last node that it
    stepped into.

    Parameters
    ----------
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    G : A list of probabilities, with names of the list being the node names in the network.
    output_dir : The local directory at which you want still PNG images to be saved.
    visited_nodes : A character vector of node names, storing the history of previous draws in the node ranking.
    S : A character vector of node names in the subset you want the network walker to find.
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
    # adj mat has to be DataFrame currently
    # TODO: potentially change method so it works with matrix and add vertex_names as argument

    ig.vs['color'] = ['white' if ig.vs['name'][i] not in S else 'green' for i in range(len(G))]

    if use_labels:
        ig.vs['label'] = ig.vs['name']
    else:
        ig.vs['label'] = ['' for i in range(len(ig.vs['name']))]

    ig.vs['cluster'] = [1 if node in visited_nodes else 0 for node in ig.vs['name']]
    communities = igraph.VertexClustering.FromAttribute(ig, 'cluster')
    palette = igraph.PrecalculatedPalette(['dark red', 'white'])

    out_fig_name = f'{output_dir}/netWalkMovie{S.index(visited_nodes[0])}_{img_num}.png'

    visual_style = {
        'vertex_color' : ig.vs['color'],
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
    G : A list of probabilities, with names of the list being the node names in the network.
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
    # adj mat has to be DataFrame currently
    # TODO: potentially change method so it works with matrix and add vertex_names as argument

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
        'vertex_label_dist': [3 for v in G.values()],
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
    adj_matAfter : The adjacency matrix where the start_node is now connected to its unvisited "extended" neighbors.
    An extended neighbor is the neighbor of a neighbor.

    """

    #start_node_nbors = list(adj_mat[abs(adj_mat.loc[start_node, :]) > 0].index)
    ind = np.nonzero(adj_mat[start_node].values)
    start_node_nbors = adj_mat.index[ind]

    start_node_unvisited_nbors = [node for node in start_node_nbors if node not in visited_nodes]
    vN = [node for node in visited_nodes if node != start_node]  # visited nodes excluding start node

    ext_connections = {}  # or None

    if vN and not start_node_unvisited_nbors:
        adj_mat_after = adj_mat.drop(labels=vN).drop(columns=vN)
        connections = adj_mat[start_node]
        conn_yes = connections[abs(connections) > 0].to_dict()
        conn_no = connections[connections == 0]
        conn_no = conn_no.drop([node for node in conn_no.index if node in vN or node == start_node]).to_dict()

        if conn_no:
            for n1 in conn_no:
                if conn_yes:
                    for n2 in conn_yes:
                        if abs(adj_mat.loc[n2, n1]) > 0:
                            conn_no[n1] = adj_mat.loc[n2, n1]
                            ext_connections[n1] = conn_no[n1]

        if ext_connections:
            adj_mat_after.loc[start_node, list(ext_connections.keys())] = list(ext_connections.values())
            adj_mat_after.loc[list(ext_connections.keys()), start_node] = list(ext_connections.values())

    else:
        adj_mat_after = adj_mat

    return adj_mat_after


def diffuse_p1(p1, start_node, G, visited_nodes, threshold_diff, adj_mat, verbose=False, out_dir='', r_level=1,
               coords=None):

    """
    Diffuse Probability P1 from a starting node

    Recursively diffuse probability from a starting node based on the
    connectivity of the network, representing the likelihood that a
    variable is most influenced by a perturbation in the starting node.


    Parameters
    ----------
    p1 : The probability being dispersed from the starting node, start_node, which is preferentially distributed between network
     nodes by the probability diffusion algorithm based solely on network connectivity.
    start_node : "Start node", or the node most recently visited by the network walker, from which p1 gets dispersed.
    G : A list of probabilities, with names of the list being the node names in the network.
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

    if verbose:
        logging.debug(f'{n_tabs}prob. to diffuse:{p1} start_node: {start_node}, visitedNodes: {visited_nodes}')

    if out_dir:
        diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                            visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

    adj_mat2 = connect_to_ext(adj_mat=adj_mat, start_node=start_node, visited_nodes=visited_nodes)
    #start_node_nbors = list(adj_mat2[abs(adj_mat2.loc[start_node, :]) > 0].index)
    ind = np.nonzero(adj_mat2[start_node].values)
    start_node_nbors = adj_mat2.index[ind]
    start_node_unvisited_nbors = [node for node in start_node_nbors if node not in visited_nodes]

    if start_node_unvisited_nbors:
        w_edges_sum = sum(adj_mat2.loc[start_node_unvisited_nbors, start_node])

        for z in range(len(start_node_unvisited_nbors)):
            i_prob = p1 * abs(adj_mat2.loc[start_node_unvisited_nbors[z], start_node]) / w_edges_sum  # inherited prob
            G[start_node_unvisited_nbors[z]] = G[start_node_unvisited_nbors[z]] + i_prob

            if verbose:
                logging.debug(f'{n_tabs}child#{z} {start_node_unvisited_nbors[z]} got {i_prob}')

            if out_dir:
                diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                    visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

            #unv = list(adj_mat2[abs(adj_mat2.loc[:, start_node_unvisited_nbors[z]]) > 0].index)
            #n_nbors = [node for node in G.keys() if node in unv]
            unv_ind = np.nonzero(adj_mat2[start_node_unvisited_nbors[z]].values)
            unv = adj_mat2.index[unv_ind]
            n_nbors = [x for x in unv if x in G]

            if n_nbors and i_prob / 2 > threshold_diff and len(visited_nodes) + 1 < len(G):
                G[start_node_unvisited_nbors[z]] = G[start_node_unvisited_nbors[z]] - i_prob / 2
                if verbose:
                    logging.debug(f'{n_tabs}took {i_prob / 2} from child#{z}:{start_node_unvisited_nbors[z]} to send')

                if out_dir:
                    diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                        visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

                G = diffuse_p1(p1=i_prob / 2, start_node=start_node_unvisited_nbors[z], G=G, visited_nodes=visited_nodes + [start_node_unvisited_nbors[z]],
                               threshold_diff=threshold_diff, adj_mat=adj_mat, verbose=verbose, out_dir=out_dir,
                               coords=coords, r_level=r_level + 1)

        if out_dir:
            diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)
    else:

        if verbose:
            logging.debug(f'{start_node} is singleton or stranded by visited n.bors')
            logging.debug('Diffuse p1 uniformly amongst all unvisited nodes.')

        for node, val in G.items():
            if node not in visited_nodes:
                G[node] = val + p1 / (len(G) - len(visited_nodes))

        #G.update({node: val + p1 / (len(G) - len(visited_nodes)) for node, val in G.items() if node not in visited_nodes})

        if out_dir:
            diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

    return G


def diffuse_p1_iterative(p1, G, start_node, visited_nodes, adj_mat, threshold_diff=0.01, verbose=False, alpha=0.5):
    """
    Diffuse Probability P1 from a starting node

    Iteratively diffuse probability from a starting node based on the connectivity of the network, representing the
    likelihood that a variable is most influenced by a perturbation in the starting node.

    Parameters
    ----------
    p1 : The probability being dispersed from the starting node, start_node, which is preferentially distributed between network
    nodes by the probability diffusion algorithm based solely on network connectivity.
    start_node : "Start node", or the node most recently visited by the network walker, from which p1 gets dispersed.
    G : A list of probabilities, with names of the list being the node names in the network.
    visited_nodes : "Visited nodes", or the history of previous draws in the node ranking sequence.
    threshold_diff : When the probability diffusion algorithm exchanges this amount (threshold_diff) or less between
    nodes, the algorithm returns up the call stack.
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    verbose : If debugging or tracking a diffusion event, verbose=True will activate print statements. Default is False.
    alpha :

    Returns
    -------
    G : A dictionary of returned probabilities after the diffusion of probability has truncated, with the keys being the
    node names in the network.

    Examples
    --------
    """

    queue = [(start_node, p1, visited_nodes)]
    path_end = False

    while queue:

        current_node, current_probability, visited_nodes = queue.pop(0)

        lvl = len(visited_nodes) - 1

        if verbose and path_end:
            n_tabs = '\t' * (lvl - 1)
            logging.debug(f'{n_tabs}child {current_node} got {2 * current_probability}')
            logging.debug(f'{n_tabs}took {current_probability} from child {current_node} to send')

        n_tabs = '\t' * lvl

        if verbose:
            logging.debug(
                f'{n_tabs}prob. to diffuse:{current_probability} start node: {current_node}, visited nodes: {visited_nodes}')

        path_end = False

        adj_mat2 = connect_to_ext(adj_mat=adj_mat, start_node=current_node, visited_nodes=visited_nodes)

        unvisited_nodes = G.keys() - visited_nodes
        start_node_nbors = list(adj_mat2[abs(adj_mat2.loc[current_node, :]) > 0].index)
        start_node_unvisited_nbors = [node for node in start_node_nbors if node not in visited_nodes]

        if not start_node_unvisited_nbors:

            path_end = True

            if verbose:
                logging.debug(f'{current_node} is singleton or stranded by visited n.bors')
                logging.debug('Diffuse p1 uniformly amongst all unvisited nodes.')

            if not unvisited_nodes:
                continue

            to_add = current_probability / len(unvisited_nodes)
            for node in unvisited_nodes:
                G[node] += to_add

            continue

        sum_weights = sum(adj_mat2.loc[list(start_node_unvisited_nbors), current_node])

        k = 0

        for node in start_node_unvisited_nbors:

            inherited_prob = current_probability * adj_mat2.loc[node, current_node] / sum_weights
            G[node] += inherited_prob

            if verbose and not k:
                try:
                    logging.debug(f'{n_tabs}child {node} got {inherited_prob}')
                except IndexError:
                    pass

            #unv = list(adj_mat2[abs(adj_mat2.loc[:, node]) > 0].index)
            unv = adj_mat2[adj_mat2[node] != 0].index

            n_nbors = list(set(unv) & set(G.keys()))

            if n_nbors and inherited_prob * alpha > threshold_diff and len(visited_nodes) + 1 < len(G):
                G[node] -= inherited_prob * alpha
                queue.insert(k, (node, inherited_prob * alpha, visited_nodes + [node]))

            k += 1

        if verbose:
            try:
                logging.debug(f'{n_tabs}took {queue[0][1]} from child {queue[0][0]} to send')
            except IndexError:
                pass

    return G


def diffuse_p1_iterative_non_aligned(p1, G, start_node, visited_nodes, adj_mat, threshold_diff=0.01, alpha=0.5):
    """
    Diffuse Probability P1 from a starting node

    Iteratively diffuse probability from a starting node based on the connectivity of the network, representing the
    likelihood that a variable is most influenced by a perturbation in the starting node.

    Parameters
    ----------
    p1 : The probability being dispersed from the starting node, start_node, which is preferentially distributed between network
    nodes by the probability diffusion algorithm based solely on network connectivity.
    start_node : "Start node", or the node most recently visited by the network walker, from which p1 gets dispersed.
    G : A list of probabilities, with names of the list being the node names in the network.
    visited_nodes : "Visited nodes", or the history of previous draws in the node ranking sequence.
    threshold_diff : When the probability diffusion algorithm exchanges this amount (threshold_diff) or less between
    nodes, the algorithm returns up the call stack.
    adj_mat : The adjacency matrix that encodes the edge weights for the network, G.
    alpha :

    Returns
    -------
    G : A dictionary of returned probabilities after the diffusion of probability has truncated, with the keys being the
    node names in the network.

    Examples
    --------
    """

    queue = [(start_node, p1, set(visited_nodes))]

    while queue:

        current_node, current_probability, visited_nodes = queue.pop(0)

        unvisited_nodes = G.keys() - visited_nodes
        start_node_unvisited_nodes = set(filter(lambda x: adj_mat.loc[current_node, x] > 0, unvisited_nodes))

        if not len(start_node_unvisited_nodes):
            if not len(unvisited_nodes):
                continue
            to_add = current_probability / len(unvisited_nodes)
            for node in unvisited_nodes:
                G[node] += to_add
            continue

        sum_weights = 0
        for node in start_node_unvisited_nodes:
            sum_weights += adj_mat.loc[current_node, node]
        for node in start_node_unvisited_nodes:
            inherited_prob = current_probability * (adj_mat.loc[current_node, node] / sum_weights)
            G[node] += inherited_prob
            if inherited_prob * alpha > threshold_diff and len(
                    set(filter(lambda x: adj_mat.loc[node, x] > 0, unvisited_nodes))) > 0:
                G[node] -= inherited_prob * alpha
                queue.append((node, inherited_prob * alpha, visited_nodes.union({node})))

    return G
