from igraph import Graph
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

    ig = Graph.Weighted_Adjacency(adj_mat, mode='max')
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

    ig = Graph.Weighted_Adjacency(adj_mat, mode='max')
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

    start_node_nbors = list(adj_mat[abs(adj_mat.loc[start_node, :]) > 0].index)
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
    start_node_nbors = list(adj_mat2[abs(adj_mat2.loc[start_node, :]) > 0].index)
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

            n_nbors = {node: val for node, val in G.items() if
                       node in adj_mat2[abs(adj_mat2.loc[:, start_node_unvisited_nbors[z]]) > 0].index}

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

        G.update({node: val + p1 / (len(G) - len(visited_nodes)) for node, val in G.items() if node not in visited_nodes})

        if out_dir:
            diffusion_snap_shot(adj_mat=adj_mat, G=G, output_dir=out_dir, p1=p1, start_node=start_node,
                                visited_nodes=visited_nodes, coords=coords, recursion_level=r_level)

    return G
