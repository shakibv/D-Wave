import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.cm import get_cmap

from utilities import grouper


def draw_chimera(dwave):
    """
    Creates an image of the chimera graph.

    TODO: generalize assignment of node positions for more than chimeras in a row
    TODO: allow modification of drawing parameters
    TODO: draw nodes with colors of biases

    Parameters
    ------------
    dwave (dict): the input for D-Wave.
    """
    # Create graph
    graph = nx.Graph()
    graph.add_weighted_edges_from([
        (n1, n2, w) for (n1, n2), w in dwave.items()
    ])

    # Extract nodes and edges
    nodes = list(graph.nodes())
    edges = nx.get_edge_attributes(graph, 'weight')

    # Calculate number of chimera units and total number of nodes
    n_chimeras = (max(nodes) // 8) + 1
    n_nodes = n_chimeras * 8
    graph.add_nodes_from(range(n_nodes))

    # Assign positions to nodes
    positions = dict()
    for n_unit, unit in enumerate(grouper(range(n_nodes), 8)):
        for n_node, node in enumerate(unit):
            if n_node < (8 * n_unit + 4):
                positions[node] = (4 * n_unit + n_node, 4 * n_unit + 1.5)
            else:
                positions[node] = (4 * n_unit + 1.5, n_node - 8 * n_unit - 4)

    # Draw the graph
    fig = plt.figure()
    ax = fig.add_subplot()
    drawn_nodes = nx.draw(
        graph,
        pos=positions,
        with_labels=True,
    )
    drawn_edges = nx.draw_networkx_edges(
        graph,
        pos=positions,
        edgelist = edges.keys(),
        edge_color=list(edges.values()),
        width=3,
        edge_cmap=get_cmap('seismic'),
        edge_vmin=-1,
        edge_vmax=1,
    )
    plt.colorbar(drawn_edges)
    plt.show()


def main():
    # Sample D-Wave input
    qubit_biases = {
        (0, 0): 0.3333,
        (1, 1): -0.333,
        (4, 4): -0.333,
        (5, 5): 0.333
    }
    coupler_strengths = {
        (0, 4): 0.667,
        (0, 5): -1,
        (1, 4): 0.667,
        (1, 5): 0.667
    }
    Q = dict(qubit_biases)
    Q.update(coupler_strengths)

    draw_chimera(Q)

if __name__ == '__main__':
    main()
