import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.cm import get_cmap


def draw_chimera(DWave_n, DWave_e, legend=1):
    """
        Creates an image of the chimera graph.
        
        Parameters
        ------------
        Dwave_n (dict): a dictonary for node biases.
        Dwave_e (dict): a dictonary for edges.
        """
    # Create graph
    graph = nx.Graph()
    #graph.add_nodes_from(DWave_n)
    graph.add_weighted_edges_from([
                                   (n1, n2, w) for (n1, n2), w in DWave_e.items() if w != 0.0
                                   ])
        
    # Extract nodes and edges
    nodes = list(graph.nodes())
    edges = nx.get_edge_attributes(graph, 'weight')
                                   
    # Calculate number of chimera units and total number of nodes
    n_chimeras = (max(nodes) // 8) + 1
    n_nodes = n_chimeras * 8
                                   
    # Assign positions to nodes
    positions = dict()
    for node in range(n_nodes):
        if node not in nodes:
            continue
        if (node // 4) % 2 == 0:
            x = 4 * ((node % 128) // 8) + ((node % 128) % 8)
            y = -4 * (node // 128)
            positions[node] = (x, y)
        else:
            x = (4 * ((node % 128) // 8)) + 1.5
            y = (-4 * (node // 128)) + ((node % 8) - 5.5)
            positions[node] = (x, y)

    # Draw the graph
    #fig = plt.figure()
    #ax = fig.add_subplot(dpi=300)
    drawn_nodes = nx.draw_networkx(
                               graph,
                               pos=positions,
                               node_size=200,
                               node_color=list(DWave_n.values()),
                               cmap=get_cmap('seismic',16),
                               vmin=-1,
                               vmax=1,
                               with_labels=False,
                               )
    drawn_edges = nx.draw_networkx_edges(
                                         graph,
                                         pos=positions,
                                         edgelist=edges.keys(),
                                         edge_color=list(edges.values()),
                                         width=3,
                                         edge_cmap=get_cmap('seismic',16),
                                         edge_vmin=-1,
                                         edge_vmax=1,
                                         )
    if legend==1:
        plt.colorbar(drawn_edges)
    plt.axis('off')
    plt.show()
#,cbar_kw=dict(ticks=np.arange(-3, 4), format=fmt

def main():
    # Sample D-Wave input
    qubit_biases = {
        0: 1,
        1: -1,
        4: -.75,
        5: .75
    }
    coupler_strengths = {
        (0, 4): 1,
        (0, 5): -1,
        (1, 4): -.75,
        (1, 5): 0.75
    }
    Q = dict(qubit_biases)
    Q.update(coupler_strengths)

    draw_chimera(qubit_biases, coupler_strengths)

if __name__ == '__main__':
    main()

