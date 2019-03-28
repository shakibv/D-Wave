import numpy as np
import yaml

from copy import deepcopy
from itertools import combinations
from os import path

from .chimera_visualizer import draw_chimera
from .path_finding import fully_connected


class Chimera:
    """A class for creating and managing chimera graphs.

    The chimera graph is based on a 16x16 grid of bipartite cells with
    8 qubits each. The basic functionality of the class is outlined
    below.

    1. To create a graph, use the self.create_graph() method.

    >>> chimera = Chimera()
    >>> graph = chimera.create_graph()

    2. To increase the size of the graph, use the self.replicate()
       method.

    >>> graphs = list()
    >>> for i in range(1, 16):
    >>>     for j in range(1, 16):
    >>>         graphs.append(
    >>>             chimera.replicate(
    >>>                 graph=graph,
    >>>                 rows=i,
    >>>                 columns=j,
    >>>             )
    >>>         )

    3. The graphs can be made to occupy different regions of the D-Wave
       architecture by using the self.translate() method. For example, to
       translate the graph 3 chimera units down and 2 chimera units
       right:

    >>> translated_graph = chimera.translate(
    >>>     graph=graph,
    >>>     x=3,
    >>>     y=-2,
    >>> )

    """
    def __init__(self):
        filename = path.dirname(path.abspath(__file__)) + '/vector_indices.yaml'
        with open(filename, 'r') as file:
            self.indices = yaml.load(file)

    def assign_biases(self, graph, sampler=np.random.random_sample):
        """Assigns biases to nodes in a given graph.

        The intention is for the edges to be assigned first, and the
        biases assigned afterword.

        Parameters
        ------------
        graph  : (dict)     the output of self.assign_edges().
        sampler: (function) the sampling function to assign biases to
                            each node.

        Returns
        ---------
        (dict) the updated graph with the assigned biases.
        """
        biases = dict()
        for (n1, n2), w in graph.items():
            if biases.get((n1, n1)) is None:
                biases[(n1, n1)] = sampler()

            if biases.get((n2, n2)) is None:
                biases[(n2, n2)] = sampler()

        graph.update(biases)

        return graph

    def assign_edges(self, rows=1, columns=1, r=0.5, sampler=np.random.random):
        """Creates edges between nodes and assigns weights.

        Parameters
        ------------
        rows   : (int)      the number of rows of chimera units in the
                            graph.
        columns: (int)      the number of columns of chimera units in the
                            graph.
        r      : (float)    rate at which the edges are turned on.
        sampler: (function) the sampling function used to assign weights
                            to each edge.

        Returns
        ---------
        (dict) keys are tuples of the two edges connected and the values
        are the weight assignments.
        """
        def assign_edge(n1, n2, row, column):
            n1 += (128 * row) + (8 * column)
            n2 += (128 * row) + (8 * column)

            if column == 15:
                minum = 128 * (i + 1)
                maxum = minum + 8
                if minum <= n1 < maxum or minum <= n2 < maxum:
                    return None

            if 0 <= n1 < 2048 and 0 <= n2 < 2048:
                edge = tuple(sorted([n1, n2]))
                return edge

            return None

        internal_edges = (
            (0, 4), (0, 5), (0, 6), (0, 7),
            (1, 4), (1, 5), (1, 6), (1, 7),
            (2, 4), (2, 5), (2, 6), (2, 7),
            (3, 4), (3, 5), (3, 6), (3, 7),
        )

        horizontal_edges = (
            (4, 12), (5, 13), (6, 14), (7, 15),
        )

        vertical_edges = (
            (0, 128), (1, 129), (2, 130), (3, 131),
        )

        graph = dict()
        for row in range(rows):
            for column in range(columns):
                # Assign internal edges
                for n1, n2 in internal_edges:
                    edge = assign_edge(n1, n2, row, column)
                    if edge is not None:
                        graph[edge] = (r > np.random.random()) * sampler()

                # Ensure at least one horizontal and vertical edge
                indx = np.random.randint(4)
                edge = assign_edge(*horizontal_edges[indx], row, column)
                if edge is not None:
                    graph[edge] = sampler()

                indx = np.random.randint(4)
                edge = assign_edge(*vertical_edges[indx], row, column)
                if edge is not None:
                    graph[edge] = sampler()

                # Assign the rest of the external edges
                for n1, n2 in horizontal_edges:
                    edge = assign_edge(n1, n2, row, column)
                    if edge is not None and graph.get(edge) is None:
                        graph[edge] = (r > np.random.random()) * sampler()

                for n1, n2 in vertical_edges:
                    edge = assign_edge(n1, n2, row, column)
                    if edge is not None and graph.get(edge) is None:
                        graph[edge] = (r > np.random.random()) * sampler()

        return graph

    def create_graph(
            self,
            rows=1,
            columns=1,
            r=0.5,
            bias_sampler=np.random.random,
            edge_sampler=np.random.random,
            connected=True,
    ):
        """Creates a chimera graph in the D-Wave input format.

        The edges are assigned first, and the biases are assigned after.
        Sampler functions should be created to assign the weights and
        biases for the edges and nodes. These functions should take no
        arguments as input and return a single number.

        Parameters
        ------------
        rows        : (int)      the number of rows of chimera units in
                                 the graph.
        columns     : (int)      the number of columns of chimera units
                                 in the graph.
        r           : (float)    the rate at which the edges are turned
                                 on.
        bias_sampler: (function) the sampler to assign biases to each
                                 node.
        edge_sampler: (function) the sampler to assign weights to each
                                 edge.
        connected   : (bool)     determines if the returned graph should
                                 be fully connected.

        Returns
        ---------
        (dict) a graph with the assigned edges and biases.
        """
        # Create the initial graph
        graph = self.assign_edges(
            rows=rows,
            columns=columns,
            r=r,
            sampler=edge_sampler,
        )

        while connected:
            # Create a new graph if it is not connected
            if fully_connected(graph):
                break
            else:
                graph = self.assign_edges(
                    rows=rows,
                    columns=columns,
                    r=r,
                    sampler=edge_sampler,
                )

        graph = self.assign_biases(
            graph=graph,
            sampler=bias_sampler
        )

        return graph

    @staticmethod
    def find_n_nodes(graph):
        """Determines the maximum numbered node in a given graph.

        Parameters
        ------------
        graph: (dict) the D-Wave input or the output from self.create_graph().

        Returns
        ---------
        (int) the maximum numbered graph.
        """
        nodes = set()
        for (n1, n2) in graph.keys():
            nodes.update({n1, n2})

        return max(nodes)

    def find_shape(self, graph):
        """Determines the shape of the graph in terms of the number of
        chimera units.

        Parameters
        ------------
        graph: (dict) the D-Wave input or the output from
                      self.create_graph().

        Returns
        ---------
        (int, int) the number of rows and columns of chimera units in the
        graph.
        """
        n_nodes = self.find_n_nodes(graph)
        n_rows = (n_nodes // 128)
        n_columns = ((n_nodes % 128) // 8) + 1

        return n_rows, n_columns

    def graph_to_vector(self, graph):
        """Formats the D-Wave inputs into the corresponding adjacency
        matrix.

        The graph is represented as a vector with length 6016 (the upper
        bound on the possible number of edges on the D-Wave
        architecture).

        Parameters
        ------------
        graph  : (dict) the D-Wave input or the output from
                        self.create_graph().

        Returns
        ---------
        (np.array) a 6016 length vector.
        """
        vector = np.zeros(shape=6016)
        for (n1, n2), w in graph.items():
            indx = self.indices[tuple(sorted([n1, n2]))]
            vector[indx] = w

        return vector

    @staticmethod
    def size(graph):
        """Determines the size of the problem.

        Note that this doesn't subtract the nodes that extend out of the
        unit cell. This means that the number returned is actually higher
        than the actual size of the problem. This is a bug, not a
        feature.
        """
        nodes = set()
        for (n1, n2) in graph.keys():
            nodes.update({n1, n2})

        return len(nodes)

    def vector_to_graph(self, vector):
        """Formats the graph vector into the D-Wave input format.

        The elements of the graph vector have a direct correspondance to
        weights or biases on the D-Wave input.

        Parameters
        ------------
        matrix: (np.array) an nxn array.

        Returns
        ---------
        (dict) the graph formatted for D-Wave input.
        """
        keys = self.indices.key()

        graph = dict()
        for indx, w in enumerate(vector):
            graph[keys[indx]] = w

        return graph

    @staticmethod
    def plot(graph):
        """Plots the graph."""
        draw_chimera(graph)

    def replicate(self, graph, rows=1, columns=1):
        """Replicates the cell to increase the size of the graph.

        Parameters
        ------------
        graph  : (dict) the D-Wave input or output of
                        self.create_graph().
        rows   : (int)  the number of times the graph is replicated
                        downward, or the number of rows of graph
                        replications.
        columns: (int)  the number of times the graph is replicated
                        rightward, of the number of columns of graph
                        replications.

        Returns
        ---------
        (dict) the replicated graph in the D-Wave input format.
        """
        n_rows, n_columns = self.find_shape(graph)

        replicated = dict()
        for row in range(rows):
            for column in range(columns):
                for (n1, n2), w in graph.items():
                    # Rescale nodes to correct position
                    n1 += (128 * row * n_rows) + (8 * column * n_columns)
                    n2 += (128 * row * n_rows) + (8 * column * n_columns)

                    # Ignore if one of the nodes is out of range
                    if 0 <= n1 < 2048 and 0 <= n2 < 2048:
                        edge = tuple(sorted([n1, n2]))
                        replicated[edge] = w

        return replicated

    def translate(self, graph, x=0, y=0):
        """Translates a graph by a number of chimera units.

        The translation is based on the normal cartesian plane, where
        positive x is rightward and positive y is upward.

        Parameters
        ------------
        x    : (int)  the number of chimera units to translate the graph
                      in the horizontal direction. This can also be
                      negative.
        y    : (int)  the number of chimera units to translate the graph
                      in the vertical direction. This can also be
                      negative.
        graph: (dict) the D-Wave input or output of self.create_graph().
        """
        translated = dict()
        for (n1, n2), w in graph.items():
            # Rescale nodes to correct position
            n1 += (8 * x) - (128 * y)
            n2 += (8 * x) - (128 * y)

            # Ignore if one of the nodes is out of range
            if 0 <= n1 < 2048 and 0 <= n2 < 2048:
                edge = tuple(sorted([n1, n2]))
                translated[edge] = w

        return translated


def main():
    import matplotlib.pyplot as plt

    from scipy.linalg import svd

    from chimera_visualizer import draw_chimera

    chimera = Chimera()
    graph = chimera.create_graph(rows=2, columns=3)
    chimera.plot(graph)
    chimera.plot(chimera.replicate(graph, 2, 2))

    matrix = chimera.graph_to_matrix(graph)
    U, D, Vt = svd(matrix)
    D = np.diag(D)
    F = U @ D
    vector = F.shape[1]
    rows, columns = matrix.shape
    for i in range(rows):
        for j in range(columns):
            vector += matrix[i, j] * (F[i] + F[j])


if __name__ == '__main__':
    main()
