import numpy as np

from copy import deepcopy
from itertools import combinations


class Chimera:
    """Creates a chimera graph for the D-Wave 2000Q.

    The chimera graph is based on a 16x16 grid of bipartite cells with
    8 qubits each.
    """
    def __init__(self):
        self.graph = self.create_graph()

    def assign_biases(self, sampler=np.random.random_sample):
        """Assigns biases to each node.
        """
        weights = dict()
        for node in range(0, 8):
            weights[(node, node)] = sampler()

        return weights

    def assign_edges(self, rows=1, columns=1, r=0.5, sampler=np.random.random):
        """Creates edges between nodes and assigns weights.

        TODO: expand assignment for more than a single unit
        TODO: get a better way to assign edges

        Parameters
        ------------
        r: (float) rate at which the edges are turned on.
        sampler: (function) sampling function.

        Returns
        ---------
        (dict) keys are tuples of the two edges connected and the values
        are the weight assignments.
        """
        def assign_edge(n1, n2, row, column):
            n1 += (128 * row) + (8 * column)
            n2 += (128 * row) + (8 * column)

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
            edge_sampler=np.random.random
    ):
        graph = self.assign_edges(
            rows=rows,
            columns=columns,
            r=r,
            sampler=edge_sampler
        )
        graph.update(self.assign_biases(bias_sampler))

        return graph

    def replicate(self, n_rows=1, n_columns=1, graph=None):
        """Replicates the cell to increase the size of the graph.

        Note that this method assumes the smallest size is a single cell
        using nodes 0-8.
        Note that this method assumes a 16x16 grid, the size of the
        D-Wave 2000Q.
        """
        if graph is None:
            graph = deepcopy(self.graph)

        replicated = dict()
        for row in range(n_rows):
            for column in range(n_columns):
                for (n1, n2), w in graph.items():
                    # Rescale nodes to correct position
                    n1 += (128 * row) + (8 * column)
                    n2 += (128 * row) + (8 * column)

                    # Ignore if one of the nodes is out of range
                    if 0 <= n1 < 2048 and 0 <= n2 < 2048:
                        edge = tuple(sorted([n1, n2]))
                        replicated[edge] = w

        return replicated

    def translate(self, x=0, y=0, graph=None):
        """Translates a graph by a number of cells."""
        if graph is None:
            graph = deepcopy(self.graph)

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
    from chimera_visualizer import draw_chimera
    from path_finding import fully_connected

    chimera = Chimera()
    while True:
        graph = chimera.create_graph(rows=1, columns=3)
        if fully_connected(graph):
            draw_chimera(graph)
            break


if __name__ == '__main__':
    main()
