"""
1. Make test data with no weights, only 1 and 0 as edge values.
2. Score graphs with certain shapes as consistently higher.
3. Use clustering algorithm with distance function as Hamming distance. Note
   that this is not perfect since it does not take into account possible
   interactions that are not through direct connections.
4. From the clustering, use a GAN to generate new samples in the cluster.
"""


import matplotlib.pyplot as plt
import numpy as np
import pickle

from collections import deque
from itertools import groupby
from scipy.cluster.hierarchy import fclusterdata
from sklearn.cluster import KMeans

from graphs import assign_edges, assign_biases, draw_chimera, fully_connected


def create_data():
    """Creates test data.

    The graphs are created with varying number of edges determined by
    the value `r` and with 1000 samples with each r.

    The edge weights are only 1 or 0.
    """
    graphs = deque()
    scores = deque()

    for r in range(1, 10):
        r /= 10.0

        n = 0
        while n < 1000:
            graph = assign_edges(r=r, sampler=sampler)
            graph.update(assign_biases(sampler=sampler))

            if fully_connected(graph):
                score = score_graph(graph)

                graphs.append(graph)
                scores.append(score_graph(graph))

                n += 1
                print(n, score)

    with open('./test_data/test2.dat', 'wb') as file:
        pickle.dump((graphs, scores), file)

    return graphs, scores


def sampler():
    """Distribution of the weights.

    In this case, it's just a 1 or 0.
    """
    return 1.0


def score_graph(graph):
    """Gives a score to the graph based on the connected edges.

    For this test, the following scheme is used to score the graph.
    - Graphs with edges connecting the rim nodes (0, 7, 3, 4) or (1, 6, 2, 5) are given 10 points
    - Each node that has three connections is given 5 points
    """
    score = 0

    outer_rim = ((0, 7), (3, 7), (3, 4), (0, 4))
    inner_rim = ((1, 6), (2, 6), (2, 5), (1, 5))

    if all(graph.get(edge, 0.0) == 1.0 for edge in outer_rim):
        score += 10 + np.random.random() - 0.5

    if all(graph.get(edge, 0.0) == 1.0 for edge in inner_rim):
        score += 10 + np.random.random() - 0.5

    for n1 in range(0, 8):
        edges = 0

        for n2 in range(0, 8):
            if n1 != n2:
                edges += graph.get((n1, n2), 0.0)

        if edges == 3.0:
            score += 5 + np.random.random() - 0.5

    return score


def load_data():
    """Loads the test data."""
    with open('./test_data/test2.dat', 'rb') as file:
        graphs, scores = pickle.load(file)

    return graphs, scores


def format_graphs(graphs):
    """Formats the D-Wave inputs into a 1-D vector."""
    matrices = deque()
    for graph in graphs:
        matrix = np.zeros(shape=(8, 8))
        for (n1, n2), w in graph.items():
            matrix[n1, n2] = w
            matrix[n2, n1] = w

        matrices.append(matrix.flatten())

    return matrices


def train_clustering():
    pass


def hamming_distance(x, y):
    """The Hamming distance between two equal-length vectors."""
    return sum(e1 != e2 for e1, e2 in zip(x, y))


def cosine_similarity(x, y):
    """The cosine similarity between two equal-length vectors that
    measures the angle between them.
    """
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


def train_gan():
    pass


def main():
    graphs, scores = create_data()
    graphs, scores = load_data()
    matrices = format_graphs(graphs)
    x = np.array([
        np.append(matrix, [score]) for matrix, score in zip(matrices, scores)
    ])

    # Show groups of scores
    # plt.hist(scores)
    # plt.show()

    # Using k-means from scikit-learn with euclidean distance
    model = KMeans(n_clusters=12)
    model.fit(x)
    y = model.predict(x)

    zipper = sorted(zip(x, y, graphs, scores), key=lambda x: x[1])
    x, y, graphs, scores = list(zip(*zipper))
    x = np.array(x)
    indx = 0
    for component, components in groupby(y):
        components = len(list(components))
        plt.plot(
            sorted(x[indx:indx+components, -1]),
            label='{}'.format(component)
        )
        for i in range(5):
            print(scores[indx+i])
            draw_chimera(graphs[indx+i])
        indx += components

    plt.legend(loc='best')
    plt.show()

    # Using agglomerative clustering from scipy (doesn't converge)
    # fcluster = fclusterdata(
    #     X=x,
    #     t=10,
    #     metric=cosine_similarity,
    # )
    # plt.hist(fcluster)
    # plt.show()


if __name__ == '__main__':
    main()
