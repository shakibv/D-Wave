"""
1. Given an nxn symmetric adjacency matrix, perform singular value
   decomposition (or another method of spectral decomposition) so

        A = USV*

   However, A should be a matrix representative of the entire dataset. Perhaps A should be some combinatorial formulation of all the graphs at their maximum size. This is to avoid the arbitrary location of small graphs.

2. Since A is a real symmetric matrix, the decomposition simplifies to

        A = USUT

3. The vectors representing every node is then computed by multiplying
   the diagonal matrix S with eigenvalues with their respective
   eigenvectors as the column vectors of Q to get

        F = US

4. Each graph is then represented by some combination of the node
   vectors. A simple example is the sum

        G = sum(wi + fi)

"""


import numpy as np

from scipy.linalg import svd


class GraphVector:
    def __init__(self):
        self.node_vectors = None

    def fit(self, data):
        U, D, Vt = svd(data)
        F = U @ D

        self.node_vectors = F

    def predict(self, graph):
        rows, columns = graph.shape

        vector = np.zeros(self.node_vectors.shape[1])
        for i in rows:
            for j in columns:
                vector += self.node_vectors[i] + self.node_vectors[j]

        return vector
