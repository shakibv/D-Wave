"""
Input
-------
We have training data consisting of chimera graphs G and their
corresponding time to solutions T.

P(G, T)

Output
--------
We want to sample from distribution D, created from clustering graphs G
and their corresponding time to solutions T, to generate graph g.

P(g | G, T)

Model
-------
We know P(T |  G), and from the model we want to generate P(G | T).

Try with no weights, just 1 and 0
Try with data where more geometric properties are higher targets
The distance between graphs should be taken into account for clustering, instead of cartesian, calculate hamming distance?
May have to take into account interactions that are not through direct connections
The cluster that we find interesting can be given to GAN to generate new samples
"""


import matplotlib.pyplot as plt
import numpy as np
import pickle

from itertools import groupby
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from generator import assign_edges, assign_biases
from utilities import fully_connected


def create_data():
    n = 0
    graphs = dict()
    while n < 10000:
        unit_cell = assign_edges(r=n/10000+0.1)
        unit_cell.update(assign_biases())
        if fully_connected(unit_cell):
            graphs[calculate_speed(unit_cell)] = unit_cell
            n += 1
            print(n)

    with open('./test_data/test1.dat', 'wb') as file:
        pickle.dump(graphs, file)


def calculate_speed(graph):
    # speed = sum(graph.values())
    speed = sum([abs(n1 - n2) * w for (n1, n2), w in graph.items()])

    return speed


def load_data():
    with open('./test_data/test.dat', 'rb') as file:
        data = pickle.load(file)

    x = list()
    for speed, graph in data.items():
        matrix = np.zeros(shape=(8, 8))
        for (n1, n2), w in graph.items():
            matrix[n1, n2] = w
            matrix[n2, n1] = w

        matrix = matrix.flatten()
        np.concatenate([matrix, [speed]])
        x.append(matrix)

    return x


def train_model(x):
    # Determine best number of components
    # Determined 53 was best number of components through AIC
    # n_components = range(50, 55)
    # models = [GaussianMixture(n_components=n) for n in n_components]
    # aics = [model.fit(x).aic(x) for model in models]
    # print(aics)

    # plt.plot(n_components, aics)
    # plt.show()

    model = GaussianMixture(n_components=53)
    model.fit(x)

    if not model.converged_:
        print('Model did not converge')
        exit()
    else:
        with open('./test_data/model.dat', 'wb') as file:
            pickle.dump(model, file)

    return model


def load_model():
    with open('./test_data/model.dat', 'rb') as file:
        model = pickle.load(file)

    return model


def main():
    # create_data()

    # Training data has shape [n_samples, n_features]
    x = np.array(load_data())

    # Reduces the dimensionality of training data
    # In this case, keeping 99% of components
    pca = PCA(n_components=0.99, whiten=True)
    x_train = pca.fit_transform(x)

    # model = train_model(x_train)
    model = load_model()
    y = model.predict(x_train)

    zipper = sorted(zip(x, y), key=lambda x: x[-1])
    x, y = list(zip(*zipper))
    x = np.array(x)
    indx = 0
    for component, components in groupby(y):
        components = len(list(components))
        plt.plot(
            sorted(x[indx:indx+components, -1]),
            label='{}'.format(component)
        )
        indx += components

    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    main()
