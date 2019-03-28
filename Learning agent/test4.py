import matplotlib.pyplot as plt
import numpy as np
import pickle

from collections import deque
from itertools import groupby
from sklearn.cluster import KMeans, AgglomerativeClustering

from generator.dbm import train_dbm
from graphs.chimera import Chimera
from graphs.create_training_data import create_graphs, create_data, save_data, load_data, DWAVE_TOKEN


def train_cluster(x, scores, algorithm=KMeans, plot=False, **kwargs):
    model = algorithm(**kwargs)
    y = model.fit_predict(x)

    if plot:
        print('Plotting clusters')
        plot_cluster(x, y, scores)

    return model, y


def plot_cluster(x, y, scores):
    zipper = sorted(zip(x, y, scores), key=lambda x: x[1])
    x, y, scores = list(zip(*zipper))

    indx = 0
    for component, components in groupby(y):
        components = len(list(components))
        plt.plot(
            sorted(scores[indx:indx+components]),
            label='{}'.format(component),
        )

        indx += components

    plt.legend(loc='best')
    plt.show()


def save_model(filename, model):
    with open(filename, 'wb') as file:
        pickle.dump(model, file)


def load_mode(filename):
    with open(filename, 'rb') as file:
        model = pickle.load(file)

    return model


def main():
    # CREATE TRAINING DATA
    # graphs = create_graphs(max_rows=3, max_columns=3)

    filename = '/home/augh/D-Wave/Learning agent/test_data/graphs.dat'
    # with open(filename, 'wb') as file:
    #     pickle.dump(graphs, file)

    print('Loading data')
    with open(filename, 'rb') as file:
        graphs = pickle.load(file)
    print('Loading finished')

    settings = {
        'dwave': False,
        'sa': True,
        'dwave_params': {
            'token': DWAVE_TOKEN,
        },
        'sa_params': {
            '-s': 100,
            '-r': 200,
        }
    }

    chimera = Chimera()
    for graph in graphs:
        sized_graphs = deque()
        for rows in range(4):
            for columns in range(4):
                sized_graphs.append(chimera.replicate(
                    graph=graph,
                    rows=rows+1,
                    columns=columns+1,
                ))

        scores = run_instances(
            instances=sized_graphs,
            settings=settings,
            verbose=True,
        )
        print(scores)
        exit()

    exit()

    # -------------------
    # Get training data
    # -------------------

    # graphs, scores = create_data(max_rows=4, max_columns=4)
    # save_data('./test_data/test4.dat', graphs, scores)
    graphs, scores = load_data('./test_data/test3.dat')

    chimera = Chimera()
    vectors = [chimera.graph_to_vector(graph) for graph in graphs]

    x = np.array([
        np.append(vector, [score]) for vector, score in zip(vectors, scores)
    ])

    print('Data loaded and transformed into vectors')

    # ---------
    # Cluster
    # ---------

    print('Performing K-means clustering')
    kmeans_model, y = train_cluster(
        x=x,
        scores=scores,
        algorithm=KMeans,
        plot=True,
        n_clusters=12,
    )

    print('Performing agglomerative clustering')
    # agglomerative_model = train_cluster(
    #     x=x,
    #     scores=scores,
    #     algorithm=AgglomerativeClustering,
    #     plot=True,
    #     n_clusters=12,
    #     affinity='euclidean',
    #     linkage='complete',
    # )

    # -------------------------------
    # Choose cluster and split data
    # -------------------------------
    chosen = int(input('Which cluster should I take?'))

    x_training = list()
    x_validation = list()
    for vector, n_cluster in zip(vectors, y):
        if n_cluster == chosen:
            if np.random.random() < 0.66:
                x_training.append(vector)
            else:
                x_validation.append(vector)

    x_training = np.array(x_training)
    x_validation = np.array(x_validation)

    # -----------
    # Generator
    # -----------

    # Try GAN
    # Try RBM or DBM
    print('Running deep Boltzmann machine')
    dbm = train_dbm(x_training, x_validation)


if __name__ == '__main__':
    main()

