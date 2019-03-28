import numpy as np
import pickle

from collections import deque
from progress.bar import Bar

from graphs.chimera import Chimera
from Performance_metrics import run_instances


# This is Austin's D-Wave token
DWAVE_TOKEN = 'ARL-e59e6801f04fb470e7e13f10756e97fc93720a63'



def create_graphs(n=10000, max_rows=1, max_columns=1):
    graphs = deque()
    chimera = Chimera()

    bar = Bar('Processing', max=n*max_rows*max_columns)
    for rows in range(max_rows):
        for columns in range(max_columns):
            for n in range(n):
                graphs.append(
                    chimera.create_graph(
                        rows=rows+1,
                        columns=columns+1,
                        r=np.random.power(17),
                        bias_sampler=sampler,
                        edge_sampler=sampler,
                        connected=True,
                    )
                )

                bar.next()

    bar.finish()

    return graphs


def create_data(max_rows=1, max_columns=1, settings=None):
    """Creates test data.

    The graphs are created with varying number of edges determined by
    the value r and with 1000 samples with each r.

    Parameters
    ------------
    max_rows   : (int)  the maximum number of rows of chimera units to
                        create problem instances.
    max_columns: (int)  the maximum number of columns of chimera units to
                        create problem instances.
    settings   : (dict) the settings for the SA and D-Wave solvers. By
                        default, both solvers are used.

    """
    graphs = deque()
    scores = deque()

    chimera = Chimera()
    for r in range(1, 10):
        r /= 10.0

        for rows in range(max_rows):
            for columns in range(max_columns):
                n = 0
                while n < 1000:
                    graphs.append(
                        chimera.create_graph(
                            rows=rows,
                            columns=columns,
                            r=r,
                            bias_sampler=sampler,
                            edge_sampler=sampler,
                            connected=True,
                        )
                    )

                    n += 1

    if settings is None:
        settings = {
            'dwave': True,
            'sa': True,
            'dwave_params': {
                'token': DWAVE_TOKEN,
            },
            'sa_params': {
                '-s': 100,
                '-r': 200,
            }
        }

    scores = run_instances(
        instances=graphs,
        settings=settings,
        verbose=True,
    )

    return graphs, scores


def sampler():
    """Randomly chooses a value from the 16-bit possible values for
    weights and biases.

    There are in total 17 possible values, but since we are sampling
    based on non-zero values, there are only 16 in this sampler.
    """
    values = [
        -1.0, -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125,
        +0.125, +0.25, +0.375, +0.5, +0.625, +0.75, +0.875, +1.0,
    ]

    indx = np.random.randint(0, 15)

    return values[indx]


def save_data(filename, graphs, scores):
    """Writes the training data to disk.

    Parameters
    ------------
    filename: (str)  path to the file where the data is to be saved.
    graphs  : (iter) iterable of graphs. This should be the same length
                     as the scores and should be pickle-able.
    scores  : (iter) iterable of the scores or times. This should be the
                     same length as the graphs and should be pickle-able.

    """
    with open(filename, 'wb') as file:
        pickle.dump((graphs, scores['sa']), file)

    return graphs, scores


def load_data(filename):
    """Loads the training data from disk.

    Parameters
    ------------
    filename: (str) path to the file where the data is saved.

    Returns
    ---------
    (tuple) the un-pickled graphs and scores.

    """
    with open(filename, 'rb') as file:
        graphs, scores = pickle.load(file)

    return graphs, scores


def main():
    graphs = create_graphs(max_rows=3, max_columns=3)

    filename = '/home/augh/D-Wave/Learning agent/test_data/graphs.dat'
    with open(filename, 'wb') as file:
        pickle.dump(graphs, file)


if __name__ == '__main__':
    main()
