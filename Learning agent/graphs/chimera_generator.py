import numpy as np

from itertools import combinations


def assign_edges(r=0.5, sampler=np.random.random):
    """Creates edges between nodes and assigns weights.

    TODO: expand assignment for more than a single unit
    TODO: get a better way to assign edges

    Parameters
    ------------
    r: (float) rate at which the edges are turned on.
    sampler: (function) sampling function.

    Returns
    ---------
    (dict) keys are tuples of the two edges connected and the values are the weight assignments.
    """
    edges = (
        (0, 4), (0, 5), (0, 6), (0, 7),
        (1, 4), (1, 5), (1, 6), (1, 7),
        (2, 4), (2, 5), (2, 6), (2, 7),
        (3, 4), (3, 5), (3, 6), (3, 7),
    )

    unit_cell = dict()
    for edge in edges:
        unit_cell[edge] = (r > np.random.random()) * sampler()

    return unit_cell


def assign_biases(sampler=np.random.random_sample):
    """Assigns biases to each node.
    """
    weights = dict()
    for node in range(0, 8):
        weights[(node, node)] = sampler()

    return weights

