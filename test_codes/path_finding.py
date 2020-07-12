from collections import defaultdict, deque
from itertools import combinations


def path_exists(graph, start, end):
    """Checks if there is a path that exists between the start and end nodes.

    Taken from https://en.wikipedia.org/wiki/Breadth-first_search.

    Parameters
    ------------
    graph:    (dict) there should be a key for each node in the graph and its
              corresponding value should be an iterable of connections.
    start:    (int) the node to start from.
    end:      (int) the node to find a path to.

    Returns
    ---------
    (bool) True if a path between the start and end nodes exists, False
    otherwise.
    """
    queue = deque([[start]])
    visited = set()

    while queue:
        path = queue.popleft()
        vertex = path[-1]

        if vertex == end:
            return True
        elif vertex not in visited:
            for neighbour in graph.get(vertex, []):
                new_path = list(path)
                new_path.append(neighbour)
                queue.append(new_path)

            visited.add(vertex)

    return False


def fully_connected(dwave):
    """Checks if the chimera graph has all edges connected to each other.

    Not all nodes have to be connected to each other, but this simply ensures
    that every node with an existing edge is not isolated.

    Parameters
    ------------
    dwave:    (dict) input for the D-Wave.

    Returns
    ---------
    (bool) True if the chimera graph is fully connected, False otherwise.
    """
    graph = defaultdict(deque)
    for (n1, n2), weight in dwave.items():
        if weight != 0:
            graph[n1].append(n2)
            graph[n2].append(n1)

    for n1, n2 in combinations(graph.keys(), 2):
        if not path_exists(graph, n1, n2):
            return False
    else:
        return True

def main():
    qubit_biases = {(0, 0): 0.3333, (1, 1): -0.333, (4, 4): -0.333, (5, 5): 0.333}
    coupler_strengths = {(0, 4): 0.667, (0, 5): -1, (1, 4): 0.667, (1, 5): 0.667}
    Q = dict(qubit_biases)
    Q.update(coupler_strengths)

    print(fully_connected(Q) == True)


if __name__ == '__main__':
    main()
