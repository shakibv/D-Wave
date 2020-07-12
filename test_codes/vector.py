import yaml


edges = {
    (0, 4), (0, 5), (0, 6), (0, 7),
    (1, 4), (1, 5), (1, 6), (1, 7),
    (2, 4), (2, 5), (2, 6), (2, 7),
    (3, 4), (3, 5), (3, 6), (3, 7),
    (4, 12), (5, 13), (6, 14), (7, 15),
    (0, 128), (1, 129), (2, 130), (3, 131),
    (0, 0), (1, 1), (2, 2), (3, 3),
    (4, 4), (5, 5), (6, 6), (7, 7),
}

all_edges = set()
for (n1, n2) in edges:
    for i in range(16):

        for j in range(16):
            n1p = (128 * i) + (8 * j) + n1
            n2p = (128 * i) + (8 * j) + n2

            if j == 15:
                minum = 128 * (i + 1)
                maxum = minum + 8
                if minum <= n1p < maxum or minum <= n2p < maxum:
                    continue

            if 0 <= n1p < 2048 and 0 <= n2p < 2048:
                edge = tuple(sorted([n1p, n2p]))
                all_edges.add(edge)


indices = dict()
for indx, edge in enumerate(sorted(all_edges)):
    indices[edge] = indx


with open('./vector_indices.yaml', 'w') as file:
    yaml.dump(indices, file)

with open('./vector_indices.yaml', 'r') as file:
    indices = yaml.load(file)

print(indices)
