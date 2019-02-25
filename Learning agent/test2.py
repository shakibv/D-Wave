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
import sys

from collections import deque
from itertools import groupby
from keras.layers import Input, Dense, Reshape, Flatten, Dropout
from keras.layers import BatchNormalization, Activation, ZeroPadding2D
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import UpSampling2D, Conv2D
from keras.models import Sequential, Model
from keras.optimizers import Adam
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering

from graphs import assign_edges, assign_biases, draw_chimera, fully_connected


class GAN:
    def __init__(
            self,
            n_rows,
            n_columns,
            n_channels,
            n_latent,
            optimizer=Adam(0.0002, 0.5),
            loss_function='binary_crossentropy',
    ):
        self.n_rows = n_rows
        self.n_columns = n_columns
        self.n_channels = n_channels
        self.shape = (self.n_rows, self.n_columns, self.n_channels)
        self.n_latent = n_latent

        self.optimizer = optimizer
        self.loss_function = loss_function

        # Build and compile the discriminator
        self.discriminator = self.build_discriminator()
        self.discriminator.compile(
            loss=self.loss_function,
            optimizer=self.optimizer,
            metrics=['accuracy'],
        )

        # Build the generator
        self.generator = self.build_generator()

        # The generator takes noise as input and generates matrices
        z = Input(shape=(self.n_latent,))
        matrix = self.generator(z)

        # For the combined model we will only train the generator
        self.discriminator.trainable = False

        # The discriminator takes generated images as input and determines validity
        validity = self.discriminator(matrix)

        # The combined model  (stacked generator and discriminator)
        # Trains the generator to fool the discriminator
        self.combined = Model(z, validity)
        self.combined.compile(
            loss=self.loss_function,
            optimizer=self.optimizer,
        )

    def build_generator(self):
        """The generator takes a vector of random numbers as input and a
        matrix as output.

        From the vector of random numbers, the model has four layers,
        increasing the size of the output until it is reshaped into the
        desired matrix shape, which is a sample within the set.
        """
        inputs = Input(shape=(self.n_latent,))

        model = Sequential()
        model.add(Dense(256, input_dim=self.n_latent))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(1024))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(np.prod(self.shape), activation='tanh'))
        model.add(Reshape(self.shape))
        model.summary()

        outputs = model(inputs)

        return Model(inputs, outputs)

    def build_discriminator(self):
        """The discriminator takes the matrix as input and a single
        number as output.

        From the flattened matrix, the model has three layers, reducing
        the output size until it is a single number, which determines
        whether or not the input is within the set.
        """
        inputs = Input(shape=self.shape)

        model = Sequential()
        model.add(Flatten(input_shape=self.shape))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(256))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(1, activation='sigmoid'))
        model.summary()

        outputs = model(inputs)

        return Model(inputs, outputs)

    def train(self, x_train, epochs, batch_size=128, sample_interval=50):
        # Adversarial ground truths
        valid = np.ones((batch_size, 1))
        fake = np.zeros((batch_size, 1))

        for epoch in range(epochs):

            # ---------------------
            #  Train Discriminator
            # ---------------------

            # Select a random batch of matrices
            indices = np.random.randint(0, x_train.shape[0], batch_size)
            matrices = x_train[indices]

            noise = np.random.normal(0, 1, (batch_size, self.n_latent))

            # Generate a batch of new images
            generated = self.generator.predict(noise)

            # Train the discriminator
            d_loss_real = self.discriminator.train_on_batch(matrices, valid)
            d_loss_fake = self.discriminator.train_on_batch(generated, fake)
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)

            # ---------------------
            #  Train Generator
            # ---------------------

            # Why is noise generated from normal distribution?
            noise = np.random.normal(0, 1, (batch_size, self.n_latent))

            # Train the generator (to have the discriminator label samples as valid)
            g_loss = self.combined.train_on_batch(noise, valid)

            # Plot the progress
            print("%d [D loss: %f, acc.: %.2f%%] [G loss: %f]" % (epoch, d_loss[0], 100*d_loss[1], g_loss))

            # If at save interval => sample graphs to check progress
            if epoch % sample_interval == 0:
                self.sample(epoch)

    def sample(self, epoch):
        r, c = 5, 5
        noise = np.random.normal(0, 1, (r * c, self.n_latent))
        generated = self.generator.predict(noise)

        # Rescale and plot


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


def hamming_distance(x, y):
    """The Hamming distance between two equal-length vectors."""
    return sum(e1 != e2 for e1, e2 in zip(x, y))


def cosine_similarity(x, y):
    """The cosine similarity between two equal-length vectors that
    measures the angle between them.
    """
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


def train_clustering(cluster, kwargs, x, graphs, scores):
    # Using k-means from scikit-learn with euclidean distance
    model = cluster(**kwargs)
    y = model.fit_predict(x)

    zipper = sorted(zip(x, y, graphs, scores), key=lambda x: x[1])
    x, y, graphs, scores = list(zip(*zipper))
    x = np.array(x)
    indx = 0
    for component, components in groupby(y):
        components = len(list(components))
        plt.plot(
            sorted(scores[indx:indx+components]),
            label='{}'.format(component)
        )
        # for i in range(5):
        #     print(scores[indx+i])
        #     draw_chimera(graphs[indx+i])
        indx += components

    plt.legend(loc='best')
    plt.show()

    return model


def main():
    # graphs, scores = create_data()
    graphs, scores = load_data()
    matrices = format_graphs(graphs)
    x = np.array([np.append(matrix, [score]) for matrix, score in zip(matrices, scores)])
    # x = np.array(matrices)

    # Show groups of scores
    # plt.hist(scores)
    # plt.show()

    cluster = train_clustering(
        cluster=KMeans,
        kwargs={'n_clusters': 12},
        x=x,
        graphs=graphs,
        scores=scores,
    )

    with open('./test_data/model2.dat', 'wb') as file:
        pickle.dump(cluster, file)

    x_train = np.array([
        np.reshape(vector[:-1], (8, 8, 1)) for vector in x if vector[-1] > 22
    ])

    gan = GAN(8, 8, 1, 100)
    gan.train(x_train, epochs=30000, batch_size=32, sample_interval=200)

    with open('./test_data/gan2.dat', 'wb') as file:
        pickle.dump(gan, file)


if __name__ == '__main__':
    main()

