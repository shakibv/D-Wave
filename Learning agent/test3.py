"""
1. Create training data
    - Create graphs of size 1x1
    - Evaluate using SA solver
    - Store training data

2. Run clustering algorithm to cluster graphs that do best with SA solver

3. Run GAN to generate graphs and evaluate performance

4. Repeat with alterations:
    - Create graphs with varying size (maybe up to 3x3)
    - Evaluate on D-Wave solver as well and try using raw inputs and
      maximal difference
    - Use different methods of clustering
        - scipy hierarchial with dendrogram
    - Use restricted Boltzmann machines instead of GAN
"""


import matplotlib.pyplot as plt
import numpy as np
import pickle

from collections import deque
from itertools import groupby
from keras.layers import Input, Dense, Reshape, Flatten, Dropout
from keras.layers import BatchNormalization, Activation, ZeroPadding2D
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import UpSampling2D, Conv2D
from keras.models import Sequential, Model
from keras.optimizers import Adam
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering

from graphs import Chimera
from Performance_metrics import run_instances


DWAVE_TOKEN = 'ARL-e59e6801f04fb470e7e13f10756e97fc93720a63'


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
        n_samples = 5
        noise = np.random.normal(0, 1, (n_samples, self.n_latent))
        generated = self.generator.predict(noise)
        chimera = Chimera()

        # Rescale and plot
        for matrix in generated:
            graph = chimera.matrix_to_graph(matrix[0])
            score = run_instances(
                instances=graph,
                settings={
                    'dwave': False,
                    'sa': True,
                    'dwave_params': {},
                    'sa_params': {
                        '-s': 100,
                        '-r': 200,
                    },
                },
                verbose=True,
            )
            print(score)



def create_data():
    """Creates test data.

    The graphs are created with varying number of edges determined by
    the value `r` and with 1000 samples with each r.
    """
    graphs = deque()
    scores = deque()

    chimera = Chimera()
    for r in range(1, 10):
        r /= 10.0

        n = 0
        while n < 1000:
            graphs.append(
                chimera.create_graph(
                    rows=1,
                    columns=1,
                    r=r,
                    bias_sampler=sampler,
                    edge_sampler=sampler,
                    connected=True,
                )
            )

            n += 1

    scores = run_instances(
        instances=graphs,
        settings={
            'dwave': False,
            'sa': True,
            'dwave_params': {},
            'sa_params': {
                '-s': 100,
                '-r': 200,
            },
        },
        verbose=True,
    )

    with open('./test_data/test3.dat', 'wb') as file:
        pickle.dump((graphs, scores['sa']), file)

    return graphs, scores


def sampler():
    """Samples between -1 and +1."""
    return (np.random.random() * 2.0) - 1.0


def load_data():
    """Loads the training data."""
    with open('./test_data/test3.dat', 'rb') as file:
        graphs, scores = pickle.load(file)

    return graphs, scores


def train_clustering(kwargs, x, graphs, scores):
    # Using k-means from scikit-learn with euclidean distance
    model = KMeans(**kwargs)
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
        indx += components

    plt.legend(loc='best')
    plt.show()

    return model


def main():
    # create_data()
    graphs, scores = load_data()

    chimera = Chimera()
    matrices = [chimera.graph_to_matrix(graph) for graph in graphs]

    x = np.array([
        np.append(matrix, [score]) for matrix, score in zip(matrices, scores)
    ])

    cluster = train_clustering(
        kwargs={'n_clusters': 12},
        x=x,
        graphs=graphs,
        scores=scores,
    )

    cluster_number = int(input('Which cluster to generate?'))
    x_train = np.array([
        vector.reshape((1, vector.shape[0], 1)) for vector in x if cluster.predict([vector])[0] == cluster_number
    ])

    gan = GAN(1, x[0].shape[0], 1, 100)
    gan.train(x_train, epochs=10000, batch_size=64, sample_interval=500)
    gan.combined.save('./test_data/generator.h5')
    gan.discriminator.save('./test_data/discriminator.h5')


if __name__ == '__main__':
    main()

