import matplotlib.pyplot as plt
from numpy as np

from keras.layers import Input, Dense, Reshape, Flatten, Dropout
from keras.layers import BatchNormalization, Activation, ZeroPadding2D
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import UpSampling2D, Conv2D
from keras.models import Sequential, Model
from keras.optimizers import Adam

from graphs import Chimera


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

