import numpy as np

from boltzmann_machines.rbm import BernoulliRBM
from boltzmann_machines import DBM


def train_rbm1(x):
    rbm = BernoulliRBM(
        n_visible=6016,
        n_hidden=1024,
        W_init=0.01,
        vb_init=0,
        hb_init=0,
        n_gibbs_steps=1,
        learning_rate=0.05,
        momentum=np.geomspace(0.5, 0.9, 8),
        max_epoch=2,
        batch_size=48,
        l2=1e-3,
        sample_v_states=True,
        sample_h_states=True,
        dropout=None,
        sparsity_target=0.1,
        sparsity_cost=1e-5,
        sparsity_damping=0.9,
        dbm_first=True,
        metrics_config={
            'msre': True,
            'pl1': True,
            'train_metrics_every_iter': 500,
        },
        verbose=True,
        display_filters=30,
        display_hidden_activations=24,
        v_shape=(1, 6016,),
        random_seed=1337,
        dtype='float32',
        tf_saver_params={'max_to_keep': 1},
        model_path='./test_data/rbm1/'
    )
    rbm.fit(x)

    return rbm


def train_rbm2(x):
    rbm = BernoulliRBM(
        n_visible=1024,
        n_hidden=2048,
        W_init=0.05,
        vb_init=0,
        hb_init=0,
        n_gibbs_steps=np.repeat([1.0], 20),
        learning_rate=np.repeat([0.05], 20),
        momentum=np.geomspace(0.5, 0.9, 8),
        max_epoch=20,
        batch_size=48,
        l2=2e-4,
        sample_v_states=True,
        sample_h_states=True,
        dropout=None,
        sparsity_target=0.1,
        sparsity_cost=1e-5,
        sparsity_damping=0.9,
        dbm_last=True,
        metrics_config={
            'msre': True,
            'pl1': True,
            'train_metrics_every_iter': 500,
        },
        verbose=True,
        display_filters=0,
        display_hidden_activations=24,
        v_shape=(1, 1024,),
        random_seed=1111,
        dtype='float32',
        tf_saver_params={'max_to_keep': 1},
        model_path='./test_data/rbm2/'
    )
    rbm.fit(x)

    return rbm


def train_dbm(x_training, x_validation):
    x = np.concatenate((x_training, x_validation))

    rbm1 = train_rbm1(x)
    q = rbm1.transform(x)

    rbm2 = train_rbm2(q)
    g = rbm2.transform(q)

    dbm = DBM(
        rbms=(rbm1, rbm2),
        n_particles=100,
        v_particle_init=x_training[:100].copy(),
        h_particles_init=(q[:100].copy(), g[:100].copy()),
        n_gibbs_steps=1,
        max_mf_updates=50,
        mf_tol=1e-7,
        learning_rate=np.geomspace(1e-7, 5e-6, 400),
        momentum=np.geomspace(0.5, 0.9, 10),
        max_epoch=500,
        batch_size=100,
        l2=1e-7,
        max_norm=6.0,
        sample_v_states=True,
        sample_h_states=(True, True),
        sparsity_target=0.1,
        sparsity_cost=5e-5,
        sparsity_damping=0.9,
        train_metrics_every_iter=400,
        val_metrics_every_epoch=2,
        random_seed=2222,
        verbose=True,
        display_filters=10,
        display_particles=20,
        v_shape=(1, 6016),
        dtype='float32',
        tf_saver_params={'max_to_keep': 1},
        model_path='./test_data/dbm/',
    )
    dbm.fit(x_training, x_validation)

    return dbm
