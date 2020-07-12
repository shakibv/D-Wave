#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import pickle

from collections import deque
#from progress.bar import Bar

from dwave.system.samplers import DWaveSampler
from chimera import Chimera
import chimera_visualizer as vs
from dwave.cloud import Client
import send_to_SA as sim


def results(max_size, max_rows, max_columns):
    """
    Generates a .txt file with the description of a problem building block, TTS from D-Wave and SA for increasing problem size.
    Parameters
    ------------
    max_size: (int) the maximum size of a problem class
    max_rows: (int) the number of rows of unit cells in the building block
    max_columns: (int) the number of columns of unit cells in the building block
        
    Returns
    ---------
    
    """
    chimera = Chimera()
    Dsampler = DWaveSampler()
    B = chimera.create_graph(
                             rows=max_rows,
                             columns=max_columns,
                             r=np.random.power(17),
                             bias_sampler=sampler,
                             edge_sampler=sampler,
                             connected=True,
                             )
    P = []
    TTSD = []
    TTSS = []
    for i in range(max_size):
        C = chimera.replicate(B,i+1,1,connector_sampler=sampler)
        # input format for D-Wave solver
        h = C[0].copy()
        J = C[1].copy()
        # input format for SA solver
        size = 8*(i+1)
        f = open("bin/instance.txt","w+")
        f.write("%d\n" %size)
        for (n1,n2) in C[1].keys():
            f.write("%d %d %f\n" %(n1, n2, C[1][(n1,n2)]))
            if i==0:
                P.append(C[1][(n1,n2)])
            if (n1,n2) not in Dsampler.edgelist:
                J.pop((n1,n2))
        for n1 in C[0].keys():
            f.write("%d %d %f\n" %(n1, n1, C[0][n1]))
            if i==0:
                P.append(C[0][n1])
            if n1 not in Dsampler.nodelist:
                h.pop(n1)
        f.close()
        TTSS.append(sim.SA(r"/Users/archie/Google Drive/D-Wave-master/test_codes/bin", "instance.txt", size))
        with Client.from_config() as client:  # doctest: +SKIP
            solver = client.get_solver()
            if solver.check_problem(h,J):
                TTSD.append(TTS_DWave(h, J, num_reads=1000))
    f = open("test_data.txt","a+")
    f.write("Problem: [")
    for item in P:
        f.write("%f, " % item)
    f.write("]\n")
    f.write("TTSs from DWave: [")
    for item in TTSD:
        f.write("%f, " % item)
    f.write("]\n")
    f.write("TTSs from SA: [")
    for item in TTSS:
        f.write("%f, " % item)
    f.write("]\n")
    f.close()
#return B[1], TTSD, TTSS

def TTS_DWave(h, J, num_reads=1000):
    Dsampler = DWaveSampler()
    response = Dsampler.sample_ising(h, J, num_reads=1000)
    mine = float('inf')
    freq = 0
    for (sample, energy, num) in response.data():
        if energy<mine:
            mine = energy
            freq = num
    p_s = freq/num_reads
    p_d = 0.99
    if p_s == 1:
        T = 0
    else:
        T = 20*np.log2(1-p_d)/np.log2(1-p_s) #Using annealing time = 20 microsec, can put this as function parameter for generalization
    return T

def create_data(n_samples=1, max_rows=1, max_columns=1):
    """
    Generates a .txt file with the description of a problem building block, TTS from D-Wave and SA for increasing problem size, each for "n_samples" number of problem classes.
    Parameters
    ------------
    n_samples: (int) the number of different problem classes
    max_rows: (int) the number of rows of unit cells in the building block
    max_columns: (int) the number of columns of unit cells in the building block
        
    Returns
    ---------
    """
    for i in range(n_samples):
        R = results(2, 1, 1)
#        print(R[1])
#        print(R[2])


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


