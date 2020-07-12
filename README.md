**Todo**
* [Learning agent](#learning-agent)
* [Performance metrics](#performance-metrics)

Learning agent
--------------
Graph clustering: K-means
Graph generation: Generative adversarial networks


Performance metrics
-------------------
Ising oracle: Belief propagation algorithm, weighted MAX-2SAT solver or QUBO solver, time to target instead of time to solution, planted solution
Optimal time to solution: ?
Performance computer: ...   

## test_codes
-------------
Changed the following files: chimera_visualizer.py, chimera.py and create_training_data.py  
`chimera_visualizer.py`: now draws the problem graph in Chimera structure with coloured edges and nodes.
`chimera.py`: can now generate problems of different block sizes and increase the problem size by connecting blocks with random weights in either horizontal or vertical direction
`create_training_data.py`: it generates a text file with the description of the problem along with TTS of solving it using D-Wave and SA on increasing its size
