# Graphs

Contains the `Chimera` class, which handles:

- Random generation of chimera graphs of specified size (specified by the number
  of cells, not the number of nodes) by using `Chimera.create_graph()`
- The generated graphs can be checked to be a connected graph by specifying
  `connected=True` when creating the graph
- Replication of graphs to increase their size (limited to the size of the
  D-Wave architecture) by using `Chimera.replicate()`
- Translation of graphs to different regions of the D-Wave architecture
- Visualization of graphs by using `Chimera.translate()`


# Test Data

Previous scripts, models, and training data. Most of the files in this folder
are legacy and will be outdated.


# Utilities

Basic functions unrelated to graph generation or machine learning.
