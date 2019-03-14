#import sys
#sys.path.append(r"C:\Users\nolan\OneDrive\Documents\GitHub\D-Wave")
#print(sys.path)

#from Learning_agent import graphs

#l = list()
#g = graphs.chimera.Chimera()
#a = g.create_graph()
#for i in range(1, 3):
#    for j in range(1, 3):
#        l.append(g.replicate(graph=a, rows=i, columns=j))
#g.plot(graph=l[len(l)-1])
#print(l)
#g.plot(graph={(2, 6): 0.22325651460425322, (13, 13): 0.4190657375908027, (6, 6): 0.6449733125529578, (1, 6): 0.0, (12, 12): 0.36984407182820034, (7, 7): 0.4626979957557906, (0, 7): 0.4553935911379925, (5, 13): 0.9573705612476301, (0, 128): 0.8506612996311672, (0, 5): 0.5057135639431842, (1, 4): 0.07159584524918838, (128, 128): 0.27277600871972496, (2, 4): 0.4688041355524043, (3, 7): 0.0, (2, 5): 0.0, (14, 14): 0.0681987939964529, (6, 14): 0.0, (2, 7): 0.0, (7, 15): 0.0, (129, 129): 0.8877803411227139, (131, 131): 0.18010719513083873, (1, 129): 0.0, (5, 5): 0.9789899616495775, (1, 7): 0.0, (4, 4): 0.24024811390313472, (3, 131): 0.9090634762663303, (0, 0): 0.6565189955075825, (15, 15): 0.7025385375526844, (0, 6): 0.9870773453000441, (2, 130): 0.0, (4, 12): 0.0, (3, 3): 0.2787077403951821, (3, 6): 0.0, (2, 2): 0.8330261710766741, (0, 4): 0.7087790759645399, (1, 5): 0.41186835495804364, (3, 5): 0.0, (3, 4): 0.7925899118249323, (130, 130): 0.38269735903973345, (1, 1): 0.18741263558643972})

# Functional relative imports
#from solvers.simulated_annealing.send_to_SA import SA
#from dwave_interface.Problem_Post import post_problem, convert_to_text

# Broken Relative Imports
#from Learning_agent.graphs.chimera as Chimera
#from Learning_agent.graphs.chimera_visualizer import draw_chimera
#from Learning_agent.graphs.path_finding import path_exists, fully_connectede

# Relative import of copies
#from temp.chimera_visualizer import draw_chimera
#from temp.chimera import Chimera
#from temp.path_finding import path_exists, fully_connected

#def run_metrics():
    # Assuming that all necessary information is in the variable data
 #   print(SA(r"C:\Users\nolan\OneDrive\Documents\GitHub\D-Wave\Performance metrics\solvers\simulated_annealing\bin", "instance.txt"))
  #  return None
    
#.run_metrics()

from Performance_metrics.run import run_instances

x = {(0, 128): 0.9634326790980382, (13, 13): 0.05884838050338814, 
    (6, 6): 0.8254353772059342, (1, 6): 0.0, (12, 12): 0.42888892824804614, 
    (7, 7): 0.13954571203557697, (0, 7): 0.426274214580713, 
    (5, 13): 0.4296085972603224, (2, 6): 0.0, (0, 5): 0.0, 
    (1, 4): 0.18536386687975215, (128, 128): 0.07984593877937396, 
    (2, 4): 0.8121249644712547, (3, 7): 0.47313140959380473, 
    (2, 5): 0.5823642000984334, (14, 14): 0.46563741383959745, 
    (6, 14): 0.0, (2, 7): 0.0, (7, 15): 0.0394906423857051, 
    (129, 129): 0.7500484228206831, (131, 131): 0.6626937065196934, 
    (1, 129): 0.1862002193279464, (5, 5): 0.669401067957829, (1, 7): 0.0, 
    (4, 4): 0.2869378734769994, (3, 131): 0.9153825258929245, 
    (0, 0): 0.2276063274962139, (15, 15): 0.17977118796812597, 
    (0, 6): 0.44354427888710846, (2, 130): 0.0, (4, 12): 0.09861194280484942, 
    (3, 3): 0.4274950012817952, (3, 6): 0.8480613521078071, 
    (2, 2): 0.34538793954976255, (0, 4): 0.0, (1, 5): 0.13550825337241645, 
    (3, 5): 0.0, (3, 4): 0.7656510529136723, (130, 130): 0.5454203562887191, 
    (1, 1): 0.7294036921536221}
    
print(run_instances(instances = x, settings = {"dwave" : True, "sa" : True,
    "dwave_params" : {"num_reads" : 100, "solver" : "DW_2000Q_VFYC_2_1"},
    "sa_params" : {"-r" : 100, "-s" : 200}}))