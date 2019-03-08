def run_instance(instance={(0, 128): 0.9634326790980382, (13, 13): 0.05884838050338814, (6, 6): 0.8254353772059342, (1, 6): 0.0, (12, 12): 0.42888892824804614, (7, 7): 0.13954571203557697, (0, 7): 0.426274214580713, (5, 13): 0.4296085972603224, (2, 6): 0.0, (0, 5): 0.0, (1, 4): 0.18536386687975215, (128, 128): 0.07984593877937396, (2, 4): 0.8121249644712547, (3, 7): 0.47313140959380473, (2, 5): 0.5823642000984334, (14, 14): 0.46563741383959745, (6, 14): 0.0, (2, 7): 0.0, (7, 15): 0.0394906423857051, (129, 129): 0.7500484228206831, (131, 131): 0.6626937065196934, (1, 129): 0.1862002193279464, (5, 5): 0.669401067957829, (1, 7): 0.0, (4, 4): 0.2869378734769994, (3, 131): 0.9153825258929245, (0, 0): 0.2276063274962139, (15, 15): 0.17977118796812597, (0, 6): 0.44354427888710846, (2, 130): 0.0, (4, 12): 0.09861194280484942, (3, 3): 0.4274950012817952, (3, 6): 0.8480613521078071, (2, 2): 0.34538793954976255, (0, 4): 0.0, (1, 5): 0.13550825337241645, (3, 5): 0.0, (3, 4): 0.7656510529136723, (130, 130): 0.5454203562887191, (1, 1): 0.7294036921536221}, settings={'sa':True,'dwave':True}):
    """
    Given an input problem instance, run_instance runs and returns all of the 
    statistics from the various solvers.
    
    Parameters:
    -----------
    instance: a dictionary containing all of the information necessary to fully
              specify the problem instance
              
    Returns:
    --------
    result: 
    """
    # Import solvers
    from simulated_annealing.send_to_SA import SA
    from dwave_interface.Problem_Post import post_problem, convert_to_text
    
    # Initialize output
    result = {}
    
    # Run (only) selected solvers
    
    # Run on D-Wave
    if settings['dwave']:
        import base64
        dwave = post_problem(instance)
        #result['dwave'] = {'anneal_time_per_run':dwave['anneal_time_per_run']}
        print(dwave)
        print(base64.b64decode(dwave['answer']['solutions']))
        
    # Run Simulated Annealing
    if settings["sa"]:
        # Data first needs to be saved before run
        import os 
        dir_path = os.path.dirname(os.path.realpath(__file__))
        
        instance_file = open("temp_file.txt", "w")
        instance_file.write("This file is written automatically for the purposes "+
        "of evaluation by the Simulated Annealing algorithm and is safe to delete.\n")
        for key in instance:
            for value in key:
                instance_file.write(str(value)+" ")
            instance_file.write(str(instance[key])+"\n")
        instance_file.close()
        
        for root, dirs, files in os.walk(dir_path):
            if "an_ms_r1_fi.h" in files:
                 solver_path = os.path.join(root)
                 
        print(solver_path)
                
        sa = SA(dir_path, "temp_file.txt", verbose=True, solverDir=solver_path, solver="an_ss_ge_fi_vdeg")
        print(sa)
        
    return result
    
run_instance()