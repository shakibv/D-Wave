﻿def run_instances(instances, settings={'sa':True,'dwave':True,'dwave_params':{'num_reads':100},'sa_params':{}}, verbose=False, plot_settings=None):
    """
    Given an input problem instance(s), run_instances runs and returns the time to
    solution (TTS) statistics from the various solvers.
    
    Parameters:
    -----------
    instances: a dictionary containing all of the information necessary to fully
              specify the problem instance OR a list of dictionaries containing
              all of the information necessary to fully specify several problem
              instances. Fully compatible with output from chimera.create_graph()
              method.
              
    settings: (optional, though highly recommended) a dictionary used to 
              configure the solvers (see below)
    
    verbose: boolean (optional) if True, more solution details are printed to the
             screen
              
    Returns:
    --------
    result: a dictionary with a key for each utilized solver ("dwave" and/or "sa") 
            paired with a list of TTS values for each instance in instances, 
            respecting the input order
            
    =====================
    Quick Start/Tutorial:
    =====================
    
    --------------------------
    Example 1 (Basic Example):
    --------------------------
    
    Suppose that you want to evaluate the time to solution for the following 
    instance both on D-Wave and with the simulated annealer:
        
    >>> x = {(0, 128): 0.9634326790980382, (13, 13): 0.05884838050338814, 
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
    
    The following function call will evaluate that instance with both solvers
    using the "DW_2000Q_VFYC_2_1" solver on D-Wave and the "an_ss_ge_fi_vdeg" 
    simulated annealing solver. Each function will sample 100 times, the 
    simulated annealer will perform 200 sweeps per sample, and fully detailed
    output of the solution will be displayed.
    
    >>> run_instances(instances = x, settings = {"dwave" : True, "sa" : True,
    "dwave_params" : {"num_reads" : 100, "solver" : "DW_2000Q_VFYC_2_1"},
    "sa_params" : {"-r" : 100, "-s" : 200}})
    
    Should output (exact values will vary):
        
    {'dwave': [268.922...], 'sa': [559.876...]}
    
    -------------------------
    Example 2 (Full Example):
    -------------------------
    
    Suppose that you want to evaluate the time to solution for the following 
    list of instances of the same problem, each with an increasing size on both
    D-Wave and with the simulated annealer:
        
    >>> y = [{(2, 6): 0.5160852615777896, (13, 13): 0.9319646681708734, 
    (6, 6): 0.8606380663808766, (1, 6): 0.0, (12, 12): 0.6513698529343085, 
    (7, 7): 0.16563375072880127, (0, 7): 0.0, (5, 13): 0.0, 
    (4, 12): 0.7338340760527127, (0, 5): 0.23660222273755804, 
    (1, 4): 0.7323245483844127, (128, 128): 0.045783443174777894, 
    (2, 4): 0.0, (3, 7): 0.9038168616058962, (2, 5): 0.0, 
    (14, 14): 0.487573139345321, (6, 14): 0.0, (2, 7): 0.9426353078362628, 
    (7, 15): 0.6455592536238158, (0, 128): 0.7406443781286747, 
    (129, 129): 0.7132101466390699, (131, 131): 0.9854602920895158, 
    (1, 5): 0.25902842868967213, (5, 5): 0.4570409226555594, 
    (1, 7): 0.24653882575442454, (4, 4): 0.6055189513596845, (3, 131): 0.0, 
    (0, 0): 0.16268264406083688, (15, 15): 0.009710795948105533, 
    (0, 6): 0.12062024006460492, (3, 3): 0.9640744750887139, 
    (3, 6): 0.15441807354502235, (2, 2): 0.5217999096926298, 
    (0, 4): 0.7663816443072197, (2, 130): 0.0, (3, 5): 0.0, (3, 4): 0.0, 
    (130, 130): 0.7916268753073616, (1, 129): 0.7862934888901292, 
    (1, 1): 0.316204105625695}, {(11, 139): 0.0, (8, 13): 0.23660222273755804, 
    (13, 13): 0.4570409226555594, (9, 9): 0.316204105625695, 
    (6, 6): 0.8606380663808766, (12, 12): 0.6055189513596845, (10, 12): 0.0, 
    (7, 7): 0.16563375072880127, (0, 7): 0.0, (137, 137): 0.7132101466390699, 
    (10, 15): 0.9426353078362628, (15, 23): 0.6455592536238158, 
    (4, 12): 0.7338340760527127, (1, 6): 0.0, (128, 128): 0.045783443174777894, 
    (11, 13): 0.0, (3, 7): 0.9038168616058962, (2, 5): 0.0, (14, 22): 0.0, 
    (13, 21): 0.0, (21, 21): 0.9319646681708734, (8, 14): 0.12062024006460492, 
    (0, 128): 0.7406443781286747, (12, 20): 0.7338340760527127, 
    (7, 15): 0.6455592536238158, (9, 15): 0.24653882575442454, 
    (4, 4): 0.6055189513596845, (3, 131): 0.0, (0, 0): 0.16268264406083688, 
    (8, 8): 0.16268264406083688, (20, 20): 0.6513698529343085, 
    (1, 5): 0.25902842868967213, (136, 136): 0.045783443174777894, 
    (3, 6): 0.15441807354502235, (14, 14): 0.8606380663808766, 
    (3, 3): 0.9640744750887139, (6, 14): 0.0, (129, 129): 0.7132101466390699, 
    (8, 15): 0.0, (1, 1): 0.316204105625695, (8, 136): 0.7406443781286747, 
    (11, 14): 0.15441807354502235, (2, 6): 0.5160852615777896, (9, 14): 0.0, 
    (22, 22): 0.487573139345321, (10, 13): 0.0, (2, 2): 0.5217999096926298, 
    (5, 5): 0.4570409226555594, (1, 4): 0.7323245483844127, (5, 13): 0.0, 
    (130, 130): 0.7916268753073616, (9, 137): 0.7862934888901292, (3, 5): 0.0, 
    (0, 5): 0.23660222273755804, (11, 15): 0.9038168616058962, 
    (0, 4): 0.7663816443072197, (8, 12): 0.7663816443072197, (11, 12): 0.0, 
    (131, 131): 0.9854602920895158, (2, 7): 0.9426353078362628, 
    (9, 13): 0.25902842868967213, (1, 129): 0.7862934888901292, 
    (10, 10): 0.5217999096926298, (11, 11): 0.9640744750887139, 
    (23, 23): 0.009710795948105533, (15, 15): 0.16563375072880127, 
    (0, 6): 0.12062024006460492, (138, 138): 0.7916268753073616, 
    (10, 14): 0.5160852615777896, (139, 139): 0.9854602920895158, 
    (1, 7): 0.24653882575442454, (2, 130): 0.0, (10, 138): 0.0, (3, 4): 0.0, 
    (2, 4): 0.0, (9, 12): 0.7323245483844127}, {(132, 132): 0.6055189513596845, 
    (13, 13): 0.9319646681708734, (6, 6): 0.8606380663808766, 
    (259, 259): 0.9854602920895158, (12, 12): 0.6513698529343085, 
    (7, 7): 0.16563375072880127, (0, 7): 0.0, (130, 134): 0.5160852615777896, 
    (1, 6): 0.0, (142, 142): 0.487573139345321, (128, 128): 0.16268264406083688, 
    (131, 135): 0.9038168616058962, (3, 7): 0.9038168616058962, (2, 5): 0.0, 
    (130, 133): 0.0, (6, 14): 0.0, (135, 135): 0.16563375072880127, 
    (143, 143): 0.009710795948105533, (0, 128): 0.7406443781286747, 
    (7, 15): 0.6455592536238158, (129, 133): 0.25902842868967213, 
    (128, 134): 0.12062024006460492, (4, 4): 0.6055189513596845, 
    (130, 135): 0.9426353078362628, (3, 131): 0.0, (0, 0): 0.16268264406083688, 
    (1, 5): 0.25902842868967213, (129, 135): 0.24653882575442454, (131, 259): 0.0, 
    (3, 6): 0.15441807354502235, (14, 14): 0.487573139345321, (3, 3): 0.9640744750887139, 
    (129, 257): 0.7862934888901292, (3, 5): 0.0, (134, 142): 0.0, (129, 129): 0.316204105625695, 
    (129, 134): 0.0, (1, 1): 0.316204105625695, (128, 132): 0.7663816443072197, 
    (2, 6): 0.5160852615777896, (131, 134): 0.15441807354502235, 
    (256, 256): 0.045783443174777894, (135, 143): 0.6455592536238158, 
    (2, 2): 0.5217999096926298, (128, 133): 0.23660222273755804, 
    (5, 5): 0.4570409226555594, (1, 4): 0.7323245483844127, (5, 13): 0.0, 
    (130, 130): 0.5217999096926298, (0, 5): 0.23660222273755804, (128, 135): 0.0, 
    (258, 258): 0.7916268753073616, (133, 133): 0.4570409226555594, 
    (0, 4): 0.7663816443072197, (131, 131): 0.9640744750887139, 
    (2, 7): 0.9426353078362628, (130, 258): 0.0, (1, 129): 0.7862934888901292, 
    (140, 140): 0.6513698529343085, (131, 132): 0.0, (132, 140): 0.7338340760527127, 
    (257, 257): 0.7132101466390699, (130, 132): 0.0, (15, 15): 0.009710795948105533, 
    (141, 141): 0.9319646681708734, (0, 6): 0.12062024006460492, 
    (128, 256): 0.7406443781286747, (4, 12): 0.7338340760527127, (133, 141): 0.0, 
    (131, 133): 0.0, (1, 7): 0.24653882575442454, (134, 134): 0.8606380663808766, 
    (2, 130): 0.0, (3, 4): 0.0, (2, 4): 0.0, (129, 132): 0.7323245483844127}, 
    {(139, 267): 0.0, (132, 132): 0.6055189513596845, (8, 13): 0.23660222273755804, 
    (11, 11): 0.9640744750887139, (259, 259): 0.9854602920895158, 
    (12, 12): 0.6055189513596845, (133, 141): 0.0, (0, 7): 0.0, 
    (266, 266): 0.7916268753073616, (10, 15): 0.9426353078362628, 
    (143, 143): 0.16563375072880127, (1, 6): 0.0, (142, 142): 0.8606380663808766, 
    (11, 14): 0.15441807354502235, (150, 150): 0.487573139345321, 
    (3, 7): 0.9038168616058962, (2, 5): 0.0, (133, 133): 0.4570409226555594, 
    (13, 13): 0.4570409226555594, (149, 149): 0.9319646681708734, 
    (131, 131): 0.9640744750887139, (12, 20): 0.7338340760527127, 
    (138, 142): 0.5160852615777896, (3, 131): 0.0, (139, 140): 0.0, 
    (129, 133): 0.25902842868967213, (137, 140): 0.7323245483844127, 
    (138, 140): 0.0, (139, 142): 0.15441807354502235, (134, 134): 0.8606380663808766, 
    (0, 0): 0.16268264406083688, (136, 141): 0.23660222273755804, (129, 134): 0.0, 
    (131, 259): 0.0, (131, 133): 0.0, (128, 132): 0.7663816443072197, 
    (3, 3): 0.9640744750887139, (3, 5): 0.0, (129, 129): 0.316204105625695, 
    (8, 15): 0.0, (11, 13): 0.0, (130, 135): 0.9426353078362628, 
    (2, 6): 0.5160852615777896, (130, 134): 0.5160852615777896, (9, 14): 0.0, 
    (131, 134): 0.15441807354502235, (136, 264): 0.7406443781286747, 
    (135, 143): 0.6455592536238158, (10, 13): 0.0, (2, 2): 0.5217999096926298, 
    (5, 5): 0.4570409226555594, (1, 4): 0.7323245483844127, 
    (135, 135): 0.16563375072880127, (21, 21): 0.9319646681708734, 
    (138, 138): 0.5217999096926298, (128, 135): 0.0, (258, 258): 0.7916268753073616, 
    (128, 134): 0.12062024006460492, (0, 4): 0.7663816443072197, 
    (137, 141): 0.25902842868967213, (267, 267): 0.9854602920895158, 
    (138, 143): 0.9426353078362628, (8, 12): 0.7663816443072197, 
    (10, 10): 0.5217999096926298, (142, 150): 0.0, (9, 9): 0.316204105625695, 
    (1, 1): 0.316204105625695, (141, 141): 0.4570409226555594, 
    (0, 6): 0.12062024006460492, (128, 256): 0.7406443781286747, 
    (8, 136): 0.7406443781286747, (139, 139): 0.9640744750887139, 
    (136, 143): 0.0, (1, 7): 0.24653882575442454, (11, 15): 0.9038168616058962, 
    (257, 257): 0.7132101466390699, (139, 141): 0.0, (3, 4): 0.0, (2, 4): 0.0, 
    (9, 12): 0.7323245483844127, (265, 265): 0.7132101466390699, 
    (137, 265): 0.7862934888901292, (11, 139): 0.0, (6, 6): 0.8606380663808766, 
    (7, 7): 0.16563375072880127, (137, 137): 0.316204105625695, 
    (4, 12): 0.7338340760527127, (131, 135): 0.9038168616058962, 
    (148, 148): 0.6513698529343085, (137, 142): 0.0, (13, 21): 0.0, 
    (10, 14): 0.5160852615777896, (14, 22): 0.0, (8, 14): 0.12062024006460492, 
    (0, 128): 0.7406443781286747, (151, 151): 0.009710795948105533, (5, 13): 0.0, 
    (15, 23): 0.6455592536238158, (9, 15): 0.24653882575442454, 
    (4, 4): 0.6055189513596845, (10, 12): 0.0, (141, 149): 0.0, 
    (20, 20): 0.6513698529343085, (1, 5): 0.25902842868967213, (130, 133): 0.0, 
    (136, 136): 0.16268264406083688, (138, 141): 0.0, (3, 6): 0.15441807354502235, 
    (14, 14): 0.8606380663808766, (15, 15): 0.16563375072880127, (6, 14): 0.0, 
    (134, 142): 0.0, (1, 129): 0.7862934888901292, (129, 257): 0.7862934888901292, 
    (23, 23): 0.009710795948105533, (128, 128): 0.16268264406083688, 
    (136, 142): 0.12062024006460492, (129, 135): 0.24653882575442454, 
    (128, 133): 0.23660222273755804, (256, 256): 0.045783443174777894, 
    (264, 264): 0.045783443174777894, (130, 130): 0.5217999096926298, 
    (136, 140): 0.7663816443072197, (140, 148): 0.7338340760527127, 
    (0, 5): 0.23660222273755804, (129, 132): 0.7323245483844127, 
    (139, 143): 0.9038168616058962, (11, 12): 0.0, (7, 15): 0.6455592536238158, 
    (2, 7): 0.9426353078362628, (9, 13): 0.25902842868967213, (138, 266): 0.0, 
    (143, 151): 0.6455592536238158, (22, 22): 0.487573139345321, 
    (140, 140): 0.6055189513596845, (131, 132): 0.0, (137, 143): 0.24653882575442454, 
    (132, 140): 0.7338340760527127, (9, 137): 0.7862934888901292, (130, 132): 0.0, 
    (8, 8): 0.16268264406083688, (130, 258): 0.0, (2, 130): 0.0, (10, 138): 0.0}]
    
    The following function call will evaluate that set of instances with both 
    D-Wave's "DW_2000Q_2_1" and the "an_ss_ge_fi_vdeg" simulated annealing 
    solver. Output will be verbose.
    
    The D-Wave system will be configured to run with an annealing time
    of 25 microseconds instead of 20 microseconds, and sample 150 times for each
    problem, which will be treated as a QUBO problem as opposed to an ising problem.
    
    The simulated annealer will sample each problem 1000 times, running 250 sweeps
    each time. It will also make use of 2 threads in parallel. Also, it will follow
    an exponential annealing schedule.
    
    >>> run_instances(instances = y, settings = {"dwave" : True, "sa" : True,
    "dwave_params" : {"num_reads" : 150, "solver" : "DW_2000Q_2_1", 
    "annealing_time": 25, "t": "qubo"}, "sa_params" : {"-r" : 1000, "-s" : 250,
    "-sched": "exp", "-t": 2, "solver": "an_ss_ge_fi_vdeg"}, verbose=True})
    
    The output should include many lines detailing the various energies obtained
    by the solver as well as something similar to (exact values will vary):
        
    {'dwave': [183.14932975230417, 350.4652338429121, 340.81424841408057, 
    3395.987763365066], 'sa': [1363.8740295509374, 6559.503667008985, 
    2656.384930395255, 306550.59810938564]}
            
    ===================
    Settings Parameter:
    ===================
    
    The settings parameter is a dictionary and can contain a wide variety of 
    configuration options for the solvers.
    
        =================
        Activate Solvers:
        =================
        
        The settings parameter must specify which solvers to solve the given problem
        with. (Present version only supports D-Wave's solvers and a classical
        simulated annealing solver). 
        
        If you want to solve the instances on D-Wave, include:
            
            >>> settings = {"dwave" : True, "dwave_params": []}
            
        If you want to solve the instances with the simulated annealer, include:
            
            >>> settings = {"sa" : True, "sa_params": []}
            
        If you want to solve the instances with both solvers, include both:
            
            >>> settings = {"dwave" : True, "sa" : True, "dwave_params": [], "sa_params": []}
            
        ==================
        Configure Solvers:
        ==================
        
        Each solver has a unique list of parameters that can be configured to 
        change how the solver handles the problem instances.
        
            =========================
            D-Wave Solver Parameters:
            =========================
            
            All of the following are parameters that can be passed to the D-Wave
            solver. None are required. For more information about these parameters,
            see https://docs.dwavesys.com/docs/latest/c_rest_api_4.html
            
            ╔════════════════════════════════╦═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
            ║ Key                            ║ Value                                                                                                         ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ anneal_offsets                 ║ JSON array of length equal to                                                                                 ║
            ║                                ║   the solver property anneal_offset_ranges                                                                    ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ anneal_schedule                ║ JSON array of two element                                                                                     ║
            ║                                ║   arrays of floating-point numbers                                                                            ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ annealing_time                 ║ Integer                                                                                                       ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ answer_mode                    ║ [“raw” | “histogram”]                                                                                         ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ auto_scale                     ║ [true | false]                                                                                                ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ beta                           ║ Float                                                                                                         ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ chains                         ║ JSON list of lists of qubit                                                                                   ║
            ║                                ║   indices. The lists are the chains, and they must be disjoint. This parameter                                ║
            ║                                ║   is used only by post-processing.                                                                            ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ flux_biases                    ║ JSON array of floating point                                                                                  ║
            ║                                ║   numbers with length equal to the num_qubits property of the solver. Use 0 for unused or unavailable         ║
            ║                                ║   devices.                                                                                                    ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ flux_drift_compensation        ║ [true | false]                                                                                                ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ initial_state                  ║ JSON array of length equal                                                                                    ║
            ║                                ║   to num_qubits,                                                                                              ║
            ║                                ║   containing {1,0,3} for                                                                                      ║
            ║                                ║   QUBO problems and {1,-1,3} for Ising problems. Use 3 for inactive qubits. This format allows you to take an ║
            ║                                ║   answer from the SAPI clients and send it directly as this parameter.                                        ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ max_answers                    ║ Positive integer                                                                                              ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ num_reads                      ║ Positive integer                                                                                              ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ num_spin_reversal_transforms   ║ Zero or a positive integer                                                                                    ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ postprocess                    ║ [“”|”sampling”|”optimization”]                                                                                ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ programming_thermalization     ║ Integer bounded by solver                                                                                     ║
            ║                                ║   property programming_thermalization_range                                                                   ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ readout_thermalization         ║ Integer bounded by solver                                                                                     ║
            ║                                ║   property readout_thermalization_range                                                                       ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ reduce_intersample_correlation ║ [true | false]                                                                                                ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ reinitialize_state             ║ [true | false]                                                                                                ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ solver                         ║ a string representing the name of the solver to be used to solve the problem. Choices are listed at the       ║
            ║                                ║ bottom of the D-Wave Leap Dashboard page under the subheading: "Supported Solvers". The default value is      ║
            ║                                ║ DW_2000Q_VFYC_2_1.                                                                                            ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ t                              ║ a string indicating the problem type. Choose between "qubo" or "ising". "ising" is the default.               ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣ 
            ║ token                          ║ a string representing the unique API token assigned to the solving account. This value can be found on the    ║
            ║                                ║ left side-bar of the D-Wave Leap Dashboard page under the subheading: "API Token". The defaultvalue is the    ║
            ║                                ║ API Token of our research group's collective account (the one with two hours of computation time per month).  ║
            ╠════════════════════════════════╬═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╣ 
            ║ num_qubits                     ║ an integer (<= 2048, > 0) which indicates the number of qubits to be used in the problem. Using the default   ║
            ║                                ║ value of 2048 reserves all of the available qubits for computation, whether they are actually assigned weights║ 
            ║                                ║ or not.                                                                                                       ║
            ╚════════════════════════════════╩═══════════════════════════════════════════════════════════════════════════════════════════════════════════════╝ 
                                               
            ======================================
            Simulated Annealing Solver Parameters:
            ======================================
            
            All of these parameters can be passed to the simulated annealing
            solver. Both "-s" and "-r" are required. The rest are optional.
            
            ╔════════╦═══════════════════════════════════════════════════════════════════════════╗
            ║ Key:   ║ Value:                                                                    ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -s     ║ a positive integer representing the number of sweeps                      ║
            ║        ║                                                                           ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -r     ║ a positive integer representing the number of repetitions                 ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -b0    ║ [beta0] is initial inverse temperature. Default value: 0.1                ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -b1    ║ [beta1] is final inverse temperature. Default value: 3.0                  ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -r0    ║ [rep0] is starting repetition. Default value: 0                           ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -v     ║  if -v is set, timing                                                     ║
            ║        ║   and some other info is printed. Default value: not set                  ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -g     ║ if -g is set, only the lowest energy solution is printed.                 ║
            ║        ║   Default value: not set                                                  ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -sched ║ specifies a schedule. It can either be lin, exp or be a text              ║
            ║        ║   file on the system which contains an inverse temperature on every line. ║
            ║        ║   Default value: lin                                                      ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ -t     ║ is the number of threads to run in parallel. Default value:               ║
            ║        ║   OMP NUM THREADS                                                         ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ save   ║ (raw string, legacy) if provided, the output of the annealer will be saved║
            ║        ║ to a file with this name in directory                                     ║
            ╠════════╬═══════════════════════════════════════════════════════════════════════════╣
            ║ solver ║ a string naming the specific simulated annealing solver to use, see table ║
            ║        ║ below                                                                     ║
            ╚════════╩═══════════════════════════════════════════════════════════════════════════╝
            
            There are several included simulated annealing solvers available.
            Some of which are not suitable for problems on a Chimera graph.
            Below is a list of available solvers:
                
            ╔═════════════════════╦════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
            ║ an_ms_r1_fi         ║ Multi-spin code for range-1 interactions with magnetic field (approach one)                                        ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ms_r3_nf         ║ Multi-spin code for range-3 interactions without magnetic field (approach one)                                     ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ms_r1_nf_v0      ║ Multi-spin code for range-1 interactions without magnetic field (approach two)                                     ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_ge_fi         ║ Single-spin code for general interactions with magnetic field (fixed number of neighbors)                          ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_ge_fi_vdeg    ║ Single-spin code for general interactions with magnetic field (any number of neighbors)                            ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_ge_nf_bp      ║ Single-spin code for general interactions on bipartite lattices without magnetic field (fixed number of neighbors) ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_ge_nf_bp_vdeg ║ Single-spin code for general interactions on bipartite lattices without magnetic field (any number of neighbors)   ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_rn_fi         ║ Single-spin code for range-n interactions with magnetic field (fixed number of neighbors)                          ║
            ╠═════════════════════╬════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣
            ║ an_ss_rn_fi_vdeg    ║ Single-spin code for range-n interactions with magnetic field (any number of neighbors)</pre>                      ║
            ╚═════════════════════╩════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝
            
            
    ========================
    Interpreting the Output:
    ========================
    
    The output is a dictionary with keys "dwave" and/or "sa" (depending on which
    solvers were used to evaluate the problem instance) which have lists of 
    TTS values associated with them. The TTS values are sorted in the same order
    as the instances were in the instances input parameter.
    
    
    ========================
    Plot Settings Parameter:
    ========================
    
    run_instances accepts a parameter called plot_settings. This allows for the
    results of the various solvers to be depicted graphically. Note that none of these
    parameters are required. By default, plot_settings is set to None, thereby not
    plotting the results. Currently, it supports the following key-value pairs:
    
    ╔════════════════╦═══════════════════════════════════════════════════════════════════════════════════════╗
    ║ Key            ║ Value                                                                                 ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "x_normalize"  ║ True/False - Normalizes the x-axis to instance base units if True                     ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "y_normalize"  ║ True/False - Normalizes the y-axis (TTS values) to a maximum of 1                     ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "title"        ║ String - Title of the Graph                                                           ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "x_label"      ║ String - Label for the x-axis                                                         ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "y_label"      ║ String - Label for the y-axis                                                         ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "x_lim"        ║ Tuple - Limits for the x-axis                                                         ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "y_lim"        ║ Tuple - Limits for the y-axis                                                         ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "fig_size"     ║ Tuple - Dimensions of the Graph (default is (7,7))                                    ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "dwave_colour" ║ String - Plotting parameter for D-Wave plt.plot. Example: "r--" for a red dotted line ║
    ╠════════════════╬═══════════════════════════════════════════════════════════════════════════════════════╣
    ║ "sa_colour"    ║ String - Plotting parameter for SA plt.plot. Example: "g+" for green "+"s             ║
    ╚════════════════╩═══════════════════════════════════════════════════════════════════════════════════════╝
        
    """
    # Import solvers
    from .simulated_annealing.send_to_SA import SA
    from .dwave_interface.Problem_Post import post_problem
    from numpy import log
    
    # Initialize output
    result = {}
    
    # Run (only) selected solvers
    
    # Accomodate multiple inputs
    if type(instances) == dict:
        instances = [instances]
        
    for count, instance in enumerate(instances):
        # Run on D-Wave
        if settings['dwave']:
            
            # Initialize D-Wave results list
            if "dwave" not in result.keys():
                result["dwave"] = []
            
            # Set defaults
            if count == 0:
                solver_dwave = "DW_2000Q_VFYC_2_1"
                t = "ising"
                token = "ARL-251d10ed7bea1e1340e03681887a4e9b36169c81"
                num_qubits = 2048
            
            # Set important optional parameters if not already set
            if "num_reads" not in settings['dwave_params'].keys():
                settings['dwave_params']['num_reads'] = 1
                
            # Remove "non-standard" parameters from settings and overwrite defaults
            if "solver" in settings['dwave_params'].keys():
                solver_dwave = settings['dwave_params'].pop("solver")
            if "t" in settings['dwave_params'].keys():
                t = settings['dwave_params'].pop("t")
            if "token" in settings['dwave_params'].keys():
                token = settings['dwave_params'].pop("token")
            if "num_qubits" in settings['dwave_params'].keys():
                num_qubits = settings['dwave_params'].pop("num_qubits")
                
            # Run instance on D-Wave
            dwave = post_problem(instance, t=t, solver=solver_dwave, token=token, num_qubits=num_qubits, params=settings['dwave_params'])
            
            # Produce tabular output
            if verbose:
                print("D-Wave Solution: ", end="\n\n")
                print("Energies:\t\tCounts:\tSuccess Rate:")
                for i in range(len(dwave["answer"]["energies"])):
                    print(dwave["answer"]["energies"][i], "\t", dwave["answer"]["num_occurrences"][i],
                      "\t", dwave["answer"]["num_occurrences"][i]/settings["dwave_params"]["num_reads"])
                print("\n")
                
            # Compute TTS    
            result['dwave'].append(dwave['answer']['timing']['anneal_time_per_run']/1E6*log(0.01)/log(1-dwave["answer"]["num_occurrences"][0]/settings["dwave_params"]["num_reads"]))
                
        # Run Simulated Annealing
        if settings["sa"]:
            
            # Initialize Simulated Annealing result list
            if "sa" not in result.keys():
                result["sa"] = []
            
            # Set defaults
            if count == 0:
                solver_sa = "an_ss_ge_fi_vdeg"
                save = None
                
            # Set important optional parameters if not already set
            if "-s" not in settings['sa_params'].keys():
                settings['sa_params']["-s"] = 200
            if "-r" not in settings['sa_params'].keys():
                settings['sa_params']["-r"] = 1000
            if "-v" not in settings['sa_params'].keys():
                settings['sa_params']["-v"] = True
                
            # Remove "non-standard" parameters from settings and overwrite defaults
            if "solver" in settings['sa_params'].keys():
                solver_sa = settings['sa_params'].pop("solver")
            if "save" in settings['sa_params'].keys():
                save = settings['sa_params'].pop("save")
                
            # Data first needs to be saved before run
            from os import path, walk, remove
            dir_path = path.dirname(path.realpath(__file__))
            
            # Write instance data into a file in the same directory as this script
            instance_file = open(dir_path+"/temp_file.txt", "w")
            instance_file.write("This file is written automatically for the purposes "+
            "of evaluation by the Simulated Annealing algorithm and is safe to delete.\n")
            for key in instance:
                for value in key:
                    instance_file.write(str(value)+" ")
                instance_file.write(str(instance[key])+"\n")
            instance_file.close()
            
            # Locate the solver directory (by finding solver "an_ms_r1_fi.h")
            for root, dirs, files in walk(dir_path):
                if "an_ms_r1_fi.h" in files:
                     solver_path = path.join(root)
                    
            # Run simulated annealer
            if verbose:
                sa = SA(dir_path, "temp_file.txt", verbose=True, solverDir=solver_path, solver=solver_sa, save=save, solver_params=settings["sa_params"])

            else:
                sa = SA(dir_path, "temp_file.txt", verbose=False, solverDir=solver_path, solver=solver_sa, save=save, solver_params=settings["sa_params"])
                
            # Delete instance data file
            remove(dir_path+"/temp_file.txt")
            
            # Compute TTS
            result['sa'].append(sa[2]/sa[3]*log(0.01)/log(1-sa[1]))
            print(result)
            
    if plot_settings != None:
        
        import matplotlib.pyplot as plt
        
        if "fig_size" not in plot_settings.keys():
            plot_settings["fig_size"] = (7,7)
        if "title" not in plot_settings.keys():
            plot_settings["title"] = "'title' not set"
        if "x_label" not in plot_settings.keys():
            plot_settings["x_label"] = "'x_label' not set"
        if "y_label" not in plot_settings.keys():
            plot_settings["y_label"] = "'y_label' not set"
        if settings["dwave"] and "dwave_colour" not in plot_settings.keys():
            plot_settings["dwave_colour"] = "r-"
        if settings["sa"] and "sa_colour" not in plot_settings.keys():
            plot_settings["sa"] = "b-"
        
        
        if plot_settings["y_normalize"]:
            all_data = []
            for key in result.keys():
                all_data.extend(result[key])
            data_max = max(all_data)
            data_min = min(all_data)
            for key in result.keys():
                result[key] = (result[key]-data_min)/data_max
            
        if not plot_settings["x_normalize"]:
            ns = []
            base_size = 0
            
            if type(instances) == list:
                for instance in instances:
                    for key in instance.keys():
                        if key[0] == key[1]:
                            base_size += 1
                    ns.append(base_size)
                    
            else:
                for key in instances.keys():
                    if key[0] == key[1]:
                        base_size += 1
                ns.append(base_size)
                        
        else:
            if type(instances) == list:
                ns = list(range(len(instances)))
            else:
                ns = [1]
                
        fig = plt.figure(figsize=plot_settings["fig_size"])
        
        if settings["dwave"]:
            plt.plot(ns, result["dwave"], plot_settings["dwave_colour"])
        if settings["sa"]:
            plt.plot(ns, result["sa"], plot_settings["sa_colour"])
            
        plt.title(plot_settings["title"])
        plt.xlabel(plot_settings["x_label"])
        plt.ylabel(plot_settings["y_label"])
        
        if "x_lim" in plot_settings.keys():
            plt.xlim(plot_settings["x_lim"])
        if "y_lim" in plot_settings.keys():
            plt.ylim(plot_settings["y_lim"])
            
    return result
