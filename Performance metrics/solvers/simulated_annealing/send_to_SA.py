def send_to_SA_solver(directory, instance, solver, solver_params={"-s":200,"-r":1000}, verbose=False, save=None):
    
    """
    send_to_solver prepares a problem instance, sends it off to the C++ simulated 
    annealing solver, and handles the output of that solver for subsequent use.
    
    Parameters:
    -----------
    directory: (str) the directory is a string representing the location on disk
               of the problem instance.
    instance: (str) the instance is the file name of the text file containing the 
              the problem details.
    solver: (str) the solver is the name of the simulated annealing solver to use
            to solve the problem, see the README for a list of available solvers
            and their features.
    solver_params: (dict, optional) the solver_params is a dictionary containing 
                   all of the configuration options to pass to the simulated 
                   annealing solver, see the README for details about all of the 
                   options that can be passed to the annealer.
    verbose: (bool, optional) if True, desciptive output will be displayed in the
             Python console as the code executes
    save: (str, optional) if provided, the output of the annealer will be saved to 
          a file with this name in directory
          
    Returns:
    --------
    A text file in directory containing the solution information, and/or a Python
    variable containing all of the solution information from the annealer.
    """
    
    # Necessary imports
    from subprocess import Popen, PIPE
    from os import remove
    
    # Perfom parameter checks
    if type(directory) != str:
        raise Exception("The directory must be a string!")
    if type(instance) != str:
        raise Exception("The instance must be a string!")
    if type(solver_params) != dict:
        raise Exception("The solver parameters must come in a dictionary!")
    if type(verbose) != bool:
        raise Exception("Verbose must be a boolean.")
    if save != None and type(save) != str:
        raise Exception("In order to save the result, a file name must be provided as a string.")
    if "-s" not in solver_params.keys() or "-r" not in solver_params.keys():
        raise Exception("Solver parameters must include \"-s\" and \"-r\" as keys.")
        
    acceptable_solvers = ["an_ms_r1_nf", "an_ms_r1_fi", "an_ms_r3_nf", "an_ms_r1_nf_v0",
                          "an_ss_ge_fi", "an_ss_ge_fi_vdeg", "an_ss_ge_nf_bp", 
                          "an_ss_ge_nf_bp_vdeg", "an_ss_ge_fi_bp_vdeg", 
                          "an_ss_rn_fi", "an_ss_rn_fi_vdeg"]

    if solver not in acceptable_solvers:
        print("WARNING: Solver not recognized! Defaulting to an_ms_r1_nf. Choose one of the following solvers:", acceptable_solvers, end="\n\n")
        solver = "an_ms_r1_nf"
        
    acceptable_params = ["-s", "-r", "-b0", "-b1", "-r0", "-v", "-g", "-sched", "-t"]
    
    for key in solver_params.keys():
        if key not in acceptable_params:
            raise Exception(key, " is an undefined parameter! The following is a list of acceptable parameters: ", acceptable_params)
        if key == "-s" or key == "-r" or key == "-r0":
            if type(solver_params[key]) != int or solver_params[key] < 1:
                raise Exception("The parameter ", key, " must be a positive integer.")
        if key == "-b0" or key == "-b1":
            if (type(solver_params[key]) != int and type(solver_params[key]) != float) or solver_params[key] < 0:
                raise Exception("The parameter ", key, " must be positive integer or float.")
        if key == "-v" or key == "-g":
            solver_params[key] = None
            if verbose:
                print(key, " is set.", end="\n\n")
        if key == "-sched":
            if solver_params[key] != "lin" and solver_params[key] != "exp" and solver_params[key][:4:-1] != "txt.":
                solver_params["-sched"] = "lin"
                print("Schedule parameter not understood, defaulting to \"lin\". Choose from \"lin\", \"exp\" or the name of a "+
                      "text file which contains an inverse temperature on every line.", end="\n\n")
        if key == "-t":
            if type(solver_params[key]) != int and solver_params[key] < 1:
                raise Exception("Number of parallel threads -t must be a positive integer.")
    
        
    # Prepare instance name
    if instance[:4:-1] == "txt.":
        instance = instance[:len(instance)-3]
        
    # Write Batch File
    write_batch(directory, instance, solver, solver_params)
    
    if verbose:
        print("Batch File Created! File contents:", end="\n\n")
        file = open("solve_"+instance+"_batch.bat", "r")
        print(file.read(), end="\n\n")
        file.close()
    
    # Run the batch file and decode the output
    p = Popen("solve_"+instance+"_batch.bat", cwd=directory, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print(stderr)
    stdout = stdout.decode("utf-8")
    p_status = p.wait()
    remove("solve_"+instance+"_batch.bat")
    
    if verbose:
        print(stdout)
        print("solve_"+instance+"_batch.bat successfully deleted!", end="\n\n")
    
    # Save the output if desired
    if save != None:
        if save[:4:-1] != "txt.":
            save = save+".txt"
        f = open(save, "a+")
        f.write("Solution for Problem Instance: " + instance + "\n\n")
        f.write("    Energy:  Frequency: Success Rate:   Instance Name:\n")
        f.write(stdout+"\n\n")
        f.close()
        
    return stdout
    
def write_batch(directory, instance, solver, solver_params):
    
    """
    write_batch writes a batch file which, when run, will call the appropriate
    C++ functions and run the simulated annealing algorithm on the provided problem
    instance
    
    Parameters:
    -----------
    directory: (str) the directory is a string representing the location on disk
               of the problem instance.
    instance: (str) the instance is the file name of the text file containing the 
              the problem details.
    solver: (str) the solver is the name of the simulated annealing solver to use
            to solve the problem, see the README for a list of available solvers
            and their features.
    solver_params: (dict, optional) the solver_params is a dictionary containing 
                   all of the configuration options to pass to the simulated 
                   annealing solver, see the README for details about all of the 
                   options that can be passed to the annealer.
    Returns:
    --------
    A batch file in directory which can be called to run the appropriate C++
    simulated annealing algorithm.
    """
    
    f= open("solve_"+instance+"_batch.bat", "w+")
    f.write("@echo off\n")
    f.write("cd: "+directory+"\n")
    f.write(solver+" -l "+instance+".txt ")
    for param in solver_params.keys():
        if param != "-v" and param != "-g":
            f.write(param+" "+str(solver_params[param])+" ")
        else:
            f.write(param+" ")
    f.close()
    
send_to_SA_solver("E:\Physics 598\Run_SA","instance","an_ms_r1_nf",{"-s":200,"-r":1000}, save="output",verbose=True)
send_to_SA_solver("D:\Physics 598\Run_SA","503_pm_nf_0000","an_ms_r1_nf",{"-s":200,"-r":1000,"-t":4}, save="output",verbose=True)