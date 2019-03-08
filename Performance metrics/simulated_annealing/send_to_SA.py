def SA(directory, instance, solver="an_ms_r1_nf", solver_params={"-s":200,"-r":1000}, verbose=False, save=None, solverDir=None):
    
    """
    send_to_solver prepares a problem instance, sends it off to the C++ simulated 
    annealing solver, and handles the output of that solver for subsequent use.
    
    Parameters:
    -----------
    directory: (raw str) the directory is a string representing the location on disk
               of the problem instance.
    instance: (str) the instance is the file name of the text file containing the 
              the problem details.
    solver: (str optional) the solver is the name of the simulated annealing solver to use
            to solve the problem, see the README for a list of available solvers
            and their features. Defaults to an_ms_r1_nf.
    solver_params: (dict, optional) the solver_params is a dictionary containing 
                   all of the configuration options to pass to the simulated 
                   annealing solver, see the README for details about all of the 
                   options that can be passed to the annealer.
    verbose: (bool, optional) if True, desciptive output will be displayed in the
             Python console as the code executes
    save: (raw str, optional) if provided, the output of the annealer will be saved to 
          a file with this name in directory
    solverDir: (raw str, optional) only needed if the instance is not saved in the same
               directory as the solvers. Specifies the absolute system path to
               the solver's directory
          
    Returns:
    --------
    A text file in directory containing the solution information, and/or the minimum
    energy solution and number of sweeps as a tuple.
    
    ==============
    Example Usage:
    ==============
    
    Simple Example:
    ---------------
    
    Suppose that you had an instance named my_instance.txt in the directory 
    C:\\users\\John\\Documents. You want to solve for the minimum energy of this
    instance. You have made it as simple as possible by putting all of the 
    solver source codes in the same directory as the instance file:
    C:\\users\\John\\Documents. The following function call would solve for the energy:
        
    """
        
    #energy = SA(r"C:\users\John\Documents", "my_instance.txt")
       
    """
    
    Full Functionality Example:
    ---------------------------
    
    Suppose that you had an instance named my_instance.txt in the directory
    C:\\users\\John\\Documents\\instances that you wanted to solve using the 
    an_ss_ge_nf_bp_vdeg solver using 200 sweeps, 1000 repetitions, on an exponetial
    schedule, on 4 processor cores simultaneously and only return the lowest energy
    solution. You wanted to save the result to a text file called my_instance_solution.txt
    in the directory C:\\users\\John\\Documents\\instance_solutions. The an_ss_ge_nf_bp_vdeg
    solver was located in the folder C:\\users\\John\\Documents\\source_code. The
    following function call would solve for the energy:
        
    """
    
    #energy, sweeps = SA(r"C:\users\John\Documents\instances", "my_instance.txt", 
    #"an_ss_ge_nf_bp_vdeg", {"-s":200, "-r":1000, "-sched": "exp", "-t":4, "-g":None}, 
    #False, r"C:\users\John\Documents\instance_solutions\my_instance_solution.txt", 
    #r"C:\users\John\Documents\source_code")
    
    # Necessary imports
    from subprocess import Popen, PIPE
    from os import remove
    from re import split
    
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
    if instance[:len(instance)-5:-1] == "txt.":
        instance = instance[:len(instance)-4]
        
    # Check to see if instance is in solver directory, if not, make a copy of
    # the instance in that directory
    
    dirs_dont_match = bool(solverDir != None and solverDir != directory)
    
    if dirs_dont_match:
        instance_file = open(directory+"/"+instance+".txt", "r")
        instance_file_data = instance_file.read()
        instance_file.close()
        copy_instance = open(solverDir+"/"+instance+".txt", "w")
        copy_instance.write(instance_file_data)
        copy_instance.close()

    # Build string
    command = str(solver+" -l "+instance+".txt ")
    for param in solver_params.keys():
        if param != "-v" and param != "-g":
            command = command + str(param+" "+str(solver_params[param])+" ")
        else:
            command = command + str(param+" ")
    
    # Run the batch file and decode the output
    if dirs_dont_match:
        p = Popen(command, cwd=solverDir, stdout=PIPE, stderr=PIPE, shell=True)
    else:
        p = Popen(command, cwd=directory, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    stdout = stdout.decode("utf-8")
    p_status = p.wait()
    if dirs_dont_match:
        remove(solverDir+"/"+instance+".txt")
        
    if verbose:
        print("Output:\n\n", stdout, "Error: ", stderr)
    
    # Save the output if desired
    if save != None:
        if save[:4:-1] != "txt.":
            save = save+".txt"
        f = open(save, "a+")
        f.write("Solution for Problem Instance: " + instance + "\n\n")
        f.write("    Energy:  Frequency: Success Rate:   Instance Name:\n")
        f.write(stdout+"\n\n")
        f.close()
        
    # Format output to return only lowest energy and sweep number s used for
    # timing? To be adjusted as we figure out optimal annealing time etc.
    output = split("\r\n|\t| ", stdout)
    output = list(filter(None, output))
    
    return int(output[0]), solver_params["-s"]
    
print(SA(r"D:\Physics 598\GitHub\solvers\simulated_annealing", "instance.txt"))
