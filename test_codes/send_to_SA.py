import numpy as np
def SA(directory, instance, size, solver="an_ss_ge_fi_bp_vdeg", solver_params={"-s":2000,"-r":10000}, verbose=False, save=None, solverDir=None):
    
    
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

    # Run the C++ solver and interpret the output
    if dirs_dont_match:
        p = Popen('./'+command, cwd=solverDir, stdout=PIPE, stderr=PIPE, shell=True)
    else:
        p = Popen('./'+command, cwd=directory, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    stdout, stderr = stdout.decode("utf-8"), stderr.decode("utf-8")
    p_status = p.wait()
    if dirs_dont_match:
        remove(solverDir+"/"+instance+".txt")
    
    # Print verbose output
    if verbose:
        from re import search
        print(stdout[:search("#work done in [0123456789]\.[0123456789]* s", stdout).end()], end="\n\n")
        print("Simulated Annealing Solution:\n\nEnergies:\tCounts:\tSuccess Rate:\tInstance File:\n", stdout[search("#work done in [0123456789]\.[0123456789]* s", stdout).end():])
        if stderr != '':
            print("Error: ", stderr)

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
    p_s = float(output[4])
    p_d = 0.99
    f = solver_params['-s']*size*solver_params['-r']/float(output[0])
    if p_s == 1:
        T = 0
    else:
        T = (solver_params['-s']*size/f)*np.log2(1-p_d)/np.log2(1-p_s)
    
    #return float(output[2]), T
    return T

#print(SA(r"/Users/archie/Google Drive/D-Wave-master/test_codes/bin", "instance.txt", 128))


