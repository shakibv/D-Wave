# Import Packages (Requests takes care of http stuff)
import json
import requests

"""
Tutorial:
---------
If necessary, to better understand how the post_problem function works, try 
running the following code. This will create a problem dictionary and run it
on the D-Wave system.

qubit_biases = {(0, 0): 0.3333, (1, 1): -0.333, (4, 4): -0.333, (5, 5): 0.333}
coupler_strengths = {(0, 4): 0.667, (0, 5): -1, (1, 4): 0.667, (1, 5): 0.667}
Q = dict(qubit_biases)
Q.update(coupler_strengths)
print(post_problem(Q))

"""

def post_problem(problem, t="ising", solver="DW_2000Q_VFYC_2_1", token="ARL-251d10ed7bea1e1340e03681887a4e9b36169c81", params={}, num_qubits=2048):
    """
    post_problem is a function which poses a given problem to a D-Wave solver
    and returns the solver's response.
    
    Parameters:
    -----------
    problem: a dictionary containing all of the assignments of weights to the 
             qubits and strengths to the edges. Each key-value pair in the dictionary
             should be structured as follows:
                 
                 (first qubit index, second qubit index) : weight/strength
                 
             Key-value pairs which are used to assign weights to a particular qubit
             must use the same qubit index twice in the key portion of the key-value
             pair. For example, the following key-value pair assigns a weight of 
             0.667 to the qubit at index 7:
                
                 (7, 7) : 0.667
                
             Key-value pairs which are used to assign strengths to particular edges
             between qubits must used the indices of the two connected qubits in the
             key portion of the key-value pair. For example, the following
             key-value pair assigns a strength of -0.5 to the edge connecting
             the qubit at index 1 to the qubit at index 4:
                
                 (1, 4) : -0.5
             
    t: (optional) a string indicating the problem type. Choose between "qubo" 
       or "ising". "ising" is the default.
       
    solver: (optional) a string representing the name of the solver to be used
            to solve the problem. Choices are listed at the bottom of the D-Wave
            Leap Dashboard page under the subheading: "Supported Solvers". The 
            default value is DW_2000Q_VFYC_2_1.
            
    token: (optional) a string representing the unique API token assigned to the
           solving account. This value can be found on the left side-bar of the
           D-Wave Leap Dashboard page under the subheading: "API Token". The default
           value is the API Token of our research group's collective account
           (the one with two hours of computation time per month).
           
    params: (optional) a dictionary containing additional parameters to be passed
            to the solver. See: https://docs.dwavesys.com/docs/latest/c_rest_api_4.html
            for a complete list of parameters that can be passed to the solver
            and their definitions. The default is an empty dictionary, thereby
            passing no additional parameters to the solver.
            
    num_qubits: (optional) an integer (<= 2048, > 0) which indicates the number of 
                qubits to be used in the problem. Using the default value of 2048
                reserves all of the available qubits for computation, whether
                they are actually assigned weights or not.
                
    Returns:
    --------
    a dictionary containing the following key-value pairs:
        
        KEY:               | VALUE:
        -----------------------------------------------------------------------
        status             | a string indicating the progress that the
                           | solver has made towards solving the posed
                           | problem. COMPLETED indicates a successful
                           | computation and FAILED indicates an 
                           | unsuccessful computation.
                           |
        earliest_estimated_| a string representing the time and date of the 
        completion         | earliest estimated completion of the submitted
                           | problem
                           |
        solved_on          | a string representing the time and date of the
                           | completion of the solution to the posted problem
                           |
        solver             | a string denoting the solver used to solve the
                           | posted problem. Will be identical to the solver
                           | parameter
                           |
        submitted_on       | a string representing the time and date of the 
                           | submission of the problem
                           |
        answer             | a dictionary containing many key-value pairs which
                           | characterize the solution to the posed problem 
                           | (see below)
                           |
        latest_estimated_  | a string representing the time and date of the 
        completion         | latest estimated completion of the problem
                           |
        type               | a string indicating the problem type. Will be 
                           | identical to the t parameter
                           |
        id                 | a string containing the problem id used to log 
                           | the problem activity on the D-Wave Leap account
                           | associated with the given token
                           |
    
    The answer dictionary contains the following key-value pairs:
        
        KEY:               | VALUE:
        -----------------------------------------------------------------------
        num_variables      | an integer representing the number of qubits 
                           | reserved for the problem. Will be identical to the
                           | num_qubits parameter.
                           |
        format             | a string which will always read "qp".
                           |
        energies           | a base-64-encoded string of energies, each a 
                           | little-endian 8-byte floating-point number (doubles).
                           |
        num_occurances     | a base-64–encoded string of the number of occurrences 
                           | of each solution. The numbers are 4-byte little-endian 
                           | integers. Answers are ordered by energy.
                           |
        active_variables   | a base-64–encoded string of the indices of the 
                           | problem’s active variables. The indices are 4-byte 
                           | little-endian integers.
                           |
        solutions          | a base-64–encoded string of bit-packed solutions 
                           | (with 0 = -1 for Ising problems). Bits are in 
                           | little-endian order. Each solution is padded to 
                           | end on a byte boundary and contains values for 
                           | active qubits only.
                           |
        timing             | a solver-specific JSON object describing the 
                           | reported time that the solver took to handle the 
                           | problem.
                           |
    
    For more information about the returned information, see D-Wave's documentation:
        
    https://docs.dwavesys.com/docs/latest/c_rest_api_3.html
    """
    # Convert the problem dictionary to an appropriately formatted string
    problem = convert_to_text(problem, num_qubits)
    
    # Create the header for the url request
    my_header = {"X-Auth-Token": token}
    
    # This is the base url for the D-Wave API when submitting problems.
    url = "https://cloud.dwavesys.com/sapi/problems"
    
    # Data is the package to be sent to the D-Wave Server, convert it to JSON
    my_data = {"solver" : solver, "data" : problem, "type": t, "params": params}
    my_data = json.dumps(my_data)
    
    # Sends the request to the D-Wave server and returns its response
    r = requests.post(url, data=my_data, headers=my_header)
    return json.loads(r.text)
    
def convert_to_text(problem, num_qubits):
    """
    convert_to_text takes a dictionary of the weight and strength assignments
    within the Chimera graph and produces a string which is properly formated
    to be submitted to the D-Wave server for computation on a D-Wave quantum
    computer.
    
    Parameters:
    -----------
    problem: a dictionary containing all of the assignments of weights to the 
             qubits and strengths to the edges.
             
    num_qubits: (optional) an integer (<= 2048, > 0) which indicates the number of 
                qubits to be used in the problem. Using the default value of 2048
                reserves all of the available qubits for computation, whether
                they are actually assigned weights or not.
                
    Returns:
    --------
    master_str: a string which contains all of the information within the problem
                dictionary but it formatted such that it can be submitted to a 
                D-Wave solver.
    """
    master_str = str(num_qubits) + " " + str(len(problem)) + "\n"
    for value in problem:
        master_str = master_str + str(value[0]) + " " + str(value[1]) + " " + str(problem[value]) + "\n"
    return master_str
