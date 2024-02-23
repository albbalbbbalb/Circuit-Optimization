import numpy as np
from itertools import product
from mip import *
import os


def make_circuit_lp(truth='01010011', gates='NOT NAND NAND NAND',
                    filename=None):
    '''
    A function that creates a .lp file, a standard for use in Mixed
    Integer Linear Programming (MILP) for optimization. The resulting
    .lp file will be used to solve for the optimal number of logic gates
    in a switching circuit.
    
    The default arguments produce circuit.lp, solving the MILP problem
    will result in an optimal 4 gate 2x1 multiplexer circuit.
    
    The formulation for the MILP system we use here follows
    the results of Muroga 91[1].
    
    [1] Muroga, Saburo "Computer-Aided Logic Synthesis for VLSI chips",
        1991 Academic Press, Inc, isbn:0-12-012132-8,
        doi:10.1016/S0065-2458(08)60245-4
    
    Arguments:
    --------------------------------------------------------------------
    
    truth -- a string of truth values that represent the output of a
             binary function f:{0,1}^n ->{0,1} of n variables. For
             example the input [0,1,0,1,0,0,1,1] corresponds to the
             following truth table:
             
             X1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1
             X2 | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1
             X3 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1
             ---|---|---|---|---|---|---|---|---
             f  | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 1
             
             The above truth table corresponds to a 2x1 multiplexer.
             
    gates -- a string of logic gate names separated by spaces,
             e.g., "NOT NAND NAND NAND", that are to be included in the
             circuit. The order of the gates in the string matters for
             the success of the algorithm.
        
    filename -- name of the .lp file to save to. If given None, will
                save to circuit.lp

    
    Effect: produces a .lp file in the current working directory
    -------------------------------------------------------------------
    '''
    M = Model(sense=MINIMIZE, name=f'circuit for logic function {truth} \
                                     with gates {gates}')
    
    n = int(np.log2(len(truth))) # number of inputs
    truth = [int(i) for i in list(truth)]
    gates = gates.split(' ')
    gates = [ 0 if i == 'NOT' else -1 for i in gates]
    R = len(gates)               # number of NOR gates in the circuit
    A = n + R                    # upper bound we'll use later
    
    table = np.array(list(product([0,1], repeat=n)))
    table = np.fliplr(table)
    
    # generate the sets of indices
    U = []
    for el in range(1,n+1):
        U.append([])
        for k in range(1,R+1):
            U[el-1].append(M.add_var(name=f'u.{el}.{k}', var_type=BINARY))
    V = []
    for k in range(1,R+1):
        V.append([])
        for i in range(1,k):
            V[k-1].append(M.add_var(name=f'v.{i}.{k}', var_type=BINARY))
    P = []
    for k in range(1,R):
        P.append([])
        for j in range(1,2**n+1):
            P[k-1].append(M.add_var(name=f'p.{k}.{j}', var_type=BINARY))
    r = []
    for k in range(1,R+1):
        r.append([])
        for i in range(1,k):
            r[k-1].append([])
            for j in range(1,2**n+1):
                r[k-1][i-1].append(M.add_var(name=f'r.{i}.{k}.{j}',
                                             var_type=BINARY))
    
    # the objective function
    Uobj = xsum( xsum(k for k in el) for el in U)
    Vobj = xsum( xsum(j for j in k) for k in V)
    M.objective = Uobj + Vobj
    
    # construct constraint inequalities, there are a lot of them
    for k, j in product(range(1,R), range(1,2**n+1)):
        con = 0
        con += -xsum( table[j-1,el-1]*U[el-1][k-1] for el in range(1,n+1))
        con += -xsum(r[k-1][i-1][j-1] for i in range(1,k))
        M += con >= gates[k-1] -A*(1-P[k-1][j-1])
        con = 0
        con += xsum( table[j-1,el-1]*U[el-1][k-1] for el in range(1,n+1))
        con += xsum(r[k-1][i-1][j-1] for i in range(1,k))
        M += con >= 1-gates[k-1] -A*P[k-1][j-1]
        
    for j, val in enumerate(truth):
        con = 0
        if val == 1:
            con += -xsum( table[j,el-1]*U[el-1][R-1] for el in range(1,n+1))
            con += -xsum( r[R-1][i-1][j] for i in range(1,R) )
            M += con >= gates[k-1]
        else:
            con += xsum( table[j,el-1]*U[el-1][R-1] for el in range(1,n+1))
            con += xsum( r[R-1][i-1][j] for i in range(1,R) )
            M += con >= 1-gates[k-1]

        
    for k, j in product(range(2,R+1), range(1,2**n+1)):
        for i in range(1,k):
            M += P[i-1][j-1] + V[k-1][i-1] - r[k-1][i-1][j-1] <=1
            M += P[i-1][j-1] + V[k-1][i-1] - 2*r[k-1][i-1][j-1] >=0

    for k in range(1,R+1):
        con = 0
        con += xsum( U[el-1][k-1] for el in range(1,n+1)) 
        con += xsum( V[k-1][i-1] for i in range(1,k))
        M += con <= 1 - gates[k-1]
    
    if filename != None:
        M.write(f'{filename}.lp')
    else:
        M.write('circuit.lp')


def solve_circuit_lp(truth, gates, filename=None, verbose=False,
                     delete_lp=True, delete_log=True):
    '''
    A function that solves a Mixed Integer Linear Programming (MILP)
    problem specified by a .lp file for a circuit. This function
    requires the software suite SCIP [1]. The SCIP solver output is
    parsed to recover the relevant values.

    [1] Solving Constraint Integer Programs (SCIP)
        https://scipopt.org/

    Arguments
    -------------------------------------------------------------------
    
    truth -- a string of truth values that represent the output of a
             binary function f:{0,1}^n ->{0,1} of n variables. For
             example the input [0,1,0,1,0,0,1,1] corresponds to the
             following truth table:
             
             X1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1
             X2 | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1
             X3 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1
             ---|---|---|---|---|---|---|---|---
             f  | 0 | 1 | 0 | 1 | 0 | 0 | 1 | 1
             
             The above truth table corresponds to a 2x1 multiplexer.
             
    gates -- a string of logic gate names separated by spaces,
             e.g., "NOT NAND NAND NAND", that are to be included in the
             circuit. The order of the gates in the string matters for
             the algorithm, e.g. ending with a NOT gate frequently
             leads to an infeasible system. Currently only supports
             NAND and NOT gates.
        
    filename -- name of the .lp file to save to. If given None, will
                save to circuit.lp.

    verbose -- flag that determines whether the intermediate output
               of SCIP is printed.

    delete_lp -- flag that determines whether the .lp file is deleted
                 after the computations are finished.

    delete_log -- flag that determines whether the .log file created
                  by SCIP is deleted.

    Output
    -------------------------------------------------------------------
    connections -- a list of strings encoding the connections for the
                   solution. The following convention is used:

                   'u.x.y' = connect external variable x to an 
                             input of gate y
                   'v.x.y' = connect the output of gate x to an
                             input of gate y

                    where x and y are some integers.

                    If connections is an empty list, the problem
                    is infeasible and has no solution, check the .log
                    file for more info.
    -------------------------------------------------------------------
    '''

    if filename == None:
        filename = 'circuit'
    if filename[-3:] == '.lp':
        filename = filename[:-3]
        
    make_circuit_lp(truth, gates, filename)
    
    if os.path.exists(f"{filename}.log"):
        os.system(f"rm {filename}.log")
    
    if verbose:
        os.system(f'./scip -f {filename}.lp -l {filename}.log')
    else:
        os.system(f"./scip -q -f {filename}.lp -l {filename}.log")

    connect = []

    with open(f'{filename}.log') as f:
        for line in f:
            line = line.split(' ')
            if '\t(obj:1)\n' in line:
                connect.append((line[0], int(line[-2])))
                
    if delete_lp:
        os.system(f'rm {filename}.lp')
    if delete_log:
        os.system(f'rm {filename}.log')
        
    return connect
    

def print_connections(connections):
    '''
    Short helper function for printing the output of solve_circuit_lp
    in a readable form.

    Arguments
    -------------------------------------------------------------------
    
    connections -- the output of the function solve_circuit_lp

    Effects: prints the connections nicely
    -------------------------------------------------------------------
    '''

    for (item,num) in connections:
        v, n, m = item.split('.')
        if v == 'u':
            print(f'Input\t {n} connects to gate {m}')
        elif v == 'v':
            print(f'Gate\t {n} connects to gate {m}')


if __name__ == "__main__":

    gates = "NAND NAND NAND NAND"
    truth = '01010011'
    filename = 'test'
    
    print('Simple test of solve_circuit_lp')
    print('The function we are computing is', truth)
    print('The gates we chose are', gates)
    print('We save the .lp file and .log files at test.lp and test.log')
    print('We print the connections at the end.\n')
    
    connect1 = solve_circuit_lp(truth, gates, filename,
                                delete_lp=False, delete_log=False)
    print_connections(connect1)
