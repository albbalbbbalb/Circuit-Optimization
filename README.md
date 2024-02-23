# Circuit-Optimization

A collection of Python functions for creating MILP (Mixed Integer Linear Programming) problems the solution of which correspond to optimal combinational logic circuit designs. 

## Introduction

Combinational logic circuits or switching circuits are idealized electrical circuits that model Boolean functions, that is, functions of the form $f:\\{0,1\\}^n\to\\{0,1\\}$ where $n$ is the number of inputs. Some well known switching circuits are AND, OR, NOT, NAND, NOR, XOR, XNOR, also called gates, to name some 1 or 2 input examples. Combinational logic circuits are widely used in the design of modern technology, both at the basis of personal computers and in other specialized systems. The optimization of combinational circuits has been a priority in academia and in industry since the beginning of the digital age, since both the resources and manufacturing process can be costly. This project is a personal challenge to research this subject and find a feasible solution at least for small circuits.

The NAND gate (as well as the NOR gate) has the property of functional completeness. That is, any Boolean function can be designed as the composition of some number of NAND gates, this was proven by Henry M. Sheffer in 1913. The next natural question: what is the minimal number of NAND gates required for a specific Boolean function? In their paper [1], Muroga outlines a method for constructing any Boolean function using some number of logical gates. The code provided in this package is essentially what Muroga constructs but adapted from the NOR example to use only NAND gates.

[1]: Muroga, Saburo "Computer-Aided Logic Synthesis for VLSI chips", 1991 Academic Press, Inc, isbn:0-12-012132-8, doi:10.1016/S0065-2458(08)60245-4

## Dependencies

Make sure that the Python packages `mip`, `numpy`, `itertools` and `os` are installed. We suggest to install the [SCIP](https://scipopt.org/) suite of MILP optimization tools and leave the binary in the current working directory. The `mip` package includes its own solver however it did not seem beefy enough for these purposes.

## How to use this package?

If SCIP is installed, run `python circuit_optimization.py` to check that all is in order. To use the `mip` solver, please, refer to the relevant documentation. The function `make_circuit_lp` creates a `.lp` file in the current working directory which can then be passed on to any MILP solver that takes this standard. The function `solve_circuit_lp` makes a call to `make_circuit_lp` and uses SCIP to solve the problem. The output of `solve_circuit_lp` can then be passed to `print_connections` to display the result in a readable way. Solving for the 2x1 multiplexer circuit:

```
gates = "NAND NAND NAND NAND"
truth = '01010011'
filename = 'test'
connect1 = solve_circuit_lp(truth, gates, filename, delete_lp=False, delete_log=False)
print_connections(connect1)
```

`gates` specifies that the circuit should have 4 NAND gates, and `truth` specifies the output of the Boolean function we care about, e.g., from left to right $f(0,0,0)=0$, $f(0,0,1) =1$, and $f(0,1,0) = 0$ and so on. We also keep the `.lp` file and the `.log` file from SCIP in the directory for examination. The output of the above code is 

```
Input	 1 connects to gate 1
Input	 1 connects to gate 2
Input	 2 connects to gate 3
Input	 3 connects to gate 1
Input	 3 connects to gate 3
Gate	 1 connects to gate 2
Gate	 2 connects to gate 4
Gate	 3 connects to gate 4
```
where the order of the gates corresponds to their order from left to right in the string `gates`.

## Known limitations and TODOs

- To confidently determine that the minimum number of NAND gates for a particular Boolean function is $N$, the method of Muroga actually requires to verify that it cannot be done with $1,2,\dots, N-1$ NAND gates, so for now `solve_circuit_lp` should be called in a for loop to build up to $N$. It should be possible to bake a variable number of NAND gates into the formulation of the `.lp` file.
- Currently only NAND and NOT gates are implemented. We'd like to extend with more gate options.
- The method of Muroga assumes that the output of a gate can go into at most 1 input of another gate (for as many gates as needed), we'd like to relax that condition. For now, if you want to connect one output twice into a NAND gate, you have to just replace it with a NOT gate before running the solver, e.g., repeat the above example but with `gates="NOT NAND NAND NAND"`
- It would be nice to include other printing options, e.g., printing a diagram instead.

## Further comments

If you see ways that this project can be improved, let me know what you think. I can either add it to the TODOs or if you are willing you can contribute as well. 
