QuEST-TN Tutorial
======

**Table of Contents**
- [Background](#background)
- [Coding](#coding)
- [Compiling](#compiling)
- [Running](#running)
- [Test circuits](#Test-circuits)

# Background

## Motivation

The fundamental idea behind the tensor network extension is to save memory and time by splitting the full system of multiple qubits, represented by single a Qureg object in QuEST, into multiple independent collections of qubits, each with its own Qureg. As an example, for a system of 30 qubits, 2^30 complex doubles need to be stored in the usual QuEST implementation. If this can be split into two systems of 15 qubits, only 2*2^15 elements need to be stored. Any gates applied to one of these systems will also only need to process 2^15 elements. 

It can be useful to treat these subsystems as tensors of rank numQubits, with each tensor index ranging from {0,1} and corresponding to a single qubit.

## Entanglement

Naturally, this does not capture any entanglement between qubits in different subsystems. However, an entangling operation between qubits in different tensors can be represented by introducing an additional 'virtual' qubit to each of the two tensors, as described below. It is necessary to balance the savings due to splitting a system of qubits into multiple subsystems against the number of entangling operations between the two systems, each introducing an additional qubit. The method will be most useful for circuits with weakly entangled subsystems.

### Example of a controlledNot(system, controlQubit, targetQubit) operation between qubits in different tensors

For a system of 4 qubits split into 2 tensors, with tensor1 containing qubits (0,1) and tensor2 containing qubits (2,3), the operation controlledNot(globalQureg, 1, 2) can be represented by 

'''
controlledNot(tensor1, 1, virtualTargetQubit)
controlledNot(tensor2, virtualControlQubit, 2)
'''

Where the virtualTargetQubit is part of the tensor1 system and initialised in the |0> state, and the virtualControlQubit is part of the tensor2 system and is initialized in the (non-normalized) |0> + |1> state. 

The Qureg objects in the two tensors will not be normalized at this stage. However, a normalised system corresponding to the output that you would find doing the same operation on a single Qureg object can be recovered by performing a tensor contraction, contracting across the index corresponding to the virtual qubit in each tensor. 

This can be generalised to any control operation by initialising the virtual qubits in the same way but applying the corresponding control operation, and to multiple entangling operations by tracking pairs of virtual qubits across tensors. When contracting two tensors with multiple entangling operations between them, all pairs of virutal qubit indices are contracted out.

# Coding

To use the QuEST-TN extension in your C or C++ code, include

```C
#include <QuEST.h>
#include <QuEST_tn.h>
```

The currently available API for operations on tensor networks is available [here](https://aniabrown.github.io/QuEST-TN/QuEST__tn_8h.html). For operations on plain QuEST Qureg objects, see the API [here](https://quest-kit.github.io/QuEST/QuEST_8h.html) and the [QuEST tutorial](examples/README.md). 

Here's a very simple circuit which applies a controlled not on a two qubit state split over two
tensors and reports the resulting tensor network structure. 
 
```C
#include <QuEST.h>
#include <QuEST_tn.h>

int main(int narg, char *varg[]) {

  // load QuEST environment. Do this exactly once as the first step in your code.
  QuESTEnv env = createQuESTEnv();
  
  int numPqPerTensor[2] = {1, 1};
  int numVqPerTensor[2] = {1, 1};
  
  TensorNetwork tn = createTensorNetwork(2, numPqPerTensor, numVqPerTensor, env);
  printTensorNetwork(tn);
  tn_controlledNot(tn, 0, 1);
  printTensorNetwork(tn);
	
  // unload QuEST environment. Do this exactly one as the last step in your code. 
  destroyQuESTEnv(env); 
  return 0;
}
```
----------------------------

# Compiling

QuEST uses CMake (3.1 or higher) as its build system.

To compile, make sure your circuit code is accessible from the root QuEST directory.
In the root directory, initially build using
```bash
mkdir build
cd build
cmake -DUSER_SOURCE="myCode1.c;myCode2.c" ..
make
```
Paths to target sources are set as a semi-colon separated list of paths to said sources relative to the root QuEST directory.

If you wish your executable to be named something other than `demo`, you can set this too by using:
```bash
cmake -DOUTPUT_EXE="myExecutable" ..
```

When using the cmake command as above, the -D[VAR=VALUE] option can be passed other options to further configure your build.

You customise the precision with which the state-vector is stored.
```bash
cmake -DPRECISION=2 ..
```
Using greater precision means more precise computation but at the expense of additional memory requirements and runtime.
Checking results are unchanged when altaring the precision can be a great test that your calculations are sufficiently precise.

Please note that cmake caches these changes (per directory) so for any subsequent builds you should just type `make` from the build directory and the previously defined settings will be applied. If any parameters require changing, these can be redefined by:
```
cmake -D[VAR=VALUE] ..
```
as one would do in the initial configuration.

For a full list of available configuration parameters, use
```bash
cmake -LH ..
```

For manual configuration (not recommended) you can change the `CMakeLists.txt` in the root QuEST directory.

----------------------------

# Running

## locally

You can then call your code. From the build directory:
```bash
./myExecutable
```
If multithreading functionality was found when compiling, you can control how many threads your code uses by setting `OMP_NUM_THREADS`, ideally to the number of available cores on your machine
```bash
export OMP_NUM_THREADS=8
./myExecutable
```
QuEST will automatically allocate work between the given number of threads to speedup your simulation.

---------------------------

# Test circuits

There are several simple test circuits located in examples/timeSliceTNCircuits. To use one of these circuits, copy the circuit file to the root directory and then use

```bash
mkdir build
cd build
cmake -DUSER_SOURCE="testCircuitN.c"
make 
./demo
```

The list of test circuits is as follows:

testCircuit1.c
```
0 ----- * -----
        |
1 ----- o -----
```

testCircuit2.c: 
```
0 ----- * ----- * ----- 
        |       |       
1 ----- | ----- o ----- 
        |               
2 ----- o ------------- 
```

testCircuit3.c: 
```
0 ----- * ----- * ----- * -----
        |       |       |
1 ----- | ----- o ----- | -----
        |               |
2 ----- o ------------- o -----
```

