# [QuEST](https://quest.qtechtheory.org)

## Introduction

This is a prototype version of a tensor network extension to the [Quantum Exact Simulation Toolkit](https://github.com/QuEST-Kit/QuEST), a high performance simulator of universal quantum circuits, state-vectors and density matrices.  

The goal of this extension is to reduce the memory cost of representing a many qubit system in cases where the system is not maximally entangled.  

WARNING: this code is currently in early development and the API could change without warning

## Quick Start

Copy or clone this branch of the repository to your machine. E.g. in the desired directory, enter
```bash
git clone https://github.com/aniabrown/QuEST-TN.git
git checkout -t remotes/origin/blasExample
cd QuEST-TN
```
at terminal. You can then compile the library using
```bash
mkdir build
cd build
cmake ..
make
```

Dependencies: BLAS and Python3.4 or later to be installed on your system. 

Run the example circuit in the root directory with
```bash
export PYTHONPATH=$PYTHONPATH:utilities
python3 ./testCircuit1.py
```

The program will print information about your execution environment. It will then apply the same operations, first on a 2 qubit QuEST Qureg object and then on two tensors each containing one qubit. 

See the [TN tutorial](examples/timeSliceTNCircuits/README.md) for an introduction to the tensor network extension or the [QuEST tutorial](examples/README.md) for an introduction to using QuEST.

## Documentation

View the API [here](https://aniabrown.github.io/QuEST-TN/QuEST__tn_8h.html)

## Licence

QuEST-TN is released under a [MIT Licence](https://github.com/aniabrown/QuEST-TN/blob/master/LICENCE.txt)


