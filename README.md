# [QuEST](https://quest.qtechtheory.org)

## Introduction

This is a prototype version of a tensor network extension to the [Quantum Exact Simulation Toolkit](https://github.com/QuEST-Kit/QuEST), a high performance simulator of universal quantum circuits, state-vectors and density matrices.  

The goal of this extension is to reduce the memory cost of representing a many qubit system in cases where the system is not maximally entangled.  

WARNING: this code is currently in early development and the API could change without warning

## Quick Start

Copy or clone this repository to your machine. E.g. in the desired directory, enter
```bash
git clone https://github.com/aniabrown/QuEST-TN.git
cd QuEST-TN
```
at terminal. You can then compile the [simplest test circuit](examples/timeSliceTNCircuits/testCircuit1.c) using
```bash
mkdir build
cd build
cmake ..
make
```
then run it with
```bash
./demo
```
and afterward, clean up with
```bash
make clean
````

or, to remove the build directory entirely, from the root directory
```bash
rm -r build
```

The program will print information about your execution environment. It will then apply the same operations, first on a 2 qubit QuEST Qureg object and then on two tensors each containing one qubit. 

See the [TN tutorial](examples/timeSliceTNCircuits/README.md) for a better introduction or the [QuEST tutorial](examples/README.md) for an introduction to using QuEST.

## Documentation

View the API [here](https://aniabrown.github.io/QuEST-TN/QuEST__tn_8h.html)

## Licence

QuEST-TN is released under a [MIT Licence](https://github.com/aniabrown/QuEST-TN/blob/master/LICENCE.txt)


