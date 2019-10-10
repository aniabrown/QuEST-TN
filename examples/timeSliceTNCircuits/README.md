QuEST-TN Tutorial
======

**Table of Contents**
- [Background](#background)
- [Coding](#coding)

# Background

## Motivation

The fundamental idea behind the tensor network extension is to save memory and time by splitting the full system of multiple qubits, represented by single a Qureg object in QuEST, into multiple independent collections of qubits, each with its own Qureg. As an example, for a system of 30 qubits, 2^30 complex doubles need to be stored in the usual QuEST implementation. If this can be split into two systems of 15 qubits, only 2*2^15 elements need to be stored. Any gates applied to one of these systems will also only need to process 2^15 elements. 

It can be useful to treat these subsystems as tensors of rank numQubits, with each tensor index ranging from {0,1} and corresponding to a single qubit.

## Entanglement

Naturally, these independent data structures do not capture any entanglement between qubits in different subsystems. However, an entangling operation between qubits in different tensors can be represented by introducing an additional 'virtual' qubit to each of the two tensors, as described below. It is necessary to balance the savings due to splitting a system of qubits into multiple subsystems against the number of entangling operations between the two systems, each introducing an additional qubit. The method will be most useful for circuits with weakly entangled subsystems.

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

To use the QuEST-TN extension in Python code located in the root directory, include

```
from QuESTPy.QuESTBase import init_QuESTLib
from TNPy.TNBase import init_TNLib

QuESTPath = "build/TN/QuEST/"
TNPath = "build/TN/"

init_QuESTLib(QuESTPath)
init_TNLib(TNPath)

from QuESTPy.QuESTFunc import *
from TNPy.TNFunc import *
from TNPy.TNAdditionalGates import *
```

The currently available API for operations on tensor networks is at [QuEST_tn.h](TN/QuEST_tn.h) and the python-only functions are in [TNAdditionalGates.py](utilities/TNPy/TNAdditionalGates.py). For operations on plain QuEST Qureg objects, see the API [here](https://quest-kit.github.io/QuEST/QuEST_8h.html) and the [QuEST tutorial](examples/README.md). 

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


