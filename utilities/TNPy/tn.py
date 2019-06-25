from ctypes import *
from QuESTPy.QuESTBase import init_QuESTLib
from TNPy.TNBase import init_TNLib
from QuESTPy.QuESTLibDir import defaultQuESTDir

# If set up
QuESTPath = "build/TN/QuEST/"
TNPath = "build/TN/"

init_QuESTLib(QuESTPath)
init_TNLib(TNPath)

#!/usr/bin/env python3
from QuESTPy.QuESTFunc import *
from TNPy.TNFunc import *

def TN_singleQubitGate(gate, tn, *args):
    pq = getLocalPq(tn, args[0]) 
    gate(tn.tensors[pq.tensorIndex].qureg, pq.qIndex, *args[1:]);

env = createQuESTEnv()

print("Create Tensor Network")
numPqPerTensor = [1, 1]
numVqPerTensor = [1, 1]

numPqPerTensor = (c_int*2)(*numPqPerTensor)
numVqPerTensor = (c_int*2)(*numVqPerTensor)

tn = createTensorNetwork(2, numPqPerTensor, numVqPerTensor, env)
printTensorNetwork(tn)

u = ComplexMatrix2()
u.r0c0=Complex(0,0)
u.r0c1=Complex(1,0)
u.r1c0=Complex(1,0)
u.r1c1=Complex(0,0)
#print(u)


reportStateToScreen(tn.tensors[0].qureg, env, 0)
TN_singleQubitGate(pauliX, tn, 1)
reportStateToScreen(tn.tensors[0].qureg, env, 0)

"""
tn_controlledNot(tn, 0, 1)
tn_unitary(tn, 0, u)

printTensorNetwork(tn)
reportStateToScreen(tn.tensors[0].qureg, env, 0)
reportStateToScreen(tn.tensors[1].qureg, env, 0)

contractTensorNetwork(tn, env)

printTensorNetwork(tn)
reportStateToScreen(tn.tensors[0].qureg, env, 0)
"""





