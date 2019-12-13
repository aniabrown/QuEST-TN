from QuESTPy.QuESTBase import init_QuESTLib
from TNPy.TNBase import init_TNLib

# If set up
QuESTPath = "build/TN/QuEST/"
TNPath = "build/TN/"

init_QuESTLib(QuESTPath)
init_TNLib(TNPath)

#!/usr/bin/env python3
from QuESTPy.QuESTFunc import *
from TNPy.TNFunc import *
from TNPy.TNAdditionalGates import *

print("---------------------------------")
print("---- Create QuEST environment----")
print("---------------------------------")
env = createQuESTEnv()
reportQuESTEnv(env)

print("\n\n--------------------------------------------------------------")
print("---- QuEST example of a two qubit circuit: X(0), CNOT(0,1)----")
print("--------------------------------------------------------------")

print("Create qureg. It will be in the zero state by default")
qureg = createQureg(2, env);
print("qureg:")
print(qureg)
print("Apply PauliX to qubit 0")
pauliX(qureg, 0);
print("qureg:")
print(qureg)
print("Apply controlled not controlled by qubit 0 to qubit 1")
controlledNot(qureg, 0, 1);
print("qureg:")
print(qureg)

print("\n\n-----------------------------------------------------------")
print("---- TN example of a two qubit circuit: X(0), CNOT(0,1)----")
print("-----------------------------------------------------------")

print("Create tensor1. It will be in the zero state by default");
tensor1 = createTensor(1, 1, env)
print("Tensor1:");
print(tensor1.qureg)

print("Create tensor2. It will be in the zero state by default");
tensor2 = createTensor(1, 1, env)
print("Tensor2:");
print(tensor2.qureg)

print("Apply pauliX to global qubit 0, which is qubit 0 in tensor1");
TN_singleQubitGate(pauliX, tensor1, 0)
print("Tensor1:");
print(tensor1.qureg)

print("Apply a controlled not gate controlled by global qubit 0, to global qubit 1");
TN_controlledGateControlHalf(controlledNot, tensor1, 0, 1)
TN_controlledGateTargetHalf(controlledNot, tensor2, 1, 0)

print("Tensor1:");
print(tensor1.qureg)
print("Tensor2:");
print(tensor2.qureg)

print("Contract tensor1 and tensor2 across their virtual qubit index");
tensor1Contractions = [1]
tensor2Contractions = [1]
tensor1FreeIndices = [0]
tensor2FreeIndices = [0]

numContractions = 1
numTensor1FreeIndices = 1
numTensor2FreeIndices = 1

# Cast to correct type for c
tensor1Contractions = (c_int*numContractions)(*tensor1Contractions)
tensor2Contractions = (c_int*numContractions)(*tensor2Contractions)
tensor1FreeIndices = (c_int*numTensor1FreeIndices)(*tensor1FreeIndices)
tensor2FreeIndices = (c_int*numTensor2FreeIndices)(*tensor2FreeIndices)

outputTensor = contractIndices(tensor1, tensor2, tensor1Contractions, tensor2Contractions, numContractions,
		tensor1FreeIndices, numTensor2FreeIndices, tensor2FreeIndices, numTensor2FreeIndices, env)

print("output tensor:")
print(outputTensor.qureg)





