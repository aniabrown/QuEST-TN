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

env = createQuESTEnv()

printf("Create tensor1. It will be in the zero state by default\n");
tensor1 = createTensor(1, 1, env)
reportStateToScreen(tensor1.qureg, env, 0)

printf("Create tensor2. It will be in the zero state by default\n");
tensor2 = createTensor(1, 1, env)
reportStateToScreen(tensor2.qureg, env, 0)

printf("Apply pauliX to qubit 0, which is qubit zero in tensor1\n");
TN_singleQubitGate(pauliX, tensor1, 0)
reportStateToScreen(tensor1.qureg, env, 0)

printf("Apply a controlled not gate controlled by qubit 0, to qubit 1\n");
TN_controlledGateControlHalf(controlledNot, tensor1, 0, 1)
TN_controlledGateTargetHalf(controlledNot, tensor2, 1, 0)

printf("Tensor1:\n");
reportStateToScreen(tensor1.qureg, env, 0)
printf("Tensor2:\n");
reportStateToScreen(tensor2.qureg, env, 0)

printf("Contract tensor1 and tensor2 across their virtual qubit index\n");
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

reportStateToScreen(outputTensor, env, 0)





