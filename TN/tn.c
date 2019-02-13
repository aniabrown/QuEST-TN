#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>

#include "tn.h"

QuESTEnv env;

void contractTensorNetwork(TensorNetwork tn){
}

void contractTensors(TensorNetwork tn, int tensor1, int tensor2){
}
 
TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor){
    TensorNetwork tn;
    tn.numTensors = numTensors;

    // allocate memory in tensor object
    tn.tensors = malloc( numTensors*sizeof(*(tn.tensors)) );
    tn.adjacencyList = malloc( numTensors*sizeof(*(tn.adjacencyList)) )

    // calculate total number of qubits in system for qubits given per tensor
    int numTotalPq = 0;
    for (int i=0; i<numTensors; i++){
        numTotalPq += numPqPerTensor[i];
    }
    tn.tensorIdFromGlobalPq = malloc( numTotalPq*sizeof(*(tn.tensorIdFromGlobalPq)) );

    // populate tensorIdFromGlobalPq array. This contains duplicated values but should be
    // fast to look up
    int currentPq = 0;
    for (int i=0; i<numTensors; i++){
        int numPqInCurrentTensor = numPqPerTensor[i];
        for (int j=0; j<numPqInCurrentTensor; j++){
            tn.tensorIdFromGlobalPq[currentPq++] = i;
        }
    }

    int globalPqOffset = 0;

    // populate tensor array
    Tensor tmpTensor;
    for (int i=0; i<numTensors; i++){
        tmpTensor.numPq = numPqPerTensor[i];
        tmpTensor.numVq = numVqPerTensor[i];
        tmpTensor.qureg = createQureg(tmpTensor.numPq + tmpTensor.numVq, env);
        tmpTensor.firstGlobalPqIndex = globalPqOffset;
        tmpTensor.nextVqIndex = 0;

        tn.tensors[i] = tmpTensor;

        globalPqOffset += tmpTensor.numPq
    }

    // TODO: consider populating adjacency list with null values

}
 
void addTensorToNetwork(int numPq, int numVq){

}


void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit){


}














