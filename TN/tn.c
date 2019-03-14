#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>
#include "QuEST_debug.h"

#include "tn.h"

void contractTensorNetwork(TensorNetwork tn){
}

void contractTensors(TensorNetwork tn, int tensor1Index, int tensor2Index, QuESTEnv env){
    // always store new tensor in the smaller tensor index
    if (tensor2Index < tensor1Index){
        int tensorIndexTmp = tensor1Index;
        tensor1Index = tensor2Index;
        tensor2Index = tensorIndexTmp;
    }

    Tensor tensor1 = tn.tensors[tensor1Index];
    Tensor tensor2 = tn.tensors[tensor2Index]; 

    int totalNumPq = tensor1.numPq + tensor2.numPq;
    //int totalNumQ = getTotalContractedQubits(); //TODO
    int totalNumQ = totalNumPq; //TODO
    Qureg contractedQureg = createQureg(totalNumPq, env);
    Qureg qureg1 = tensor1.qureg;
    Qureg qureg2 = tensor2.qureg;

    //int numContractions = getNumContractions(); //TODO
    int numContractions = 1;
    //int *tensor1Contractions = getContractions(); //TODO [v0, v1]
    //int *tensor1Contractions; //TODO [v0, v1]
    //*tensor1Contractions = tensor1.numPq;
    //int *tensor2Contractions = getContractions(); //TODO [v1, v2]
    //int *tensor2Contractions; //TODO [v1, v2]
    //*tensor2Contractions = tensor2.numPq;

    long long int contractedStateVecSize, stateVec1Size, stateVec2Size;
    long long int contractedIndex, stateVec1Index, stateVec2Index;

    contractedStateVecSize = 1LL << totalNumQ;
    stateVec1Size = 1LL << tensor1.numPq;
    stateVec2Size = 1LL << tensor2.numPq;

    qreal sumReal, sumImag;

    for (stateVec1Index=0; stateVec1Index<stateVec1Size; stateVec1Index++) {
        for (stateVec2Index=0; stateVec2Index<stateVec2Size; stateVec2Index++) {
            sumReal = qureg1.stateVec.real[stateVec1Index] * 
                qureg2.stateVec.real[stateVec2Index];
            sumReal -= qureg1.stateVec.imag[stateVec1Index+stateVec1Size] * 
                qureg2.stateVec.imag[stateVec2Index+stateVec2Size];

            sumImag = qureg1.stateVec.imag[stateVec1Index] * 
                qureg2.stateVec.real[stateVec2Index];
            sumImag += qureg1.stateVec.real[stateVec1Index+stateVec1Size] * 
                qureg2.stateVec.imag[stateVec2Index+stateVec2Size];

            /*
            sumReal = qureg1.stateVec.real[stateVec1Index] * 
                qureg2.stateVec.real[stateVec2Index];
            sumReal += qureg1.stateVec.imag[stateVec1Index+stateVec1Size] * 
                qureg2.stateVec.imag[stateVec2Index+stateVec2Size];

            sumImag = -qureg1.stateVec.imag[stateVec1Index] * 
                qureg2.stateVec.real[stateVec2Index];
            sumImag += qureg1.stateVec.real[stateVec1Index+stateVec1Size] * 
                qureg2.stateVec.imag[stateVec2Index+stateVec2Size];
    */

            contractedIndex = stateVec2Index*stateVec1Size + stateVec1Index;
            contractedQureg.stateVec.real[contractedIndex] = sumReal;
            contractedQureg.stateVec.imag[contractedIndex] = sumImag;
        }
    }
    
    destroyQureg(tensor1.qureg, env);
    destroyQureg(tensor2.qureg, env);

    tn.tensors[tensor1Index].qureg = contractedQureg;
    tn.tensors[tensor1Index].numPq = totalNumPq;
    tn.tensors[tensor1Index].numVq = 0;

    //TODO: remove adjacencies
}
 
TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor, 
        QuESTEnv env){

    TensorNetwork tn;
    tn.numTensors = numTensors;

    // allocate memory in tensor object
    tn.tensors = malloc( numTensors*sizeof(*(tn.tensors)) );
    tn.adjacencyList = malloc( numTensors*sizeof(*(tn.adjacencyList)) );
    tn.numAdjacencies = malloc( numTensors*sizeof(*(tn.numAdjacencies)) );

    // calculate total number of qubits in system for qubits given per tensor
    int numTotalPq = 0;
    int numTotalVq = 0;
    for (int i=0; i<numTensors; i++){
        numTotalPq += numPqPerTensor[i];
        numTotalVq += numVqPerTensor[i];
    }
    tn.tensorIndexFromGlobalPq = malloc( numTotalPq*sizeof(*(tn.tensorIndexFromGlobalPq)) );

    // populate tensorIndexFromGlobalPq array. This contains duplicated values but should be
    // fast to look up
    int currentPq = 0;
    for (int i=0; i<numTensors; i++){
        int numPqInCurrentTensor = numPqPerTensor[i];
        for (int j=0; j<numPqInCurrentTensor; j++){
            tn.tensorIndexFromGlobalPq[currentPq++] = i;
        }
    }

    int globalPqOffset = 0;

    // populate tensor array
    Tensor tmpTensor;
    for (int i=0; i<numTensors; i++){
        tmpTensor.numPq = numPqPerTensor[i];
        tmpTensor.numVq = numVqPerTensor[i];

        // create qureg object and initialize to zero
        tmpTensor.qureg = createQureg(tmpTensor.numPq + tmpTensor.numVq, env);
        initZeroState(tmpTensor.qureg);
        //if (i==0) hadamard(tmpTensor.qureg, 1);
        if (i==0) pauliX(tmpTensor.qureg, 0);

        // for quick conversion to local index from global index
        tmpTensor.firstGlobalPqIndex = globalPqOffset;

        // the next virtual qubit to use
        tmpTensor.nextVqIndex = 0;

        tn.tensors[i] = tmpTensor;
        // The adjacency list per tensor will never be larger than numTotalVq so allocate
        // that much memory for now. 
        // For large numbers of virtual qubits may be more efficient to reallocate memory
        // each time a new element is added to this list instead.
        tn.adjacencyList[i] = malloc( numTotalVq*sizeof(*(tn.adjacencyList[i])) );
        tn.numAdjacencies[i] = 0;

        globalPqOffset += tmpTensor.numPq;
    }

    // TODO: consider populating adjacency list with null values
    
    return tn;
}
 
void addTensorToNetwork(int numPq, int numVq){

}

// ----- operations ------------------------------------------------------------


void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit){
    QCoord controlPqLocal = getLocalPq(tn, controlQubit);
    QCoord targetPqLocal = getLocalPq(tn, targetQubit);

    Tensor *controlTensor = &(tn.tensors[controlPqLocal.tensorIndex]);

    if (controlPqLocal.tensorIndex == targetPqLocal.tensorIndex){
        // qubits are in same tensor/qureg. Just do a normal qureg operation
        controlledNot(controlTensor->qureg, controlPqLocal.qIndex, targetPqLocal.qIndex);
    } else {
        // qubits are in different tensor/qureg
        // do control tensor half
        // TODO: note this won't work if doing multiple different operations on the same tensor in parallel 
        int virtualTargetIndex = controlTensor->nextVqIndex;
        // TODO -- shouldn't have to initialize virtual target when virtual qubits are
        // automatically initialized in the zero state
        initVirtualTarget(*controlTensor, virtualTargetIndex + controlTensor->numPq);
        controlTensor->nextVqIndex = controlTensor->nextVqIndex + 1;
        controlledNot(controlTensor->qureg, controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq);
        
        // do target tensor half
        // TODO: note this won't work if doing multiple different operations on the same tensor in parallel 
        Tensor *targetTensor = &(tn.tensors[targetPqLocal.tensorIndex]);
        int virtualControlIndex = targetTensor->nextVqIndex;
        initVirtualControl(*targetTensor, virtualControlIndex + targetTensor->numPq);
        targetTensor->nextVqIndex = targetTensor->nextVqIndex + 1;
        controlledNot(targetTensor->qureg, virtualControlIndex + targetTensor->numPq, targetPqLocal.qIndex);
      
        // update adjacency list. 
        // Only store the connection between virtualTargetIndex to virtualControlIndex, not the other direction
        // TODO: May want to store both diretions 
        // virtualTargetIndex -> virtualControlIndex corresponds to the tensor with the control pq pointing to
        // the tensor with the target pq
        QCoord adjacency;
        adjacency.tensorIndex = targetPqLocal.tensorIndex;
        adjacency.qIndex = virtualControlIndex; 
        tn.adjacencyList[controlPqLocal.tensorIndex][virtualTargetIndex] = adjacency;
        printf("target: %d\n", tn.adjacencyList[controlPqLocal.tensorIndex][virtualTargetIndex].tensorIndex);
        printf("%d %d\n", controlPqLocal.tensorIndex, virtualTargetIndex);
        tn.numAdjacencies[controlPqLocal.tensorIndex] = tn.numAdjacencies[controlPqLocal.tensorIndex] + 1;
    }

}

void initVirtualTarget(Tensor tensor, int vqIndex){
    // place target virtual qubit in the zero state
    
     collapseToOutcome(tensor.qureg, vqIndex, 0);
}

void initVirtualControl(Tensor tensor, int vqIndex){
    // initialize qubit vqIndex to |0> + |1> 
    // NOTE: not (1/sqrt(2))(|0> + |1>)

    // TODO: add API funtions to do this instead?

    Qureg qureg = tensor.qureg;

    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    // TODO: This won't work in parallel
    stateVecSize = 1LL << vqIndex;

    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, stateVecReal, stateVecImag) \
    private  (index)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            stateVecReal[index+stateVecSize] = stateVecReal[index];
            stateVecImag[index+stateVecSize] = stateVecImag[index];
        }
    }
}


// ----- utility --------------------------------------------------------------
QCoord getLocalPq(TensorNetwork tn, int globalPq){
    QCoord localPq;
    localPq.tensorIndex = tn.tensorIndexFromGlobalPq[globalPq];
    localPq.qIndex = globalPq - tn.tensors[localPq.tensorIndex].firstGlobalPqIndex;
    return localPq;
}

// ----- reporting ------------------------------------------------------------

void printTensorNetwork(TensorNetwork tn){
    printf("\n----- TensorNetwork -------------------- \n");
    printf("%d tensors\n", tn.numTensors);
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        printf("Tensor %d: %d physical qubits, %d virtual qubits\n", i, tensor.numPq, tensor.numVq);
        printf("\t%d Adjacencies:\n", tn.numAdjacencies[i]);
        QCoord *adjacencies = tn.adjacencyList[i];
        for (int j=0; j<tn.numAdjacencies[i]; j++){
            printf("\t\tVirtual qubit %d -> tensor %d, virtual qubit %d\n", j, adjacencies[j].tensorIndex,
                    adjacencies[j].qIndex);
        }
    }
    printf("\n--------------------------------------- \n");
}

















