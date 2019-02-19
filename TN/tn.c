#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>

#include "tn.h"

void contractTensorNetwork(TensorNetwork tn){
}

void contractTensors(TensorNetwork tn, int tensor1, int tensor2){
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
        controlTensor->nextVqIndex = controlTensor->nextVqIndex + 1;
        // shouldn't have to initialize virtual target as virtual qubits are currently 
        // automatically initialized in the zero state
        //initVirtualTarget(controlTensor, virtualTargetIndex + controlTensor->numPq);
        controlledNot(controlTensor->qureg, controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq);
        
        // do target tensor half
        // TODO: note this won't work if doing multiple different operations on the same tensor in parallel 
        Tensor *targetTensor = &(tn.tensors[targetPqLocal.tensorIndex]);
        int virtualControlIndex = targetTensor->nextVqIndex;
        targetTensor->nextVqIndex = targetTensor->nextVqIndex + 1;
        initVirtualControl(*targetTensor, virtualControlIndex + targetTensor->numPq);
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
    int numFilledQubits = tensor.numPq + tensor.nextVqIndex;
    stateVecSize = 1LL << numFilledQubits;

    // Can't use qureg->stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    // initialise the state-vector to all-zeroes
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

















