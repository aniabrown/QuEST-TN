#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>
#include "QuEST_debug.h"

#include "tn.h"

void contractTensorNetwork(TensorNetwork tn){
}

int getNumContractions(TensorNetwork tn, int tensor1Index, int tensor2Index){
    int numContractions = 0;
    // adjacencies are stored on both tensors so just search tensor1
    for (int i=0; i<tn.numAdjacencies[tensor1Index]; i++){
        QCoord adjacency = tn.adjacencyList[tensor1Index][i];
        if (adjacency.tensorIndex == tensor2Index) numContractions++;
    }
    return numContractions;
}

void getContractionVqIndices(TensorNetwork tn, int tensor1Index, int tensor2Index, int numContractions,
        int **tensor1Contractions, int **tensor2Contractions,
        int *numTensor1UnusedContractions, int *numTensor2UnusedContractions,
        int **tensor1UnusedContractions, int **tensor2UnusedContractions){

    // malloc space for array of virtual qubits to contract
    *tensor1Contractions = malloc(numContractions*sizeof(int));
    *tensor2Contractions = malloc(numContractions*sizeof(int));
    
    // malloc space for array of virtual qubits on each tensor that will not be contracted
    *numTensor1UnusedContractions = tn.numAdjacencies[tensor1Index] - numContractions;
    *numTensor2UnusedContractions = tn.numAdjacencies[tensor2Index] - numContractions;

    *tensor1UnusedContractions = malloc(*numTensor1UnusedContractions *
            sizeof(int));
    *tensor2UnusedContractions = malloc(*numTensor1UnusedContractions *
            sizeof(int));

    int count=0;
    int unusedCount=0;

    for (int i=0; i<tn.numAdjacencies[tensor1Index]; i++){
        QCoord adjacency = tn.adjacencyList[tensor1Index][i];
        if (adjacency.tensorIndex == tensor2Index) {
            *tensor1Contractions[count++]=i;
        } else {
            *tensor1UnusedContractions[unusedCount++]=i;
        }
    }
    count=0;
    unusedCount=0;
    for (int i=0; i<tn.numAdjacencies[tensor2Index]; i++){
        QCoord adjacency = tn.adjacencyList[tensor2Index][i];
        if (adjacency.tensorIndex == tensor1Index) {
            *tensor2Contractions[count++]=i;
        } else {
            *tensor2UnusedContractions[unusedCount++]=i;
        }
    }
    printf("tensor1contracitons in function: %d\n", *tensor1Contractions[0]);
}

long long int getContractedIndex(long long int activeStateVec1Index, long long int activeStateVec2Index,
     long long int stateVec1PhysicalSize, long long int stateVec2PhysicalSize, int numTensor1Pq, int numTensor2Pq,
     int *tensor1UnusedContractions, int *tensor2UnusedContractions, 
     int numTensor1UnusedContractions, int numTensor2UnusedContractions){
            
    long long int contractedIndex;
    long long int stateVec1IndexInPq, stateVec2IndexInPq;

    long long int sizeOutputHalfBlock, sizeActiveHalfBlock;
    int activeIsLower;
    long long int stateVectorIndex = 0;
    long long int activeIndexInBlock;

    // activeStateVec1Index % sateVec1PhysicalSize = activeStateVec1Index&(stateVec1PhysicalSize-1)
    stateVec1IndexInPq = activeStateVec1Index&(stateVec1PhysicalSize-1);
    stateVec2IndexInPq = activeStateVec2Index&(stateVec2PhysicalSize-1);

    contractedIndex = stateVec2IndexInPq*stateVec2PhysicalSize + stateVec1IndexInPq;

    // TODO -- missing unused tensor 2 virtual qubits
    // TODO -- seems inefficient given the only thing changing between this and
    // activeStateVectorIndex is number of Pq (if tensor2 has no unused vq?)
    activeIndexInBlock = activeStateVec1Index;
    for (int i=numTensor1UnusedContractions-1; i>=0; i--){
        int vqIndex = tensor1UnusedContractions[i];
        // the size of a block of elements in the state vector for an unused virtual qubit
        sizeOutputHalfBlock = 1LL << (numTensor1Pq + numTensor2Pq + vqIndex);
        sizeActiveHalfBlock = 1LL << (numTensor1Pq + i);
        activeIsLower = (activeIndexInBlock >= sizeActiveHalfBlock); // TODO -- better to get value of qubit?
        if (activeIsLower) {
            contractedIndex += sizeOutputHalfBlock; // TODO -- replace with predicate
            activeIndexInBlock -= sizeActiveHalfBlock;
        }
    }

    return contractedIndex;
}

long long int getStateVectorIndexFromActiveIndex(long long int activeStateIndex,
        int numPq, int *unusedContractions, int numUnusedContractions){

    long long int sizeHalfBlock, sizeActiveHalfBlock;
    int activeIsLower;
    long long int activeStatePhysicalSize = 1LL << numPq;
    long long int indexInPhysicalQubitBlock = activeStateIndex&(activeStatePhysicalSize-1); 
    long long int stateVectorIndex = indexInPhysicalQubitBlock;
    long long int activeIndexInBlock = activeStateIndex;
    for (int i=numUnusedContractions-1; i>=0; i--){
        int vqIndex = unusedContractions[i];
        // the size of a block of elements in the state vector for an unused virtual qubit
        sizeHalfBlock = 1LL << (numPq + vqIndex);
        sizeActiveHalfBlock = 1LL << (numPq + i);
        activeIsLower = (activeIndexInBlock >= sizeActiveHalfBlock); // TODO -- better to get value of qubit?
        if (activeIsLower) {
            stateVectorIndex += sizeHalfBlock; // TODO -- replace with predicate
            activeIndexInBlock -= sizeActiveHalfBlock;
        }
   }
    return stateVectorIndex;
}

Complex recursiveContract(Tensor tensor1, Tensor tensor2, long long int tensor1Offset, 
        long long int tensor2Offset, int *tensor1Contractions, int *tensor2Contractions, 
        int numContractions, int vqIndex){

    printf("recursive step. vqIndex: %d\n", vqIndex);
    //TODO -- numPq will need to be numPq + numUnusedLowVq
    long long int tensor1OffsetNew = tensor1Offset + (1LL << (tensor1.numPq+tensor1Contractions[vqIndex]) );
    long long int tensor2OffsetNew = tensor2Offset + (1LL << (tensor2.numPq+tensor2Contractions[vqIndex]) );

    printf("tensor1OffsetNew %lld tensor2OffsetNew %lld\n", tensor1OffsetNew, tensor2OffsetNew);
    printf("num pq %d contraction index %d\n", tensor1.numPq, tensor1Contractions[vqIndex]);
   
    if (vqIndex==numContractions-1){
        Qureg qureg1 = tensor1.qureg;
        Qureg qureg2 = tensor2.qureg;

        qreal *stateVec1Real = qureg1.stateVec.real;
        qreal *stateVec1Imag = qureg1.stateVec.imag;
        qreal *stateVec2Real = qureg2.stateVec.real;
        qreal *stateVec2Imag = qureg2.stateVec.imag;

        qreal sumReal=0;
        qreal sumImag=0;

        // vqIndex == 0 
        // Real component
        sumReal += stateVec1Real[tensor1Offset] * 
            stateVec2Real[tensor2Offset];
        sumReal += -stateVec1Imag[tensor1Offset] * 
            stateVec2Imag[tensor2Offset];

        // Imag component
        sumImag += stateVec1Imag[tensor1Offset] * 
            stateVec2Real[tensor2Offset];
        sumImag += stateVec1Real[tensor1Offset] * 
            stateVec2Imag[tensor2Offset];

        // vqIndex == 1 
        // Real component
        sumReal += stateVec1Real[tensor1OffsetNew] * 
            stateVec2Real[tensor2OffsetNew];
        sumReal += -stateVec1Imag[tensor1OffsetNew] * 
            stateVec2Imag[tensor2OffsetNew];

        // Imag component
        sumImag += stateVec1Imag[tensor1OffsetNew] * 
            stateVec2Real[tensor2OffsetNew];
        sumImag += stateVec1Real[tensor1OffsetNew] * 
            stateVec2Imag[tensor2OffsetNew];

       
        Complex sum;
        sum.real = sumReal; sum.imag = sumImag;
        return sum;
    } else {
        Complex sumPart1, sumPart2, sum;
            sumPart1 = recursiveContract(tensor1, tensor2, tensor1Offset, tensor2Offset, 
                tensor1Contractions, tensor2Contractions, numContractions, vqIndex+1);
            sumPart2 = recursiveContract(tensor1, tensor2, tensor1OffsetNew, tensor2OffsetNew, 
                tensor1Contractions, tensor2Contractions, numContractions, vqIndex+1);
        sum.real = sumPart1.real + sumPart2.real;
        sum.imag = sumPart1.imag + sumPart2.imag;
        return sum;
    }
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

    int numContractions = getNumContractions(tn, tensor1Index, tensor2Index);
    printf("numContractions: %d\n", numContractions);

    int *tensor1Contractions, *tensor2Contractions; 
    int *tensor1UnusedContractions, *tensor2UnusedContractions; 
    int numTensor1UnusedContractions=0, numTensor2UnusedContractions=0;

    getContractionVqIndices(tn, tensor1Index, tensor2Index, numContractions, 
            &tensor1Contractions, &tensor2Contractions,
            &numTensor1UnusedContractions, &numTensor2UnusedContractions,
            &tensor1UnusedContractions, &tensor2UnusedContractions);
    printf("tensor1Contractions %d\n", tensor1Contractions[0]);
    printf("numUnusedContractions t1: %d t2: %d\n", numTensor1UnusedContractions, numTensor2UnusedContractions);

    int totalNumPq = tensor1.numPq + tensor2.numPq;
    int totalNumVq = numTensor1UnusedContractions + numTensor2UnusedContractions;
    int totalNumQ = totalNumPq + totalNumVq;

    long long int contractedStateVecSize, stateVec1Size, stateVec2Size, activeStateVec1Size, activeStateVec2Size;
    long long int contractedIndex, stateVec1Index, stateVec2Index, activeStateVec1Index, activeStateVec2Index;
    long long int stateVec1PhysicalSize, stateVec2PhysicalSize;

    contractedStateVecSize = 1LL << totalNumQ;
    stateVec1Size = 1LL << tensor1.numPq + tensor1.numVq;
    stateVec2Size = 1LL << tensor2.numPq + tensor1.numVq;
    stateVec1PhysicalSize = 1LL << tensor1.numPq;
    stateVec2PhysicalSize = 1LL << tensor2.numPq;

    activeStateVec1Size = stateVec1Size >> numContractions;
    activeStateVec2Size = stateVec2Size >> numContractions;

    Qureg contractedQureg = createQureg(totalNumQ, env);
    Qureg qureg1 = tensor1.qureg;
    Qureg qureg2 = tensor2.qureg;

    qreal sumReal, sumImag;
    qreal *stateVec1Real = qureg1.stateVec.real;
    qreal *stateVec1Imag = qureg1.stateVec.imag;
    qreal *stateVec2Real = qureg2.stateVec.real;
    qreal *stateVec2Imag = qureg2.stateVec.imag;

    for (activeStateVec2Index=0; activeStateVec2Index<activeStateVec2Size; activeStateVec2Index++) {
        for (activeStateVec1Index=0; activeStateVec1Index<activeStateVec1Size; activeStateVec1Index++) {
            stateVec1Index = getStateVectorIndexFromActiveIndex(activeStateVec1Index, tensor1.numPq,
                    tensor1UnusedContractions, numTensor1UnusedContractions);
            stateVec2Index = getStateVectorIndexFromActiveIndex(activeStateVec2Index, tensor2.numPq,
                    tensor2UnusedContractions, numTensor2UnusedContractions);
            printf("-------------------------\n");
            printf("t1: ai %lld si %lld\n", activeStateVec1Index, stateVec1Index);
            printf("t2: ai %lld si %lld\n", activeStateVec2Index, stateVec2Index);

            Complex sum;

            sum = recursiveContract(tensor1, tensor2, stateVec1Index, stateVec2Index, 
                tensor1Contractions, tensor2Contractions, numContractions, 0);

            contractedIndex = getContractedIndex(activeStateVec1Index, activeStateVec2Index, 
                    stateVec2PhysicalSize, stateVec1PhysicalSize, tensor1.numPq, tensor2.numPq,
                    tensor1UnusedContractions, tensor2UnusedContractions,
                    numTensor1UnusedContractions, numTensor2UnusedContractions);
            printf("OUTPUT: %lld \n", contractedIndex);
            contractedQureg.stateVec.real[contractedIndex] = sum.real;
            contractedQureg.stateVec.imag[contractedIndex] = sum.imag;
        }
    }

    // TODO -- Update tn adjacencies 

    // Update first tensor 
    // The output of the contraction is stored in the location of the first tensor 
    // in the tensor array.
    destroyQureg(tensor1.qureg, env);
    tn.tensors[tensor1Index].qureg = contractedQureg;
    tn.tensors[tensor1Index].numPq = totalNumPq;
    tn.tensors[tensor1Index].numVq = totalNumVq;
    tn.tensors[tensor1Index].nextVqIndex = totalNumVq;
    tn.tensors[tensor1Index].numAdjacencies = totalNumVq;
    
    // Update second tensor
    // The order of tensors in the tensor array is not changed when two tensors are contracted
    // -- for now, just leave a gap where the second tensor was by setting numPq to 0.
    destroyQureg(tensor2.qureg, env);
    tn.tensors[tensor2Index].numPq = 0;
    tn.tensors[tensor2Index].numVq = 0;

    // Update tn mapping between global qubit indices and (tensor, local qubit index)
    remapTensorIndexFromGlobalPq(tn);
    remapFirstGlobalPqIndex(tn);
    
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

// Assumes each tensor has the correct numPq and updates tensorIndexFromGlobalPq to be correct
void remapTensorIndexFromGlobalPq(TensorNetwork tn){
    int currentPq = 0;
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        for (int j=0; j<tensor.numPq; j++){
            tn.tensorIndexFromGlobalPq[currentPq++] = i;
        }
    }

}

void remapFirstGlobalPqIndex(TensorNetwork tn){
    globalPqOffset = 0;
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        tensor.firstGlobalPqIndex = globalPqOffset;
        globalPqOffset += tensor.numPq;
    }
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
        printf("controlled not 1: %d %d\n", controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq);
        controlledNot(controlTensor->qureg, controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq);
        
        // do target tensor half
        // TODO: note this won't work if doing multiple different operations on the same tensor in parallel 
        Tensor *targetTensor = &(tn.tensors[targetPqLocal.tensorIndex]);
        int virtualControlIndex = targetTensor->nextVqIndex;
        initVirtualControl(*targetTensor, virtualControlIndex + targetTensor->numPq);
        targetTensor->nextVqIndex = targetTensor->nextVqIndex + 1;
        printf("controlled not 2: %d %d\n", virtualControlIndex + targetTensor->numPq, targetPqLocal.qIndex);
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
        adjacency.tensorIndex = controlPqLocal.tensorIndex;
        adjacency.qIndex = virtualTargetIndex; 
        tn.adjacencyList[targetPqLocal.tensorIndex][virtualControlIndex] = adjacency;
        tn.numAdjacencies[controlPqLocal.tensorIndex]++;
        tn.numAdjacencies[targetPqLocal.tensorIndex]++;
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
        if (tensor.numPq==0){
            printf("Tensor %d: DELETED BY CONTRACTION\n", i);
        } else {
            printf("Tensor %d: %d physical qubits, %d virtual qubits\n", i, tensor.numPq, tensor.numVq);
            printf("\t%d Adjacencies:\n", tn.numAdjacencies[i]);
            QCoord *adjacencies = tn.adjacencyList[i];
            for (int j=0; j<tn.numAdjacencies[i]; j++){
                printf("\t\tVirtual qubit %d -> tensor %d, virtual qubit %d\n", j, adjacencies[j].tensorIndex,
                        adjacencies[j].qIndex);
            }
        }
    }
    printf("\n");
//    printf("\n--------------------------------------- \n");
}

















