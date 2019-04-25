#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>
#include "QuEST_debug.h"

#include "tn.h"

void contractTensorNetwork(TensorNetwork tn, QuESTEnv env){
    // The simplest strategy -- contract tensor 0 with every other tensor in turn
    // Output will be stored as tensor 0
    // This is definitely not efficient in memory or time
    for (int i=1; i<tn.numTensors; i++){
        contractTensors(tn, 0, i, env);
    }
}

int getNumContractions(TensorNetwork tn, int tensor1Index, int tensor2Index){
    int numContractions = 0;
    VqVertex *tensor1Vertex, *entangledPair;
    tensor1Vertex = tn.tensorHeadVqVertex[tensor1Index];
    while (tensor1Vertex != NULL){
        entangledPair = tensor1Vertex->entangledPair;  
        if (entangledPair->tensorIndex == tensor2Index) numContractions++;
        tensor1Vertex = tensor1Vertex->nextInTensor;
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
    *numTensor1UnusedContractions = tn.numEntanglements[tensor1Index] - numContractions;
    *numTensor2UnusedContractions = tn.numEntanglements[tensor2Index] - numContractions;

    *tensor1UnusedContractions = malloc(*numTensor1UnusedContractions *
            sizeof(int));
    *tensor2UnusedContractions = malloc(*numTensor1UnusedContractions *
            sizeof(int));

    int contractionCount=0;
    int unusedCount=0;
    int vertexCount=0;

    VqVertex *tensorVertex, *entangledPair;
    tensorVertex = tn.tensorHeadVqVertex[tensor1Index];
    while (tensorVertex != NULL){
        entangledPair = tensorVertex->entangledPair;  
        if (entangledPair->tensorIndex == tensor2Index) {
            (*tensor1Contractions)[contractionCount++]=vertexCount++;
        } else {
            (*tensor1UnusedContractions)[unusedCount++]=vertexCount++;
        }
        tensorVertex = tensorVertex->nextInTensor;
    }

    contractionCount=0;
    unusedCount=0;
    vertexCount=0;
    tensorVertex = tn.tensorHeadVqVertex[tensor2Index];
    while (tensorVertex != NULL){
        entangledPair = tensorVertex->entangledPair;  
        if (entangledPair->tensorIndex == tensor1Index) {
            (*tensor2Contractions)[contractionCount++]=vertexCount++;
        } else {
            (*tensor2UnusedContractions)[unusedCount++]=vertexCount++;
        }
        tensorVertex = tensorVertex->nextInTensor;
    }
    printf("tensor1contracitons in function: %d\n", *tensor1Contractions[0]);
}

long long int getContractedIndex(long long int activeStateVec1Index, long long int activeStateVec2Index,
     long long int stateVec1PhysicalSize, long long int stateVec2PhysicalSize, int numTensor1Pq, int numTensor2Pq,
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

    printf("getOutputIndex: t1psize %lld t2psize %lld\n", stateVec1PhysicalSize, stateVec2PhysicalSize);
    contractedIndex = stateVec2IndexInPq*stateVec2PhysicalSize + stateVec1IndexInPq;
    printf("getOuputIndex: contractedIndex %lld\n", contractedIndex);

    // TODO -- missing unused tensor 2 virtual qubits
    // TODO -- seems inefficient given the only thing changing between this and
    // activeStateVectorIndex is number of Pq (if tensor2 has no unused vq?)
    activeIndexInBlock = activeStateVec1Index;
    for (int i=numTensor1UnusedContractions-1; i>=0; i--){
        // the size of a block of elements in the state vector for an unused virtual qubit
        sizeOutputHalfBlock = 1LL << (numTensor1Pq + numTensor2Pq + i);
        sizeActiveHalfBlock = 1LL << (numTensor1Pq + i);
        activeIsLower = (activeIndexInBlock >= sizeActiveHalfBlock); // TODO -- better to get value of qubit?
        printf("sizeOutput %lld sizeActive %lld isLower %d\n", sizeOutputHalfBlock, sizeActiveHalfBlock, activeIsLower);
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

    printf("recursive step. vqIndex: %d t1vqToContract %d\n", vqIndex, tensor1Contractions[vqIndex]);
    //TODO -- numPq will need to be numPq + numUnusedLowVq
    long long int tensor1OffsetNew = tensor1Offset + (1LL << (tensor1.numPq+tensor1Contractions[vqIndex]) );
    long long int tensor2OffsetNew = tensor2Offset + (1LL << (tensor2.numPq+tensor2Contractions[vqIndex]) );

    printf("tensor1Offset %lld tensor2Offset %lld\n", tensor1Offset, tensor2Offset);
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
    Tensor outputTensor;

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
                    stateVec1PhysicalSize, stateVec2PhysicalSize, tensor1.numPq, tensor2.numPq,
                    numTensor1UnusedContractions, numTensor2UnusedContractions);
            printf("OUTPUT: %lld VALUE: %f\n", contractedIndex, sum.real);
            contractedQureg.stateVec.real[contractedIndex] = sum.real;
            contractedQureg.stateVec.imag[contractedIndex] = sum.imag;
        }
    }

    // Update vqVertices in tn lists
    // first tensor
    VqVertex *tail=NULL;
    int foundHead=0;
    if (numTensor1UnusedContractions > 0){
        removeContractedVqVertices(tn, tensor1Index, tn.tensorHeadVqVertex[tensor1Index],
                tensor1UnusedContractions, numTensor1UnusedContractions, &tail, &foundHead);
    }
    if (numTensor2UnusedContractions > 0){
        tail->nextInTensor = tn.tensorHeadVqVertex[tensor2Index];
        removeContractedVqVertices(tn, tensor2Index, tail,
                tensor2UnusedContractions, numTensor2UnusedContractions, &tail, &foundHead);
    } 
    tn.tensorTailVqVertex[tensor1Index] = tail;
    if (tail != NULL) tail->nextInTensor = NULL;
   
    // second tensor 
    removeAllVqVertices(tn, tensor2Index);


    // Update first tensor 
    // The output of the contraction is stored in the location of the first tensor 
    // in the tensor array.
    destroyQureg(tensor1.qureg, env);
    outputTensor.qureg = contractedQureg;
    outputTensor.numPq = totalNumPq;
    outputTensor.numVq = totalNumVq;
    outputTensor.nextVqIndex = totalNumVq;
    tn.numEntanglements[tensor1Index] = totalNumVq;
    tn.tensors[tensor1Index] = outputTensor;
    
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
    tn.tensorHeadVqVertex = malloc( numTensors*sizeof(*(tn.tensorHeadVqVertex)) );
    tn.tensorTailVqVertex = malloc( numTensors*sizeof(*(tn.tensorTailVqVertex)) );
    tn.numEntanglements = malloc( numTensors*sizeof(*(tn.numEntanglements)) );

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
        //if (i==0)statevec_initStateDebugFromOffset(tmpTensor.qureg, 0);
        //if (i==1) hadamard(tmpTensor.qureg, 0);
        //if (i==1) rotateZ(tmpTensor.qureg, 0, 0.3);
        //if (i==0) hadamard(tmpTensor.qureg, 0);
        if (i==0) pauliX(tmpTensor.qureg, 0);
        if (i==0) rotateZ(tmpTensor.qureg, 0, 0.3);
        //if (i==1)statevec_initStateDebugFromOffset(tmpTensor.qureg, 2);

        // for quick conversion to local index from global index
        tmpTensor.firstGlobalPqIndex = globalPqOffset;

        // the next virtual qubit to use
        tmpTensor.nextVqIndex = 0;

        tn.tensors[i] = tmpTensor;
        tn.tensorHeadVqVertex[i] = NULL;
        tn.tensorTailVqVertex[i] = NULL;
        tn.numEntanglements[i] = 0;

        globalPqOffset += tmpTensor.numPq;
    }

    return tn;
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
    int globalPqOffset = 0;
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        tensor.firstGlobalPqIndex = globalPqOffset;
        globalPqOffset += tensor.numPq;
    }
}

void removeAllVqVertices(TensorNetwork tn, int tensorIndex){
    VqVertex *prevVqVertex;
    VqVertex *vqVertex = tn.tensorHeadVqVertex[tensorIndex];
    while (vqVertex != NULL){
        prevVqVertex = vqVertex;
        vqVertex = vqVertex->nextInTensor;
        free(prevVqVertex);
    }
}

void removeContractedVqVertices(TensorNetwork tn, int tensorIndex, VqVertex *startingVqVertex,
        int *unusedContractions, int numUnusedContractions, VqVertex **tail, int *foundHead){

    VqVertex *vqVertex, *prevVqVertex, *prevVqVertexToKeep=NULL;
    int vqVertexCount = 0;
    int unusedContractionCount = 0;
    vqVertex = startingVqVertex;
    // Search through vertices and keep those which are in unusedContractions.
    // Free vertices which no longer need to be stored
    while (vqVertex != NULL){
        if (unusedContractionCount<numUnusedContractions) {
            if (vqVertexCount++ == unusedContractions[unusedContractionCount]) {
                // keep 
                if (!(*foundHead)) {
                    tn.tensorHeadVqVertex[tensorIndex] = vqVertex;
                    *foundHead = 1;
                }
                if (prevVqVertexToKeep!=NULL){
                    prevVqVertexToKeep->nextInTensor = vqVertex;
                }
                prevVqVertexToKeep = vqVertex;
                vqVertex = vqVertex->nextInTensor;
                unusedContractionCount++;
            } else {
                // discard
                prevVqVertex = vqVertex;
                vqVertex = vqVertex->nextInTensor;
                free(prevVqVertex);
            }
        } else {
            // discard
            prevVqVertex = vqVertex;
            vqVertex = vqVertex->nextInTensor;
            free(prevVqVertex);
        }
    }

    *tail = prevVqVertexToKeep;
    prevVqVertexToKeep->nextInTensor=NULL;
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
     
        // TODO -- change t1, t2 to control/target 
        // update adjacency list. 
        // tensor1
        VqVertex *tensor1VqVertex = malloc(sizeof(*tensor1VqVertex));
        tensor1VqVertex->nextInTensor = NULL;
        tensor1VqVertex->tensorIndex = targetPqLocal.tensorIndex;
        if (tn.tensorTailVqVertex[targetPqLocal.tensorIndex] != NULL){
                tn.tensorTailVqVertex[targetPqLocal.tensorIndex]->nextInTensor = tensor1VqVertex;
        }
        tn.tensorTailVqVertex[targetPqLocal.tensorIndex] = tensor1VqVertex;
        if (tn.numEntanglements[targetPqLocal.tensorIndex]==0) {
            tn.tensorHeadVqVertex[targetPqLocal.tensorIndex] = tensor1VqVertex;
        }

        // tensor2
        VqVertex *tensor2VqVertex = malloc(sizeof(*tensor2VqVertex));
        tensor2VqVertex->nextInTensor = NULL;
        tensor2VqVertex->tensorIndex = controlPqLocal.tensorIndex;
        if (tn.tensorTailVqVertex[controlPqLocal.tensorIndex] != NULL){
                tn.tensorTailVqVertex[controlPqLocal.tensorIndex]->nextInTensor = tensor2VqVertex;
        }
        tn.tensorTailVqVertex[controlPqLocal.tensorIndex] = tensor2VqVertex;
        if (tn.numEntanglements[controlPqLocal.tensorIndex]==0) {
            tn.tensorHeadVqVertex[controlPqLocal.tensorIndex] = tensor2VqVertex;
        }

        tensor1VqVertex->entangledPair = tensor2VqVertex;
        tensor2VqVertex->entangledPair = tensor1VqVertex;

        tn.numEntanglements[controlPqLocal.tensorIndex]++;
        tn.numEntanglements[targetPqLocal.tensorIndex]++;
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

int getVqVertexIndex(TensorNetwork tn, VqVertex *vqVertex){
    int verticesUntilEndInclusive=0;
    int totalVertices = tn.numEntanglements[vqVertex->tensorIndex];
    if (vqVertex == NULL) return -1;
    while (vqVertex != NULL){
        verticesUntilEndInclusive++;
        vqVertex = vqVertex->nextInTensor;
    }
    return totalVertices - verticesUntilEndInclusive;
}

void printTensorNetwork(TensorNetwork tn){
    printf("\n----- TensorNetwork -------------------- \n");
    printf("%d tensors\n", tn.numTensors);
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        if (tensor.numPq==0){
            printf("Tensor %d: DELETED BY CONTRACTION\n", i);
        } else {
            printf("Tensor %d: %d physical qubits, %d virtual qubits\n", i, tensor.numPq, tensor.numVq);
            printf("\t%d Entanglements:\n", tn.numEntanglements[i]);
            
            VqVertex *vqVertex = tn.tensorHeadVqVertex[i];
            VqVertex *entangledPair;
            int vqCount=0;
            while (vqVertex != NULL){
                entangledPair = vqVertex->entangledPair;  
                int entangledIndex = getVqVertexIndex(tn, entangledPair);
                printf("\t\tVirtual qubit %d -> tensor %d, virtual qubit %d\n", 
                        vqCount++, entangledPair->tensorIndex, entangledIndex);
                vqVertex = vqVertex->nextInTensor;
            }
        }
    }
    printf("\n");
//    printf("\n--------------------------------------- \n");
}

















