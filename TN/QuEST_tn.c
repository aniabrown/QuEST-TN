#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>
#include "QuEST_debug.h"
#include "QuEST_tn_internal.h"
#include "QuEST_tn.h"

#define DEBUG 0

#if DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif


void contractTensorNetwork(TensorNetwork tn, QuESTEnv env){
    for (int i=1; i<tn.numTensors; i++){
        contractTensors(tn, 0, i, env);
    }
}

/** Get the number of indices that will be contracted over when the two input
 * tensors are contracted. 
 * These correspond to all the virtual qubit indices that represent entanglement
 * between the two tensors.
 *
 * @returns the number of indices to contract
 * @param[in] tn the tensor network that the two tensors belong to
 * @param[in] tensor1Index the first tensor to contract. This should have the
 *      lower index of the two tensors.
 * @param[in] tensor2Index the second tensor to contract
 */
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

/** For a particular contraction of two tensors, generate the lists of virtual qubits that will
 * be contracted and that will remain uncontracted for each tensor.
 * The number of virtual qubits to contract must be the same for each tensor but the number
 * of uncontracted qubits may differ
 *
 * @param[in] tn the tensor network that the two tensors belong to
 * @param[in] tensor1Index the first tensor to contract. This should have the
 *      lower index of the two tensors.
 * @param[in] tensor2Index the second tensor to contract
 * @param[in] numContractions the number of indices to contract over
 * @param[out] tensor1Contractions tensor1 virtual qubits to contract
 * @param[out] tensor2Contractions tensor2 virtual qubits to contract
 * @param[out] numTensor1UncontractedVqs number of virtual qubits in tensor1 that represent
 *      entanglement with tensors other than tensor2
 * @param[out] numTensor2UncontractedVqs number of virtual qubits in tensor2 that represent
 *      entanglement with tensors other than tensor1
 * @param[out] tensor1UncontractedVqs virtual qubit indices that will remain uncontracted
 *      in tensor1
 * @param[out] tensor1UncontractedVqs virtual qubit indices that will remain uncontracted
 *      in tensor2
 */
void getContractionVqIndices(TensorNetwork tn, int tensor1Index, int tensor2Index, int numContractions,
        int **tensor1Contractions, int **tensor2Contractions,
        int *numTensor1UncontractedVqs, int *numTensor2UncontractedVqs,
        int **tensor1UncontractedVqs, int **tensor2UncontractedVqs){

    // malloc space for array of virtual qubits to contract
    *tensor1Contractions = malloc(numContractions*sizeof(int));
    *tensor2Contractions = malloc(numContractions*sizeof(int));
    
    // malloc space for array of virtual qubits on each tensor that will not be contracted
    *numTensor1UncontractedVqs = tn.numEntanglements[tensor1Index] - numContractions;
    *numTensor2UncontractedVqs = tn.numEntanglements[tensor2Index] - numContractions;

    *tensor1UncontractedVqs = malloc(*numTensor1UncontractedVqs *
            sizeof(int));
    *tensor2UncontractedVqs = malloc(*numTensor1UncontractedVqs *
            sizeof(int));

    int contractionCount=0;
    int uncontractedCount=0;
    int vertexCount=0;

    VqVertex *tensorVertex, *entangledPair;
    tensorVertex = tn.tensorHeadVqVertex[tensor1Index];
    while (tensorVertex != NULL){
        entangledPair = tensorVertex->entangledPair;  
        if (entangledPair->tensorIndex == tensor2Index) {
            (*tensor1Contractions)[contractionCount]=vertexCount++;
            // add indices to contract in second tensor in same order as for first tensor
            int entangledIndex = getVqVertexIndex(tn, entangledPair);
            (*tensor2Contractions)[contractionCount++]=entangledIndex;
        } else {
            (*tensor1UncontractedVqs)[uncontractedCount++]=vertexCount++;
        }
        tensorVertex = tensorVertex->nextInTensor;
    }

    uncontractedCount=0;
    vertexCount=0;
    tensorVertex = tn.tensorHeadVqVertex[tensor2Index];
    while (tensorVertex != NULL){
        entangledPair = tensorVertex->entangledPair;  
        if (!(entangledPair->tensorIndex == tensor1Index)) {
            (*tensor2UncontractedVqs)[uncontractedCount++]=vertexCount;
        }
        vertexCount++;
        tensorVertex = tensorVertex->nextInTensor;
    }
}

/** Get the index in the output state vector where an element resulting from a contraction between two
 * tensors should be placed.
 * 
 * @param[in] tensor1FreeIndexEl the index of the element in a state vector made up of only the free indices in tensor1 
 * @param[in] tensor2FreeIndexEl the index of the element in a state vector made up of only the free indices in tensor2 
 * @param[in] tensor1PqSize the number of probability amplitudes in a system containing just the physical qubits in tensor1 
 * @param[in] tensor2PqSize the number of probability amplitudes in a system containing just the physical qubits in tensor2 
 * @param[in] numTensor1Pq the number of physical qubits in tensor1
 * @param[in] numTensor2Pq the number of physical qubits in tensor2
 * @param[in] numTensor1UncontractedVqs the number of virtual qubits which will not be contracted in this operation
 * 	in tensor1 
 * @param[in] numTensor1UncontractedVqs the number of virtual qubits which will not be contracted in this operation
 * 	in tensor2 
 * @returns the index of the element in the output state vector 
 */
long long int getContractedIndex(long long int tensor1FreeIndexEl, long long int tensor2FreeIndexEl,
     long long int tensor1PqSize, long long int tensor2PqSize, int numTensor1Pq, int numTensor2Pq,
     int numTensor1UncontractedVqs, int numTensor2UncontractedVqs){
            
    long long int contractedIndex;
    long long int stateVec1IndexInPq, stateVec2IndexInPq;

    long long int sizeOutputHalfBlock, sizeFreeHalfBlock;
    int freeIsLower;
    long long int stateVectorIndex = 0;
    long long int freeIndexInBlock;

    // tensor1FreeIndexEl % sateVec1PhysicalSize = tensor1FreeIndexEl&(tensor1PqSize-1)
    stateVec1IndexInPq = tensor1FreeIndexEl&(tensor1PqSize-1);
    stateVec2IndexInPq = tensor2FreeIndexEl&(tensor2PqSize-1);

    DEBUG_PRINT(("getOutputIndex: t1psize %lld t2psize %lld\n", tensor1PqSize, tensor2PqSize));
    contractedIndex = stateVec2IndexInPq*tensor1PqSize + stateVec1IndexInPq;
    DEBUG_PRINT(("getOuputIndex: contractedIndex %lld\n", contractedIndex));

    // TODO -- seems inefficient given the only thing changing between this and
    // activeStateVectorIndex is number of Pq (if tensor2 has no uncontracted vq?)
    freeIndexInBlock = tensor1FreeIndexEl;
    for (int i=numTensor1UncontractedVqs-1; i>=0; i--){
        // the size of a block of elements in the state vector for an uncontracted virtual qubit
        sizeOutputHalfBlock = 1LL << (numTensor1Pq + numTensor2Pq + i);
        sizeFreeHalfBlock = 1LL << (numTensor1Pq + i);
        freeIsLower = (freeIndexInBlock >= sizeFreeHalfBlock); // TODO -- better to get value of qubit?
        DEBUG_PRINT(("sizeOutput %lld sizeActive %lld isLower %d\n", sizeOutputHalfBlock, sizeFreeHalfBlock, freeIsLower));
        if (freeIsLower) {
            contractedIndex += sizeOutputHalfBlock; // TODO -- replace with predicate
            freeIndexInBlock -= sizeFreeHalfBlock;
        }
    }
    
    freeIndexInBlock = tensor2FreeIndexEl;
    for (int i=numTensor2UncontractedVqs-1; i>=0; i--){
        // the size of a block of elements in the state vector for an uncontracted virtual qubit
        sizeOutputHalfBlock = 1LL << (numTensor1Pq + numTensor2Pq + numTensor1UncontractedVqs + i);
        sizeFreeHalfBlock = 1LL << (numTensor2Pq + i);
        freeIsLower = (freeIndexInBlock >= sizeFreeHalfBlock); // TODO -- better to get value of qubit?
        DEBUG_PRINT(("sizeOutput %lld sizeActive %lld isLower %d\n", sizeOutputHalfBlock, sizeFreeHalfBlock, freeIsLower));
        if (freeIsLower) {
            contractedIndex += sizeOutputHalfBlock; // TODO -- replace with predicate
            freeIndexInBlock -= sizeFreeHalfBlock;
        }
    }

    return contractedIndex;
}

/** Get the corresponding index in the full state vector for a tensor, given an element in a state vector made up of 
 * only the free indices in that tensor. The index in the full state vector will have all physical qubits and 
 * uncontracted virtual qubits with the same value as the input element. 
 * 
 * @param[in] freeIndexEl the index of the element in a state vector made up of only the free indices in a tensor
 * @param[in] numPq the number of physical qubits in that tensor
 * @param[in] uncontractedVqs a list of vq indices that will not be contracted
 * @param[in] the number of elements in uncontractedVqs
 * @returns the index of the element in the full state vector
*/
long long int getStateVectorIndexFromFreeIndexEl(long long int freeIndexEl,
        int numPq, int *uncontractedVqs, int numUncontractedVqs){

    long long int sizeHalfBlock, sizeFreeHalfBlock;
    int freeIsLower;
    long long int pqSize = 1LL << numPq;
    long long int indexInPhysicalQubitBlock = freeIndexEl&(pqSize-1); 
    long long int stateVectorIndex = indexInPhysicalQubitBlock;
    long long int freeIndexInBlock = freeIndexEl;
    for (int i=numUncontractedVqs-1; i>=0; i--){
        int vqIndex = uncontractedVqs[i];
        // the size of a block of elements in the state vector for an uncontracted virtual qubit
        sizeHalfBlock = 1LL << (numPq + vqIndex);
        sizeFreeHalfBlock = 1LL << (numPq + i);
        freeIsLower = (freeIndexInBlock >= sizeFreeHalfBlock); // TODO -- better to get value of qubit?
        if (freeIsLower) {
            stateVectorIndex += sizeHalfBlock; // TODO -- replace with predicate
            freeIndexInBlock -= sizeFreeHalfBlock;
        }
   }
    return stateVectorIndex;
}

/** For a particular element in tensor1 and a particular element in tensor2, recursively contract across all virtual qubits 
 * that represent entanglements between the two tensors. The elements in each tensor will always represent probability
 * amplitudes where all virtual qubits to contract are equal to zero.
 * 
 * @param[in] tensor1 the first tensor to contract. 
 * @param[in] tensor2 the second tensor to contract.
 * @param[in] tensor1Offset an element in the tensor1 state vector where all qubits to contract are equal to zero
 * @param[in] tensor2Offset an element in the tensor2 state vector where all qubits to contract are equal to zero
 * @param[in] tensor1Contractions a list of virtual qubit indices to contract with the tensor2 indices
 * @param[in] tensor2Contractions a list of virtual qubit indices to contract with the tensor1 indices 
 * @param[in] numContractions the number of indices to contract
 * @param[in] the current index in the tensor1Contractions/tensor2Contractions list 
 * @returns one element in the output tensor after contraction
 */
Complex recursiveContract(Tensor tensor1, Tensor tensor2, long long int tensor1Offset, 
        long long int tensor2Offset, int *tensor1Contractions, int *tensor2Contractions, 
        int numContractions, int vqIndex){

    DEBUG_PRINT(("recursive step. vqIndex: %d t1vqToContract %d\n", vqIndex, tensor1Contractions[vqIndex]));
    long long int tensor1OffsetNew = tensor1Offset + (1LL << (tensor1.numPq+tensor1Contractions[vqIndex]) );
    long long int tensor2OffsetNew = tensor2Offset + (1LL << (tensor2.numPq+tensor2Contractions[vqIndex]) );

    DEBUG_PRINT(("tensor1Offset %lld tensor2Offset %lld\n", tensor1Offset, tensor2Offset));
    DEBUG_PRINT(("tensor1OffsetNew %lld tensor2OffsetNew %lld\n", tensor1OffsetNew, tensor2OffsetNew));
    DEBUG_PRINT(("num pq %d contraction index %d\n", tensor1.numPq, tensor1Contractions[vqIndex]));
   
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
    DEBUG_PRINT(("numContractions: %d\n", numContractions));

    int *tensor1Contractions, *tensor2Contractions; 
    int *tensor1UncontractedVqs, *tensor2UncontractedVqs; 
    int numTensor1UncontractedVqs=0, numTensor2UncontractedVqs=0;

    getContractionVqIndices(tn, tensor1Index, tensor2Index, numContractions, 
            &tensor1Contractions, &tensor2Contractions,
            &numTensor1UncontractedVqs, &numTensor2UncontractedVqs,
            &tensor1UncontractedVqs, &tensor2UncontractedVqs);

    DEBUG_PRINT(("tensor1Contractions: "));
    for (int p=0; p<numContractions; p++) {DEBUG_PRINT(("%d ", tensor1Contractions[p])); }
    DEBUG_PRINT(("\ntensor2Contractions: "));
    for (int p=0; p<numContractions; p++) {DEBUG_PRINT(("%d ", tensor2Contractions[p])); }
    DEBUG_PRINT(("\ntensor1UncontractedVqs: "));
    for (int p=0; p<numTensor1UncontractedVqs; p++) {DEBUG_PRINT(("%d ", tensor1UncontractedVqs[p])); }
    DEBUG_PRINT(("\ntensor2UncontractedVqs: "));
    for (int p=0; p<numTensor2UncontractedVqs; p++) {DEBUG_PRINT(("%d ", tensor2UncontractedVqs[p])); }
    DEBUG_PRINT(("\n"));

    int totalNumPq = tensor1.numPq + tensor2.numPq;
    int totalNumVq = numTensor1UncontractedVqs + numTensor2UncontractedVqs;
    int totalNumQ = totalNumPq + totalNumVq;

    long long int contractedStateVecSize, stateVec1Size, stateVec2Size, tensor1FreeIndexSize, tensor2FreeIndexSize;
    long long int contractedIndex, stateVec1Index, stateVec2Index, tensor1FreeIndexEl, tensor2FreeIndexEl;
    long long int tensor1PqSize, tensor2PqSize;

    contractedStateVecSize = 1LL << totalNumQ;
    stateVec1Size = 1LL << (tensor1.numPq + tensor1.numVq);
    stateVec2Size = 1LL << (tensor2.numPq + tensor2.numVq);
    tensor1PqSize = 1LL << tensor1.numPq;
    tensor2PqSize = 1LL << tensor2.numPq;

    tensor1FreeIndexSize = stateVec1Size >> numContractions;
    tensor2FreeIndexSize = stateVec2Size >> numContractions;

    Qureg contractedQureg = createQureg(totalNumQ, env);
    Qureg qureg1 = tensor1.qureg;
    Qureg qureg2 = tensor2.qureg;

    qreal sumReal, sumImag;
    qreal *stateVec1Real = qureg1.stateVec.real;
    qreal *stateVec1Imag = qureg1.stateVec.imag;
    qreal *stateVec2Real = qureg2.stateVec.real;
    qreal *stateVec2Imag = qureg2.stateVec.imag;

    DEBUG_PRINT(("tensor1FreeIndexSize %lld, tensor2FreeIndexSize %lld\n", tensor1FreeIndexSize, tensor2FreeIndexSize));

    for (tensor2FreeIndexEl=0; tensor2FreeIndexEl<tensor2FreeIndexSize; tensor2FreeIndexEl++) {
        for (tensor1FreeIndexEl=0; tensor1FreeIndexEl<tensor1FreeIndexSize; tensor1FreeIndexEl++) {
            stateVec1Index = getStateVectorIndexFromFreeIndexEl(tensor1FreeIndexEl, tensor1.numPq,
                    tensor1UncontractedVqs, numTensor1UncontractedVqs);
            stateVec2Index = getStateVectorIndexFromFreeIndexEl(tensor2FreeIndexEl, tensor2.numPq,
                    tensor2UncontractedVqs, numTensor2UncontractedVqs);
            DEBUG_PRINT(("-------------------------\n"));
            DEBUG_PRINT(("t1: ai %lld si %lld\n", tensor1FreeIndexEl, stateVec1Index));
            DEBUG_PRINT(("t2: ai %lld si %lld\n", tensor2FreeIndexEl, stateVec2Index));

            Complex sum;

            sum = recursiveContract(tensor1, tensor2, stateVec1Index, stateVec2Index, 
                tensor1Contractions, tensor2Contractions, numContractions, 0);

            contractedIndex = getContractedIndex(tensor1FreeIndexEl, tensor2FreeIndexEl, 
                    tensor1PqSize, tensor2PqSize, tensor1.numPq, tensor2.numPq,
                    numTensor1UncontractedVqs, numTensor2UncontractedVqs);
            DEBUG_PRINT(("OUTPUT: %lld VALUE: %f\n", contractedIndex, sum.real));
            contractedQureg.stateVec.real[contractedIndex] = sum.real;
            contractedQureg.stateVec.imag[contractedIndex] = sum.imag;
        }
    }

    // Update vqVertices in tn lists
    
    // first tensor's vqVertices
    VqVertex *tail=NULL;
    int foundHead=0;
    removeContractedVqVertices(tn, tensor1Index, tensor1Index,
            NULL, tensor1UncontractedVqs, numTensor1UncontractedVqs, &tail, &foundHead);

    // second tensor's vqVertices
    if (foundHead){
        // there was at least one uncontracted index in the first tensor
        tail->nextInTensor = tn.tensorHeadVqVertex[tensor2Index];
        removeContractedVqVertices(tn, tensor2Index, tensor1Index, tail,
                tensor2UncontractedVqs, numTensor2UncontractedVqs, &tail, &foundHead);
    } else {
        // there were no uncontracted indices in the first tensor
        removeContractedVqVertices(tn, tensor2Index, tensor1Index, 
                NULL, tensor2UncontractedVqs, numTensor2UncontractedVqs, &tail, &foundHead);
    }

    if (!foundHead) tn.tensorHeadVqVertex[tensor1Index] = NULL;
    tn.tensorTailVqVertex[tensor1Index] = tail;
    if (tail != NULL) tail->nextInTensor = NULL;
   
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


/** Updates tensorIndexFromGlobalPq to be consistent with changes made to numPq on a tensor
 * @param[in,out] tn The tensor network object to update
 */ 
void remapTensorIndexFromGlobalPq(TensorNetwork tn){
    int currentPq = 0;
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        for (int j=0; j<tensor.numPq; j++){
            tn.tensorIndexFromGlobalPq[currentPq++] = i;
        }
    }

}

/** Updates firstGlobalPqIndex to be consistent with changes made to numPq on a tensor
 * @param[in,out] tn The tensor network object to update
 */ 
void remapFirstGlobalPqIndex(TensorNetwork tn){
    int globalPqOffset = 0;
    for (int i=0; i<tn.numTensors; i++){
        Tensor tensor = tn.tensors[i];
        tensor.firstGlobalPqIndex = globalPqOffset;
        globalPqOffset += tensor.numPq;
    }
}

/** Remove vqVertex objects from the linked list of virtual qubits in a tensor after those virtual
 * qubit indices have been contracted away. Free any vqVertex objects that have been removed in this way.
 * The first vertex to keep across both tensors should be marked as the head of the list. If this
 * function is called on the second tensor, the tail of the previous linked list should be linked to the
 * head of the second tensor's linked list.
 *
 * @param[in,out] TensorNetwork tn the tensor network object to update
 * @param[in] tensorIndex the tensor index of the linked list being traversed and freed. 
 * @param[in] newTensorIndex the index of the output tensor, ie the first tensor in the contraction
 * @param[in] uncontractedVqs list of virtual qubit indices which have not been contracted out
 * @param[in] numUncontractedVqs the number of elements in the uncontractedVqs list
 * @param[in] tail the location to store a pointer to the remaining vqVertex at the end of the list
 * @paramp[in,out] whether or not an uncontracted vqVertex has been found yet (eg in the previous tensor)
 */
void removeContractedVqVertices(TensorNetwork tn, int tensorIndex, int newTensorIndex,
        VqVertex *prevVqVertexToKeep,
        int *uncontractedVqs, int numUncontractedVqs, VqVertex **tail, int *foundHead){

    VqVertex *vqVertex, *prevVqVertex;
    int vqVertexCount = 0;
    int uncontractedVqCount = 0;
    vqVertex = tn.tensorHeadVqVertex[tensorIndex];
    // Search through vertices and keep those which are in uncontractedVqs.
    // Free vertices which no longer need to be stored
    while (vqVertex != NULL){
        if (uncontractedVqCount<numUncontractedVqs) {
            if (vqVertexCount++ == uncontractedVqs[uncontractedVqCount]) {
                // keep 
                if (vqVertex->tensorIndex != newTensorIndex) vqVertex->tensorIndex = newTensorIndex;
                if (!(*foundHead)) {
                    tn.tensorHeadVqVertex[tensorIndex] = vqVertex;
                    *foundHead = 1;
                }
                if (prevVqVertexToKeep!=NULL){
                    prevVqVertexToKeep->nextInTensor = vqVertex;
                }
                prevVqVertexToKeep = vqVertex;
                vqVertex = vqVertex->nextInTensor;
                uncontractedVqCount++;
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
    if (prevVqVertexToKeep != NULL){
        prevVqVertexToKeep->nextInTensor=NULL;
    }
} 

// ----- operations ------------------------------------------------------------


void tn_unitary(TensorNetwork tn, const int targetQubit, ComplexMatrix2 u){
    QCoord targetPqLocal = getLocalPq(tn, targetQubit);
    Tensor *targetTensor = &(tn.tensors[targetPqLocal.tensorIndex]);
    unitary(targetTensor->qureg, targetPqLocal.qIndex, u);
}


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
        DEBUG_PRINT(("controlled not tensor 1: %d %d\n", controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq));
        controlledNot(controlTensor->qureg, controlPqLocal.qIndex, virtualTargetIndex + controlTensor->numPq);
        
        // do target tensor half
        // TODO: note this won't work if doing multiple different operations on the same tensor in parallel 
        Tensor *targetTensor = &(tn.tensors[targetPqLocal.tensorIndex]);
        int virtualControlIndex = targetTensor->nextVqIndex;
        initVirtualControl(*targetTensor, virtualControlIndex + targetTensor->numPq);
        targetTensor->nextVqIndex = targetTensor->nextVqIndex + 1;
        DEBUG_PRINT(("controlled not tensor 2: %d %d\n", virtualControlIndex + targetTensor->numPq, targetPqLocal.qIndex));
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


/** Place target virtual qubit in the zero state
 * @param[in,out] tensor the tensor object
 * @param[in] virtual qubit to initialize. Index is local to a tensor but includes all physical qubits in the tensor
 */
void initVirtualTarget(Tensor tensor, int vqIndex){
     collapseToOutcome(tensor.qureg, vqIndex, 0);
}

/** Initialize virtual qubit to |0> + |1> 
 * NOTE: not (1/sqrt(2))(|0> + |1>)
 * @param[in,out] tensor the tensor object
 * @param[in] virtual qubit to initialize. Index is local to a tensor but includes all physical qubits in the tensor
 */
void initVirtualControl(Tensor tensor, int vqIndex){
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

/** Get the index of a physical qubit within a tensor from a global index which expects all
 * tensors to be ordered as in TensorNetwork.tensors 
 * 
 * @param[in] tn the tensor network object
 * @param[in] globalPq the index of the physical qubit in the global system
 * @returns the tensor and local index within the tensor of the qubit 
 */
QCoord getLocalPq(TensorNetwork tn, int globalPq){
    QCoord localPq;
    localPq.tensorIndex = tn.tensorIndexFromGlobalPq[globalPq];
    localPq.qIndex = globalPq - tn.tensors[localPq.tensorIndex].firstGlobalPqIndex;
    return localPq;
}

// ----- reporting ------------------------------------------------------------

/** Get the index of the virtual qubit in the list of virtual qubits in that tensor from a 
 * VqVertex object. Note that this is inefficient as it counts from the vqVertex to the end
 * of the linked list.
 *
 * @param[in] tn the tensor network object
 * @param[in] a pointer to the virtual qubit vertex to locate
 * @returns the index of the virtual qubit
 */
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

















