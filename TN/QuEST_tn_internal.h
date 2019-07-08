# ifndef QUEST_TN_INTERNAL
# define QUEST_TN_INTERNAL

# include "QuEST_tn.h"

/** @file
 * Internal QuEST-TN functions

    Naming conventions:

    See QuEST_tn.h for public API naming conventions

    Free indices are defined for a particular contraction between two tensors and 
    refer to physical qubits and virtual qubits which will not be contracted away
    in that particular contraction. freeIndexEl refers to an element in the state vector that
    could be constructed for a system containing only those qubits. This is slightly awkward and
    arbitrary but is a way to distinguish between two types of indices -- contration index, ie
    index in the list of qubits and index in a state vector (called an el).

**/


// ----- TENSOR CONTRACTIONS ------------------------------------------------------

int getNumContractions(TensorNetwork tn, int tensor1Index, int tensor2Index);

void getContractionVqIndices(TensorNetwork tn, int tensor1Index, int tensor2Index, int numContractions,
    int **tensor1Contractions, int **tensor2Contractions,
    int *numTensor1UncontractedVqs, int *numTensor2UncontractedVqs,
    int **tensor1UncontractedVqs, int **tensor2UncontractedVqs);

long long int getContractedIndex(long long int tensor1FreeIndexSize, long long int tensor2FreeIndexSize,
     long long int tensor1PqSize, long long int tensor2PqSize, int numTensor1Pq, int numTensor2Pq,
     int numTensor1UncontractedVqs, int numTensor2UncontractedVqs);

long long int getStateVectorIndexFromFreeIndexEl(long long int freeIndexEl,
    int numPq, int *uncontractedVqs, int numUncontractedVqs);

Complex recursiveContract(Tensor tensor1, Tensor tensor2, long long int tensor1Offset,
    long long int tensor2Offset, int *tensor1Contractions, int *tensor2Contractions,
    int numContractions, int vqIndex);

// ----- GATES -----------------------------------------------------------

void initVirtualTarget(Tensor tensor, int virtualTargetIndex);

void initVirtualControl(Tensor tensor, int virtualControlIndex);

void updateTNForControlGate(TensorNetwork tn, const int controlQubit, const int targetQubit);

int incrementVqIndex(TensorNetwork tn, int globalPq);

// ----- VQ GRAPH UTILITY --------------------------------------------------------------

void removeContractedVqVertices(TensorNetwork tn, int tensorIndex, int newTensorIndex,
        VqVertex *prevVqVertexToKeep,
        int *uncontractedVqs, int numUncontractedVqs, VqVertex **tail, int *foundHead);


// ----- INDEXING UTILITY --------------------------------------------------------------

QCoord getLocalPq(TensorNetwork tn, int globalPq);

Tensor getTensor(TensorNetwork tn, int globalPq);

void remapTensorIndexFromGlobalPq(TensorNetwork tn);

void remapFirstGlobalPqIndex(TensorNetwork tn);

int getVqVertexIndex(TensorNetwork tn, VqVertex *vqVertex);



# endif // QUEST_TN_INTERNAL
