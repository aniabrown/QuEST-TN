# ifndef QUEST_TN_INTERNAL
# define QUEST_TN_INTERNAL

# include "QuEST_tn.h"

// ----- TENSOR CONTRACTIONS ------------------------------------------------------

int getNumContractions(TensorNetwork tn, int tensor1Index, int tensor2Index);

void getContractionVqIndices(TensorNetwork tn, int tensor1Index, int tensor2Index, int numContractions,
    int **tensor1Contractions, int **tensor2Contractions,
    int *numTensor1UncontractedVqs, int *numTensor2UncontractedVqs,
    int **tensor1UncontractedVqs, int **tensor2UncontractedVqs);

long long int getContractedIndex(long long int activeStateVec1Index, long long int activeStateVec2Index,
     long long int stateVec1PhysicalSize, long long int stateVec2PhysicalSize, int numTensor1Pq, int numTensor2Pq,
     int numTensor1UncontractedVqs, int numTensor2UncontractedVqs);

long long int getStateVectorIndexFromActiveIndex(long long int activeStateIndex,
    int numPq, int *uncontractedVqs, int numUncontractedVqs);

Complex recursiveContract(Tensor tensor1, Tensor tensor2, long long int tensor1Offset,
    long long int tensor2Offset, int *tensor1Contractions, int *tensor2Contractions,
    int numContractions, int vqIndex);

// ----- GATES -----------------------------------------------------------

void initVirtualTarget(Tensor tensor, int vqIndex);

void initVirtualControl(Tensor tensor, int vqIndex);

// ----- VQ GRAPH UTILITY --------------------------------------------------------------

void removeAllVqVertices(TensorNetwork tn, int tensorIndex);

void removeContractedVqVertices(TensorNetwork tn, int tensorIndex, int newTensorIndex,
        VqVertex *startingVqVertex, VqVertex *prevVqVertexToKeep,
        int *uncontractedVqs, int numUncontractedVqs, VqVertex **tail, int *foundHead);


// ----- INDEXING UTILITY --------------------------------------------------------------

QCoord getLocalPq(TensorNetwork tn, int globalPq);

void remapTensorIndexFromGlobalPq(TensorNetwork tn);

void remapFirstGlobalPqIndex(TensorNetwork tn);

int getVqVertexIndex(TensorNetwork tn, VqVertex *vqVertex);



# endif // QUEST_TN_INTERNAL
