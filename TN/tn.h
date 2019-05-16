# ifndef TN_H
# define TN_H

#include "QuEST.h"

/**
    Naming conventions:

    vq = virtual qubit
    pq = physical qubit
    These are not capitalised by default to distinguish ordinary vars from structs such as
    VqVertex, which do have the first letter capitalised. 

    Unless marked as global in the name, all var names that refer to qubit indices refer 
    to local index within a tensor and qubit type. i.e. physical qubits in tensor 2 start from index 0, 
    and virtual qubits in tensor 2 also start from index 0.

    In a contraction, the two tensors are referred to as tensor1 and tensor2, where tensor1
    is the tensor with the lower index in TensorNetwork.tensors

    In a contraction, all shared (virtual qubit) indices between tensor1 and tensor2 will be contracted.
    The lists of virtual qubit indices corresponding to these contractions are referred to as 
    tensor1Contractions and tensor2Contractions. The lists of any other virtual qubits which are 
    not being contracted as they represent entanglements with other tensors, are referred to as
    tensor1UncontractedVqs and tensor2UncontractedVqs.

**/


typedef struct QCoord {
    int tensorIndex;
    int qIndex;
} QCoord;

typedef struct VqVertex{
    struct VqVertex *nextInTensor;
    struct VqVertex *entangledPair;
    int tensorIndex;
} VqVertex;

typedef struct Tensor {
    Qureg qureg;
    int numPq;
    int numVq;
    int firstGlobalPqIndex;
    int nextVqIndex;
} Tensor;

typedef struct TensorNetwork {
    int numTensors;
    Tensor *tensors;
    int *tensorIndexFromGlobalPq;
    VqVertex **tensorHeadVqVertex;
    VqVertex **tensorTailVqVertex;
    int *numEntanglements;
} TensorNetwork;



// ----- TENSOR NETWORK INITIALISATION -------------------------------------------------------

TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor,
        QuESTEnv env);

/** probably ineffient, but consider including for utility
 */
//void addTensorToNetwork(int numPq, int numVq);

// ----- TENSOR NETWORK CLEANUP -------------------------------------------------------


// ----- TENSOR CONTRACTIONS -------------------------------------------------------

void contractTensorNetwork(TensorNetwork tn, QuESTEnv env);

/**
    replaces the tensor object at index tensor1 with the contraction of tensor1 and tensor 2
*/ 
void contractTensors(TensorNetwork tn, int tensor1, int tensor2, QuESTEnv env);

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

void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit);

void tn_unitary(TensorNetwork tn, const int targetQubit, ComplexMatrix2 u);

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

// ----- reporting ------------------------------------------------------------


void printTensorNetwork(TensorNetwork tn);



#endif // TN_H

