# ifndef TN_H
# define TN_H

#include "QuEST.h"

/**
    Naming conventions:
    vq = virtual qubit
    pq = physical qubit

    Unless marked as global, all indices refer to local index within a tensor and qubit type. 
    i.e. physical qubits in tensor 2 start from index 0, and virtual qubits in
    tensor 2 also start from index 0.

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


// ----- TENSOR NETWORK -------------------------------------------------------

void contractTensorNetwork(TensorNetwork tn, QuESTEnv env);

/**
    replaces the tensor object at index tensor1 with the contraction of tensor1 and tensor 2
*/ 
void contractTensors(TensorNetwork tn, int tensor1, int tensor2, QuESTEnv env);

TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor,
        QuESTEnv env);

/** probably ineffient, but included for utility
 */
void addTensorToNetwork(int numPq, int numVq);

void remapTensorIndexFromGlobalPq(TensorNetwork tn);

void remapFirstGlobalPqIndex(TensorNetwork tn);

void removeAllVqVertices(TensorNetwork tn, int tensorIndex);

void removeContractedVqVertices(TensorNetwork tn, int tensorIndex, VqVertex *startingVqVertex,
        int *unusedContractions, VqVertex **tail, int *foundHead);


// ----- OPERATIONS -----------------------------------------------------------

void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit);

void tn_unitary(TensorNetwork tn, const int targetQubit, ComplexMatrix2 u);

void initVirtualTarget(Tensor tensor, int vqIndex);
void initVirtualControl(Tensor tensor, int vqIndex);

// ----- utility --------------------------------------------------------------

QCoord getLocalPq(TensorNetwork tn, int globalPq);

// ----- reporting ------------------------------------------------------------

void printTensorNetwork(TensorNetwork tn);



#endif // TN_H

