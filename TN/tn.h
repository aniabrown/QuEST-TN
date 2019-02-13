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


typedef struct VqCoord {
    int tensorIndex;
    int vqIndex;
} VqCoord;

typedef struct VqAdjacencyList {
    int numAdjacencies;
    VqCoord *adjacencyList;
}

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
    int *tensorIdFromGlobalPq;
    VqAdjacencyList *adjacencyList;
} TensorNetwork;


// ----- TENSOR NETWORK -------------------------------------------------------

void contractTensorNetwork(TensorNetwork tn);

/**
    replaces the tensor object at index tensor1 with the contraction of tensor1 and tensor 2
*/ 
void contractTensors(TensorNetwork tn, int tensor1, int tensor2);

TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor);

/** probably ineffient, but included for utility
 */
void addTensorToNetwork(int numPq, int numVq);


// ----- OPERATIONS -----------------------------------------------------------

void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit);

void tn_unitary(TensorNetwork tn, const int targetQubit, ComplexMatrix2 u);





#endif // TN_H

