# ifndef QUEST_TN
# define QUEST_TN

#include "QuEST.h"

/** @file
 * The public interface to QuEST-TN

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
    tensor1FreeIndices and tensor2FreeIndices.

**/

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Tensor {
    // QuEST qureg object containing the physical and virtual qubits in this tensor, including
    // currently unused virtual qubits up to the maximum number of virtual qubits in the tensor
    Qureg qureg;
    // Number of physical qubits in the tensor
    int numPq;
    // Maximum number of virtual qubits that can be used in the tensor
    int numVq;
    // The index of the first physical qubit in the tensor in the global system
    int firstGlobalPqIndex;
    // The next unused virtual qubit index
    int nextVqIndex;
} Tensor;

Tensor createTensor(int numPq, int numVq, QuESTEnv env);

// ----- TENSOR CONTRACTIONS -------------------------------------------------------


Tensor contractIndices(Tensor tensor1, Tensor tensor2,
        int *tensor1Contractions, int *tensor2Contractions, int numContractions,
        int *tensor1FreeIndices, int numTensor1FreeIndices,
        int *tensor2FreeIndices, int numTensor2FreeIndices,
        QuESTEnv env);


// ----- GATES -----------------------------------------------------------

void initVirtualTarget(Tensor tensor, int virtualTargetIndex);

void initVirtualControl(Tensor tensor, int virtualControlIndex);


#ifdef __cplusplus
}
#endif

#endif // QUEST_TN

