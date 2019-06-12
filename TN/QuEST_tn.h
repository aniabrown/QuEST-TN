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
    tensor1UncontractedVqs and tensor2UncontractedVqs.

**/


typedef struct QCoord {
    // The tensor that this qubit is located in
    int tensorIndex;
    // The local index of this qubit
    int qIndex;
} QCoord;

typedef struct VqVertex{
    // The vertex corresponding to the next virtual qubit in the tensor
    struct VqVertex *nextInTensor;
    // The vertex corresponding to the virtual qubit paired with this one in another tensor
    struct VqVertex *entangledPair;
    // The index of the tensor object that this virtual qubit is in
    int tensorIndex;
} VqVertex;

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

typedef struct TensorNetwork {
    // The number of tensors in the tensor network
    int numTensors;
    // An array of all tensor objects in the network
    Tensor *tensors;
    // For each qubit in the global system, the tensor which contains the qubit referenced by that
    // global qubit ID. Structured in this way for quick lookup during gate operations
    int *tensorIndexFromGlobalPq;
    // For each tensor, the first element in the linked list of virtual qubits in use in that tensor
    VqVertex **tensorHeadVqVertex;
    // For each tensor, the last element in the linked list of virtual qubits in use in that tensor.
    // Stored to be able to append to the list quickly.
    VqVertex **tensorTailVqVertex;
    // For each tensor, the number of virtual qubits currently in use representing entanglements
    // between that tensor and others
    int *numEntanglements;
} TensorNetwork;

// ----- TENSOR NETWORK INITIALISATION -------------------------------------------------------

/** Allocate memory for and initialize a new tensor network object given the number
 * of physical qubits and maximum number of virtual qubits in each tensor. Initialise all
 * physical qubits to the zero state and virtual qubits to an undefined state.
 *
 * @param[in] numTensors the number of tensors in the tensor network
 * @param[in] numPqPerTensor an array of values corresponding to the number of physical
 *      qubits in each tensor
 * @param[in] numVqPerTensor an array of values corresponding to the maximum number of
 *      virtual qubits that could be used in each tensor
 * @param[in] env QuEST environment object
 *
 */
TensorNetwork createTensorNetwork(int numTensors, int *numPqPerTensor, int *numVqPerTensor,
        QuESTEnv env);


// ----- TENSOR NETWORK CLEANUP -------------------------------------------------------

// TODO


// ----- TENSOR CONTRACTIONS -------------------------------------------------------


/** Contract all tensors in a tensor network one by one.
 * The simplest strategy -- contract tensor 0 with every other tensor in turn
 * Output will be stored as tensor 0.
 *
 * @param[in,out] tn tensor network to contract
 * @param[in] env QuEST environment object
 */
void contractTensorNetwork(TensorNetwork tn, QuESTEnv env);

/**
    Replaces the tensor object at index tensor1 in the tensor network with the contraction 
    of tensor1 and tensor 2. 
    All virtual qubits linking the two tensors are contracted in this operation. There may
    be virtual qubits left over after the contraction, corresponding to entanglement with 
    other tensors in the network.

    The ordering of qubits in the output tensor from lowest to highest is as follows: 
    tensor1 physical qubits, tensor2 physical qubits, tensor1 virtual qubits, 
    tensor2 virtual qubits.

    @param[in, out] tn the tensor network the two tensors belong to
    @param[in] tensor1 the index of the first tensor to contract
    @param[in] tensor1 the index of the second tensor to contract
    @param[in] env the QuEST environment object 
*/ 
void contractTensors(TensorNetwork tn, int tensor1, int tensor2, QuESTEnv env);


// ----- GATES -----------------------------------------------------------


/** Apply the controlled not (single control, single target) gate, also
 * known as the c-X, c-sigma-X, c-Pauli-X and c-bit-flip gate on qubits in
 * a tensor network.
 *
 * If the two qubits are in different tensors, also track the entanglement
 * introduced between those tensors using virtual qubits. 
 *
 * Note that the control and target qubits should be supplied as global qubit indices
 * over all physical qubits in all tensors.
 *
 * This applies pauliX to the target qubit if the control qubit has value 1.
 * This effects the two-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & & 1 \\
 * & & 1
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};
                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, -.5);

                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.5);
                \end{tikzpicture}
    }
    \f]
 *
 * @param[in,out] TensorNetwork tn the tensor network containing the control and target qubits
 * @param[in] controlQubit nots the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented), or are equal.
 */
void tn_controlledNot(TensorNetwork tn, const int controlQubit, const int targetQubit);



/** Apply a general single-qubit unitary (including a global phase factor).
 * The passed 2x2 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
 * Note that the target qubit should be supplied as a global qubit index
 * over all physical qubits in all tensors.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    }
    \f]
 *
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or matrix \p u is not unitary.
 */
void tn_unitary(TensorNetwork tn, const int targetQubit, ComplexMatrix2 u);


// ----- reporting ------------------------------------------------------------


/** Print all tensors in a tensor network and their entanglement with other tensors
 *
 * @param[in] tn the tensor network to print
 */
void printTensorNetwork(TensorNetwork tn);



#endif // QUEST_TN

