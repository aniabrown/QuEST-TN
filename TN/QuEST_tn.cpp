#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>
#include "QuEST_debug.h"
#include "QuEST_tn_internal.h"
#include "QuEST_tn.h"
#include <cblas.h>

#define DEBUG 0

#if DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif

int* getTensorIndexPermutation(int* contractionIndices, int numContractions,
        int* freeIndices, int numFreeIndices, int tensor){

    int *perm = (int*) malloc((numContractions+numFreeIndices)*sizeof(int));

    int firstSize, secondSize;
    int *firstArray, *secondArray;
    if (tensor==1){
        firstSize=numFreeIndices;
        firstArray=freeIndices;
        secondSize=numContractions;
        secondArray=contractionIndices;
    } else {
        secondSize=numFreeIndices;
        secondArray=freeIndices;
        firstSize=numContractions;
        firstArray=contractionIndices;
    }

    for (int i=0; i<firstSize; i++){
        perm[i]=firstArray[i];
    }
    for (int i=0; i<secondSize; i++){
        perm[i+firstSize]=secondArray[i];
    }
    return perm;
}

/*
 * All indices are of size 2
 */
int* getTensorSizes(int numQubits){
    int *sizes = (int*) malloc(numQubits*sizeof(int));
    for (int i=0; i<numQubits; i++) sizes[i]=2;
    return sizes;
}


Tensor createTensor(int numPq, int numVq, QuESTEnv env){
    Tensor tensor;

    //! Probably don't need to store this any more
    tensor.numPq = numPq;
    tensor.numVq = numVq;
    tensor.qureg = createQureg(tensor.numPq + tensor.numVq, env);
    //! might need to store next available vq

    return tensor;
}

/*
 * Swap corresponding qubit IDs in perm and return new index
 */
int swapBits(long long int b, int *perm, int nQubits) {
  long long int out = 0;
  for(int j=0;j<nQubits; j++) {
    out += (b >> j & 1) * (1 << perm[j]);
  }
  return out;
}

qreal* permuteArray(Qureg qureg, int *perm) {
  qreal* outArr = (qreal *) malloc(sizeof(qreal)*2*qureg.numAmpsPerChunk);
  for (long long int i = 0; i < qureg.numAmpsPerChunk; i++) {
    long long int newIndex = 2*swapBits(i, perm, qureg.numQubitsRepresented);
    outArr[newIndex] = qureg.stateVec.real[i];
    outArr[newIndex+1] = qureg.stateVec.imag[i];
  }
  return outArr;
}

/*
 * Expects tensor1Contractions to be in order from smallest to largest.
 * Expects tensor2Contractions to be ordered to match indices in tensor1.
 * Ie if tensor1 and tensor2 each have 3 indices and (tensor1, index 1) is
 * contracted with (tensor2, index 2) and (tensor1, index 2) is contracted with
 * (tensor2, index 0), use:
 * tensor1Contractions = [1, 2]
 * tensor2Contractions = [2, 0]
 */
Tensor contractIndices(Tensor tensor1, Tensor tensor2,
        int *tensor1Contractions, int *tensor2Contractions, int numContractions,
        int *tensor1FreeIndices, int numTensor1FreeIndices,
        int *tensor2FreeIndices, int numTensor2FreeIndices,
        QuESTEnv env){

    printf("Contract with BLAS\n");

    int numTensor1Qubits, numTensor2Qubits;
    numTensor1Qubits = tensor1.qureg.numQubitsRepresented;
    numTensor2Qubits = tensor2.qureg.numQubitsRepresented;
    int totalNumQ = numTensor1FreeIndices + numTensor2FreeIndices;

    int *tensor1Perm = getTensorIndexPermutation(tensor1Contractions, numContractions,
						 tensor1FreeIndices, numTensor1FreeIndices, 1);
    int *tensor2Perm = getTensorIndexPermutation(tensor2Contractions, numContractions,
						 tensor2FreeIndices, numTensor2FreeIndices, 2);

    qreal* tensor1StateVecPermuted = permuteArray(tensor1.qureg, tensor1Perm);
    qreal* tensor2StateVecPermuted = permuteArray(tensor2.qureg, tensor2Perm);
    printf("finished permuting\n");

    // reportStateToScreen(tensor2.qureg, env, 0);

    qreal* outputQureg = (qreal *) malloc(sizeof(qreal)*(1LL << (totalNumQ+1LL)));

    Qureg contractedQureg = createQureg(totalNumQ, env);
    int M, N, K;
    M = 1 << numTensor2FreeIndices;
    N = 1 << numTensor1FreeIndices;
    K = 1 << numContractions;
    qreal alpha[2] = {1.0, 0.0};
    qreal beta[2] = {0.0, 0.0};

    printf("dims: %d %d %d\n", M, N, K);

    printTensor(tensor1StateVecPermuted, 1LL<<(numTensor1Qubits+1));
    printTensor(tensor2StateVecPermuted, 1LL<<(numTensor2Qubits+1));
    // tensor2 x tensor1.
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, alpha,
    		tensor2StateVecPermuted, K, tensor1StateVecPermuted, N, beta, outputQureg, N);

    // Copy data back
    for(long long int index = 0; index < contractedQureg.numAmpsPerChunk; index++) {
      contractedQureg.stateVec.real[index] = outputQureg[2*index];
      contractedQureg.stateVec.imag[index] = outputQureg[2*index+1];
    }

    free(tensor1StateVecPermuted);
    free(tensor2StateVecPermuted);
    free(outputQureg);

    free(tensor1Perm);
    free(tensor2Perm);

    // Free old quregs
    destroyQureg(tensor1.qureg, env);
    destroyQureg(tensor2.qureg, env);

    printf("freed memory\n");

    // Create output tensor
    Tensor outputTensor;
    outputTensor.qureg = contractedQureg;
    //! Set total number of qubits here -- don't know number of virtual/physical qubits
    //outputTensor.numPq = totalNumPq;
    //outputTensor.numVq = totalNumVq;

    outputTensor.qureg = contractedQureg;
    return outputTensor;
}

// ----- operations ------------------------------------------------------------


int getControlGateIsLocal(TensorNetwork tn, const int controlQubit, const int targetQubit){
    QCoord controlPqLocal = getLocalPq(tn, controlQubit);
    QCoord targetPqLocal = getLocalPq(tn, targetQubit);
    if (controlPqLocal.tensorIndex == targetPqLocal.tensorIndex) return 1;
    else return 0;
}


/** Place target virtual qubit in the zero state
 * @param[in,out] tensor the tensor object
 */
void initVirtualTarget(Tensor tensor, int virtualTargetIndex){
         int vqIndex = virtualTargetIndex + tensor.numPq;
     printf("vqIndex: %d\n", vqIndex);

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
            stateVecReal[index+stateVecSize] = 0;
            stateVecImag[index+stateVecSize] = 0;
        }
    }
}

/** Initialize virtual qubit to |0> + |1>
 * NOTE: not (1/sqrt(2))(|0> + |1>)
 * @param[in,out] tensor the tensor object
 * @param[in] virtual qubit to initialize. Index is local to a tensor but includes all physical qubits in the tensor
 */
void initVirtualControl(Tensor tensor, int virtualControlIndex){
    int vqIndex = virtualControlIndex + tensor.numPq;
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

