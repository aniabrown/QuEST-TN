#include <stdio.h>
#include "QuEST.h"
#include "QuEST_debug.h"
#include "QuEST_tn.h"

int main (int narg, char *varg[]) {

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env = createQuESTEnv();

    printf("-----------------------------------------------------------------\n");
    printf("Running QuEST test circuit 2:\n\t Basic circuit involving a system of 3 qubits in plain\n\t QuEST Qureg and tensor network representation.\n");
    printf("-----------------------------------------------------------------\n");

    /*
     * REPORT ENVIRONMENT
     */

    printf("\nThis is our environment:\n");
    reportQuESTEnv(env);

    /*
     * PLAIN QUREG REPRESENTATION
     */

    // Initialise 
    Qureg qubits = createQureg(3, env);
    printf("\nCreate a plain 3 qubit Qureg object:\n");
    reportQuregParams(qubits);

    // Apply single qubit gates
    initZeroState(qubits);
    pauliX(qubits, 0);
    rotateZ(qubits, 0, 0.3); 
    reportStateToScreen(qubits, env, 0);

    // Apply entangling operation
    printf("\nAppy controlledNot(0,1)\n");
    controlledNot(qubits, 0, 2);
    controlledNot(qubits, 0, 1);
    controlledNot(qubits, 0, 2);
    reportStateToScreen(qubits, env, 0);


    /*
     * TENSOR NETWORK
     */

    // Initialize
    int numPqPerTensor[3] = {1, 1, 1};
    int numVqPerTensor[3] = {3, 1, 2};
    printf("\n\nCreate Tensor Network representation:\n");
    TensorNetwork tn = createTensorNetwork(3, numPqPerTensor, numVqPerTensor, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    printf("Tensor 2\n");
    reportStateToScreen(tn.tensors[2].qureg, env, 0);


    // Apply single qubit gates -- we want a different probability amplitude at each
    // position in the state vector for testing purposes
    ComplexMatrix2 uPauliX, uRotateZ;
    uPauliX.r0c0 = (Complex) {.real=0, .imag=0};
    uPauliX.r0c1 = (Complex) {.real=1, .imag=0};
    uPauliX.r1c0 = (Complex) {.real=1, .imag=0};
    uPauliX.r1c1 = (Complex) {.real=0, .imag=0};

    qreal theta = 0.3;
    uRotateZ.r0c0 = (Complex) {.real=cos(theta/2), .imag=-sin(theta/2)};
    uRotateZ.r0c1 = (Complex) {.real=0, .imag=0};
    uRotateZ.r1c0 = (Complex) {.real=0, .imag=0};
    uRotateZ.r1c1 = (Complex) {.real=cos(theta/2), .imag=sin(theta/2)};

    tn_unitary(tn, 0, uPauliX);
    tn_unitary(tn, 0, uRotateZ);

    // Apply entangling operation
    printf("\n\nApply controlledNot(0,2), controlledNot(0,1), controlledNot(0,2):\n");
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 0, 1);
    tn_controlledNot(tn, 0, 2);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    printf("Tensor 2\n");
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    // Contract


    printf("\n\nContract tensors 0 and 2:\n");
    contractTensors(tn, 0, 2, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);


    printf("\n\nContract tensors 0 and 1:\n");
    contractTensors(tn, 0, 1, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);

    /*
     * FREE MEMORY
     */

    destroyQureg(qubits, env); 


    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
