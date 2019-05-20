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
    printf("Running QuEST test circuit 1:\n\t Basic circuit involving a system of 2 qubits in plain\n\t QuEST Qureg and tensor network representation.\n");
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
    Qureg qubits = createQureg(2, env);
    printf("\nCreate a plain 2 qubit Qureg object:\n");
    reportQuregParams(qubits);

    // Apply single qubit gates
    initZeroState(qubits);
    pauliX(qubits, 0);
    rotateZ(qubits, 0, 0.3); 
    reportStateToScreen(qubits, env, 0);

    // Apply entangling operation
    printf("\nAppy controlledNot(0,1)\n");
    controlledNot(qubits, 0, 1);
    reportStateToScreen(qubits, env, 0);


    /*
     * TENSOR NETWORK
     */

    // Initialize
    int numPqPerTensor[2] = {1, 1};
    int numVqPerTensor[2] = {1, 1};
    printf("\n\nCreate Tensor Network representation:\n");
    TensorNetwork tn = createTensorNetwork(2, numPqPerTensor, numVqPerTensor, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);

    // Apply entangling operation
    printf("\n\nApply controlledNot(0,1):\n");
    tn_controlledNot(tn, 0, 1);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);

/*
    // Apply single qubit unitary gate to qubit 1
    ComplexMatrix2 u;
    u.r0c0 = (Complex) {.real=.5, .imag=.5};
    u.r0c1 = (Complex) {.real=.5, .imag=-.5};
    u.r1c0 = (Complex) {.real=.5, .imag=-.5};
    u.r1c1 = (Complex) {.real=.5, .imag=.5};
    tn_unitary(tn, 1, u);
*/

    // Contract
    printf("\n\nContract tensors:\n");
    contractTensors(tn, 0, 1, env);
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
