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
    initPlusState(qubits);
    rotateZ(qubits, 0, 0.3); 
    rotateZ(qubits, 1, 0.5); 
    rotateZ(qubits, 2, 0.7); 
    reportStateToScreen(qubits, env, 0);

    // Apply entangling operation
    printf("Apply controlledNot(0,2), controlledNot(0,1), controlledNot(0,2):\n");
    controlledNot(qubits, 0, 2);
    controlledNot(qubits, 0, 1);
    reportStateToScreen(qubits, env, 0);

    /*
     * TENSOR NETWORK
     */

    // Initialize
    int numPqPerTensor[3] = {1, 1, 1};
    int numVqPerTensor[3] = {2, 1, 1};
    printf("\n\nCreate Tensor Network representation:\n");
    TensorNetwork tn = createTensorNetwork(3, numPqPerTensor, numVqPerTensor, env);

    printf("Apply single qubit gates\n");
    // Apply single qubit gates -- we want a different probability amplitude at each
    // position in the state vector for testing purposes
    ComplexMatrix2 uRotateZ1, uRotateZ2, uRotateZ3, uHadamard;

    uHadamard.r0c0 = (Complex) {.real=1/sqrt(2), .imag=0};
    uHadamard.r0c1 = (Complex) {.real=1/sqrt(2), .imag=0};
    uHadamard.r1c0 = (Complex) {.real=1/sqrt(2), .imag=0};
    uHadamard.r1c1 = (Complex) {.real=-1/sqrt(2), .imag=0};

    qreal theta = 0.3;
    uRotateZ1.r0c0 = (Complex) {.real=cos(theta/2), .imag=-sin(theta/2)};
    uRotateZ1.r0c1 = (Complex) {.real=0, .imag=0};
    uRotateZ1.r1c0 = (Complex) {.real=0, .imag=0};
    uRotateZ1.r1c1 = (Complex) {.real=cos(theta/2), .imag=sin(theta/2)};

    theta = 0.5;
    uRotateZ2.r0c0 = (Complex) {.real=cos(theta/2), .imag=-sin(theta/2)};
    uRotateZ2.r0c1 = (Complex) {.real=0, .imag=0};
    uRotateZ2.r1c0 = (Complex) {.real=0, .imag=0};
    uRotateZ2.r1c1 = (Complex) {.real=cos(theta/2), .imag=sin(theta/2)};

    theta = 0.7;
    uRotateZ3.r0c0 = (Complex) {.real=cos(theta/2), .imag=-sin(theta/2)};
    uRotateZ3.r0c1 = (Complex) {.real=0, .imag=0};
    uRotateZ3.r1c0 = (Complex) {.real=0, .imag=0};
    uRotateZ3.r1c1 = (Complex) {.real=cos(theta/2), .imag=sin(theta/2)};

    tn_unitary(tn, 0, uHadamard);
    tn_unitary(tn, 1, uHadamard);
    tn_unitary(tn, 2, uHadamard);
    tn_unitary(tn, 0, uRotateZ1);
    tn_unitary(tn, 1, uRotateZ2);
    tn_unitary(tn, 2, uRotateZ3);

    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    printf("Tensor 2\n");
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    // Apply entangling operation
    printf("\n\nApply controlledNot(0,2), controlledNot(0,1), controlledNot(0,2):\n");
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 0, 1);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 1\n");
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    printf("Tensor 2\n");
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    // Contract

    printf("\n\nContract tensors 0 and 1:\n");
    contractTensors(tn, 0, 1, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    printf("Tensor 2\n");
    reportStateToScreen(tn.tensors[2].qureg, env, 0);


    printf("\n\nContract tensors 0 and 2:\n");
    contractTensors(tn, 0, 2, env);
    printTensorNetwork(tn);
    printf("Tensor 0\n");
    reportStateToScreen(tn.tensors[0].qureg, env, 0);

    /*
     * COMPARE PLAIN QuEST and TN RESULTS
     */

    int correctOutcome = compareStates(tn.tensors[0].qureg, qubits, 10e-10);
    printf("\n**************************************************\n");
    if (correctOutcome) printf("Passed test\n");
    else printf("Failed test\n");
    printf("**************************************************\n");


    /*
     * FREE MEMORY
     */

    destroyQureg(qubits, env); 


    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return !correctOutcome;
}
