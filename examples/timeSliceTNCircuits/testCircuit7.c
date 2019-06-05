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

    printf("-------------------------------------------------------\n");
    printf("Running QuEST tutorial:\n\t Basic circuit involving a system of 3 qubits.\n");
    printf("-------------------------------------------------------\n");

    /*
     * PREPARE QUBIT SYSTEM
     */

    Qureg qubits = createQureg(3, env);
    //initStateDebug(qubits);
    initZeroState(qubits);
    hadamard(qubits, 0);
    hadamard(qubits, 1);
    hadamard(qubits, 2);
    rotateZ(qubits, 0, 0.3);
    rotateZ(qubits, 1, 0.5);
    rotateZ(qubits, 2, 0.7);
    reportStateToScreen(qubits, env, 0);

    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\nThis is our environment:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    /*
     * APPLY CIRCUIT
     */

    controlledNot(qubits, 0, 2);
    controlledNot(qubits, 1, 2);
    controlledNot(qubits, 0, 1);
    controlledNot(qubits, 0, 2);
    controlledNot(qubits, 0, 1);
    reportStateToScreen(qubits, env, 0);


    /*
     * TENSOR NETWORK
     */
    int numPqPerTensor[3] = {1, 1, 1};
    int numVqPerTensor[3] = {4, 3, 3};
    printf("\n\n--- Create Tensor Network:\n");
    TensorNetwork tn = createTensorNetwork(3, numPqPerTensor, numVqPerTensor, env);
    printTensorNetwork(tn);

    // Apply single qubit gates -- we want a different probability amplitude at each
    // position in the state vector for testing purposes
    ComplexMatrix2 uPauliX, uRotateZ1, uRotateZ2, uRotateZ3, uHadamard;
    uPauliX.r0c0 = (Complex) {.real=0, .imag=0};
    uPauliX.r0c1 = (Complex) {.real=1, .imag=0};
    uPauliX.r1c0 = (Complex) {.real=1, .imag=0};
    uPauliX.r1c1 = (Complex) {.real=0, .imag=0};

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

    //tn_unitary(tn, 0, uPauliX);
    tn_unitary(tn, 0, uHadamard);
    tn_unitary(tn, 1, uHadamard);
    tn_unitary(tn, 2, uHadamard);
    tn_unitary(tn, 0, uRotateZ1);
    tn_unitary(tn, 1, uRotateZ2);
    tn_unitary(tn, 2, uRotateZ3);

    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    printf("\n\n--- Apply controlled not:\n");
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 1, 2);
    tn_controlledNot(tn, 0, 1);
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 0, 1);
    printTensorNetwork(tn);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    printf("\n\n--- Contract tensors 0 and 1:\n");
    contractTensors(tn, 0, 1, env);
    printTensorNetwork(tn);

    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[2].qureg, env, 0);
    
    printf("\n\n--- Contract tensors 0 and 2:\n");
    contractTensors(tn, 0, 2, env);
    printTensorNetwork(tn);

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
