#include <stdio.h>
#include "QuEST.h"
#include "QuEST_debug.h"
#include "tn.h"

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
    pauliX(qubits, 0);
    rotateZ(qubits, 0, 0.3);
    //hadamard(qubits, 0);
    //hadamard(qubits, 2);
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
    controlledNot(qubits, 0, 1);
    controlledNot(qubits, 0, 2);
    controlledNot(qubits, 0, 1);
    reportStateToScreen(qubits, env, 0);


    /*
     * TENSOR NETWORK
     */
    int numPqPerTensor[3] = {1, 1, 1};
    int numVqPerTensor[3] = {4, 2, 2};
    printf("\n\n--- Create Tensor Network:\n");
    TensorNetwork tn = createTensorNetwork(3, numPqPerTensor, numVqPerTensor, env);
    printTensorNetwork(tn);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    printf("\n\n--- Apply controlled not:\n");
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 0, 1);
    tn_controlledNot(tn, 0, 2);
    tn_controlledNot(tn, 0, 1);
    printTensorNetwork(tn);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    reportStateToScreen(tn.tensors[2].qureg, env, 0);

    printf("\n\n--- Contract tensors 0 and 2:\n");
    contractTensors(tn, 0, 2, env);
    printTensorNetwork(tn);

    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);
    
    printf("\n\n--- Contract tensors 0 and 1:\n");
    contractTensors(tn, 0, 1, env);
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
