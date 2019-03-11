#include <stdio.h>
#include "QuEST.h"
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

    Qureg qubits = createQureg(4, env);
    initZeroState(qubits);

    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\nThis is our environment:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    /*
     * APPLY CIRCUIT
     */

    controlledNot(qubits, 0, 1);
    reportStateToScreen(qubits, env, 0);


    /*
     * TENSOR NETWORK
     */
    int numPqPerTensor[2] = {2, 2};
    int numVqPerTensor[2] = {1, 1};
    printf("\n\nCreate Tensor Network:\n");
    TensorNetwork tn = createTensorNetwork(2, numPqPerTensor, numVqPerTensor, env);
    printTensorNetwork(tn);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);

    printf("\n\nApply controlled not:\n");
    tn_controlledNot(tn, 1, 2);
    printTensorNetwork(tn);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);

    printf("\n\nContract tensors:\n");
    contractTensors(tn, 0, 1, env);
    reportStateToScreen(tn.tensors[0].qureg, env, 0);
    reportStateToScreen(tn.tensors[1].qureg, env, 0);



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
