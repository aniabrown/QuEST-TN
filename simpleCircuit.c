#include <stdio.h>
#include "QuEST.h"

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

    printf("\nCircuit output:\n");

    qreal prob;
    prob = getProbAmp(qubits, 7);
    printf("Probability amplitude of |111>: %f\n", prob);

    prob = calcProbOfOutcome(qubits, 2, 1);
    printf("Probability of qubit 2 being in state 1: %f\n", prob);

    int outcome = measure(qubits, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);

    outcome = measureWithStats(qubits, 2, &prob);
    printf("Qubit 2 collapsed to %d with probability %f\n", outcome, prob);


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
