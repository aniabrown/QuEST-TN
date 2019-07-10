#!/usr/bin/env python3
from QuESTPy.QuESTFunc import *
from TNPy.TNFunc import *

def TN_singleQubitGate(gate, tn, *args):
    pq = getLocalPq(tn, args[0]) 
    gate(tn.tensors[pq.tensorIndex].qureg, pq.qIndex, *args[1:]);

def TN_controlledGate(gate, tn, *args):
    controlQubit = args[0]
    targetQubit = args[1]
    is_local = getControlGateIsLocal(tn, controlQubit, targetQubit)
    controlPq = getLocalPq(tn, controlQubit) 
    targetPq = getLocalPq(tn, targetQubit) 
    if (is_local):
        gate(tn.tensors[controlPq.tensorIndex].qureg, controlPq.qIndex, targetPq.qIndex, *args[2:])
    else:
        controlTensor = getTensor(tn, controlQubit)
        virtualTargetIndex = incrementVqIndex(tn, controlQubit)
        initVirtualTarget(controlTensor, virtualTargetIndex)
        gate(controlTensor.qureg, controlPq.qIndex, virtualTargetIndex + controlTensor.numPq)

        targetTensor = getTensor(tn, targetQubit)
        virtualControlIndex = incrementVqIndex(tn, targetQubit)
        initVirtualControl(targetTensor, virtualControlIndex)
        gate(targetTensor.qureg, targetPq.qIndex, virtualControlIndex + targetTensor.numPq)

        updateTNForControlGate(tn, controlQubit, targetQubit)


