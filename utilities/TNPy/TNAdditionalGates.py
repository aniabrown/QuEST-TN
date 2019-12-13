#!/usr/bin/env python3
from QuESTPy.QuESTFunc import *
from TNPy.TNFunc import *

def TN_singleQubitGate(gate, tensor, *args):
    gate(tensor.qureg, *args)

def TN_controlledGateLocal(gate, tensor, *args):
    controlQubit = args[0]
    targetQubit = args[1]
    
    gate(tensor.qureg, controlQubit, targetQubit, *args[2:])

def TN_controlledGateTargetHalf(gate, tensor, *args):
    controlQubit = args[0]
    targetQubit = args[1]

    numPq=1

    initVirtualControl(tensor, controlQubit-numPq)
    gate(tensor.qureg, controlQubit, targetQubit)


def TN_controlledGateControlHalf(gate, tensor, *args):
    controlQubit = args[0]
    targetQubit = args[1]
    
    numPq=1
    
    initVirtualTarget(tensor, targetQubit-numPq)
    gate(tensor.qureg, controlQubit, targetQubit)


