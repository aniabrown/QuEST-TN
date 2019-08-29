from QuESTPy.QuESTTypes import *
from .TNTypes import *

# Public API

createTensorNetwork  = TNTestee ("createTensorNetwork", retType=TensorNetwork, argType=[c_int, POINTER(c_int), POINTER(c_int), QuESTEnv])

contractTensorNetwork = TNTestee ("contractTensorNetwork", retType=None, argType=[TensorNetwork, QuESTEnv])

contractTensors = TNTestee ("contractTensors", retType=None, argType=[TensorNetwork, c_int, c_int, QuESTEnv])

contractIndices = TNTestee ("contractIndices", retType=Tensor, argType=[Tensor, Tensor, POINTER(c_int), POINTER(c_int), c_int, \
        POINTER(c_int), c_int, POINTER(c_int), c_int),
        QuESTEnv])

tn_controlledNot = TNTestee ("tn_controlledNot", retType=None, argType=[TensorNetwork, c_int, c_int])

tn_unitary = TNTestee ("tn_unitary", retType=None, argType=[TensorNetwork, c_int, ComplexMatrix2])

printTensorNetwork = TNTestee ("printTensorNetwork", retType=None, argType=[TensorNetwork])


# Internal functions

getLocalPq = TNTestee ("getLocalPq", retType=QCoord, argType=[TensorNetwork, c_int])

getTensor = TNTestee("getTensor", retType=Tensor, argType=[TensorNetwork, c_int])

incrementVqIndex = TNTestee("incrementVqIndex", retType=c_int, argType = [TensorNetwork, c_int])

initVirtualTarget = TNTestee("initVirtualTarget", retType=None, argType=[Tensor, c_int])

initVirtualControl = TNTestee("initVirtualControl", retType=None, argType=[Tensor, c_int])

updateTNForControlGate = TNTestee("updateTNForControlGate", retType=None, argType=[TensorNetwork, c_int, c_int])

getControlGateIsLocal = TNTestee("getControlGateIsLocal", retType=c_int, argType=[TensorNetwork, c_int, c_int])
