from QuESTPy.QuESTTypes import *
from .TNTypes import *

# Public API

createTensorNetwork  = TNTestee ("createTensorNetwork", retType=TensorNetwork, argType=[c_int, POINTER(c_int), POINTER(c_int), QuESTEnv])

contractTensorNetwork = TNTestee ("contractTensorNetwork", retType=None, argType=[TensorNetwork, QuESTEnv])

contractTensors = TNTestee ("contractTensors", retType=None, argType=[TensorNetwork, c_int, c_int, QuESTEnv])

tn_controlledNot = TNTestee ("tn_controlledNot", retType=None, argType=[TensorNetwork, c_int, c_int])

tn_unitary = TNTestee ("tn_unitary", retType=None, argType=[TensorNetwork, c_int, ComplexMatrix2])

printTensorNetwork = TNTestee ("printTensorNetwork", retType=None, argType=[TensorNetwork])


# Internal functions

getLocalPq = TNTestee ("getLocalPq", retType=QCoord, argType=[TensorNetwork, c_int])
