from QuESTPy.QuESTTypes import *
from .TNTypes import *

# Public API

contractIndices = TNTestee ("contractIndices", retType=Tensor, argType=[Tensor, Tensor, POINTER(c_int), POINTER(c_int), c_int, \
        POINTER(c_int), c_int, POINTER(c_int), c_int,
        QuESTEnv])

initVirtualTarget = TNTestee("initVirtualTarget", retType=None, argType=[Tensor, c_int])

initVirtualControl = TNTestee("initVirtualControl", retType=None, argType=[Tensor, c_int])

createTensor = TNTestee("createTensor", retType=Tensor, argType=[c_int, c_int, QuESTEnv])
