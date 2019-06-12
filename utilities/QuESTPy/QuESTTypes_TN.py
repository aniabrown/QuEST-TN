from .QuESTTypes import *
from .QuESTBase_TN import *

class QCoord(Structure):
    _fields_ = [("tensorIndex", c_int), 
                ("qIndex", c_int)]

class VqVertex(Structure):
    pass
# Work around to allow a self-referencial struct
VqVertex._fields_ = [("nextInTensor", POINTER(VqVertex)),
                ("entangledPair", POINTER(VqVertex)),
                ("tensorIndex", c_int)]

class Tensor(Structure):
    _fields_ = [("qureg", Qureg), 
                ("numPq", c_int),
                ("numVq", c_int),
                ("firstGlobalPqIndex", c_int),
                ("nextVqIndex", c_int)]

class TensorNetwork(Structure):
    _fields_ = [("numTensors", c_int),
                ("tensors", POINTER(Tensor)),
                ("tensorIndexFromGlobalPq", POINTER(VqVertex)),
                ("tensorHeadVqVertex", POINTER(POINTER(VqVertex))),
                ("tensorTailVqVertex", POINTER(POINTER(VqVertex))),
                ("numEntanglements", POINTER(c_int))]


class TNTestee(QuESTTestee):
    def __init__(self, funcname=None, retType=None, argType=[], defArg=[], denMat=None):
        if not TNLib[funcname]:
            raise IOError(funcname+' not found in TN API')
        self.thisFunc = TNLib[funcname]
        self.thisFunc.restype = retType
        self.thisFunc.argtypes = argType

    def __call__(self,*argsList):
        # If packed as list, otherwise receive as variables
        if len(argsList) == 1 and isinstance(argsList[0],list):
            specArg = argsList[0]
        else:
            specArg = list(argsList)

        return self.thisFunc(*specArg)


