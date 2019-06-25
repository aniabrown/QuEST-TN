from ctypes import *
import os.path
import importlib

        

TNLib = None
def init_TNLib(TNPath = ""):
    global TNLib
    TNLib = CDLL(TNPath + "libTN.so")

