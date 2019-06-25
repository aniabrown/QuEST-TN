from ctypes import *
from QuESTPy.QuESTBase import init_QuESTLib
from TNPy.TNBase import init_TNLib
from QuESTPy.QuESTLibDir import defaultQuESTDir

# If set up
QuESTPath = "build/TN/QuEST/"
TNPath = "build/TN/"

init_QuESTLib(QuESTPath)
init_TNLib(TNPath)

print("RUN INIT")


