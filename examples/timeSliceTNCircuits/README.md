QuEST-TN Tutorial
======

**Table of Contents**
- [Coding](#coding)
- [Compiling](#compiling)
- [Running](#running)


# Coding

To use the QuEST-TN extension in your C or C++ code, include

```C
#include <QuEST.h>
#include <QuEST_tn.h>
```

The currently available API for operations on tensor networks is available [here](https://aniabrown.github.io/QuEST-TN/QuEST__tn_8h.html). For operations on plain QuEST Qureg objects, see the API [here](https://quest-kit.github.io/QuEST/QuEST_8h.html) and the [QuEST tutorial](examples/README.md). 

Here's a very simple circuit which applies a controlled not on a two qubit state split over two
tensors and reports the resulting tensor network structure. 
 
```C
#include <QuEST.h>
#include <QuEST_tn.h>

int main(int narg, char *varg[]) {

  // load QuEST environment. Do this exactly once as the first step in your code.
  QuESTEnv env = createQuESTEnv();
  
  int numPqPerTensor[2] = {1, 1};
  int numVqPerTensor[2] = {1, 1};
  
  TensorNetwork tn = createTensorNetwork(2, numPqPerTensor, numVqPerTensor, env);
  printTensorNetwork(tn);
  tn_controlledNot(tn, 0, 1);
  printTensorNetwork(tn);
	
  // unload QuEST environment. Do this exactly one as the last step in your code. 
  destroyQuESTEnv(env); 
  return 0;
}
```
---------------------------

## Test circuits



----------------------------

# Compiling

QuEST uses CMake (3.1 or higher) as its build system.

To compile, make sure your circuit code is accessible from the root QuEST directory.
In the root directory, initially build using
```bash
mkdir build
cd build
cmake -DUSER_SOURCE="myCode1.c;myCode2.c" ..
make
```
Paths to target sources are set as a semi-colon separated list of paths to said sources relative to the root QuEST directory.

If you wish your executable to be named something other than `demo`, you can set this too by using:
```bash
cmake -DOUTPUT_EXE="myExecutable" ..
```

When using the cmake command as above, the -D[VAR=VALUE] option can be passed other options to further configure your build.

To compile your code to run on multi-CPU systems use
```bash
cmake -DDISTRIBUTED=1 ..
```

To compile for GPU, use
```bash
cmake -DGPUACCELERATED=1 -DGPU_COMPUTE_CAPABILITY=[COMPUTE_CAPABILITY] ..
```

Where COMPUTE_CAPABILITY is the compute cabability of your GPU, written without a decimal point. This can can be looked up at the [NVIDIA website](https://developer.nvidia.com/cuda-gpus). The default value is 30.

By default, QuEST will compile with OpenMP parallelism enabled if an OpenMP compatible compiler and OpenMP library can be found on your system (e.g. [GCC 4.9](https://gcc.gnu.org/gcc-4.9/changes.html)). Using distribution requires an MPI implementation is installed on your system, and GPU acceleration requires a CUDA compiler (`nvcc`). We've made a comprehensive list of compatible compilers which you can view [here](../tests/compilers/compatibility.md). This does not change your `COMPILER` setting - the makefile will choose the appropriate MPI and CUDA wrappers automatically.

You can additionally customise the precision with which the state-vector is stored.
```bash
cmake -DPRECISION=2 ..
```
Using greater precision means more precise computation but at the expense of additional memory requirements and runtime.
Checking results are unchanged when altaring the precision can be a great test that your calculations are sufficiently precise.


Please note that cmake caches these changes (per directory) so for any subsequent builds you should just type `make` from the build directory and the previously defined settings will be applied. If any parameters require changing, these can be redefined by:
```
cmake -D[VAR=VALUE] ..
```
as one would do in the initial configuration.

For a full list of available configuration parameters, use
```bash
cmake -LH ..
```

For manual configuration (not recommended) you can change the `CMakeLists.txt` in the root QuEST directory.

----------------------------

# Running

## locally

You can then call your code. From the build directory:
```bash
./myExecutable
```
If multithreading functionality was found when compiling, you can control how many threads your code uses by setting `OMP_NUM_THREADS`, ideally to the number of available cores on your machine
```bash
export OMP_NUM_THREADS=8
./myExecutable
```
QuEST will automatically allocate work between the given number of threads to speedup your simulation.
