## How to run

Dependencies:
* C++ compiler (GCC, Clang, etl)
* GSL (≥ 2.6)
* CMake (≥ 3.12)
* OpenBlas (disable if not found)
* CUDAToolkit (disable in default, improves performance for get the probability and event rates)

Required, built-in if not found (network or pre-downloaded source is required):
* [GLoBES, General Long Baseline Experiment Simulator](https://www.mpi-hd.mpg.de/personalhomes/globes/index.html), built-in version: [globes-cmake](https://github.com/WeiMXi/globes-cmake)  


for debian, you can run `sudo apt install libopenblas-dev cmake libgsl-dev` to install most of them
```shell
cd the-project-path
mkdir build
cd build
cmake ..
make
```

the main programs will be built at ./build/src