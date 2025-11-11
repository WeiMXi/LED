## How to run

Dependencies:
* C++ compiler (GCC, Clang, etl)
* GSL (≥ 2.6)
* CMake (≥ 3.12)
* OpenBlas (disable if not found)
* CUDAToolkit (disable if not found, import performance)

Required, built-in if not found (network or pre-downloaded source is required):
* [GLoBES, General Long Baseline Experiment Simulator](https://www.mpi-hd.mpg.de/personalhomes/globes/index.html), built-in version: [globes-cmake](https://github.com/WeiMXi/globes-cmake)  


for debian, you can run `sudo apt install libopenblas-dev cmake libgsl-dev` to install most of them
```shell
cd the-project-path
mkdir build
cd build
cmake ..
make
cd src/5d_bulk
./nu5d_mat_osc
```

after getting the output, you can run script in `build/src/data/prob/T2Kplot.py` to get the fig
