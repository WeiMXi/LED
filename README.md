## How to run

Dependencies:
* C++ compiler (GCC, Clang, etl)
* GSL
* OpenBlas

Required, built-in if not found (network or pre-downloaded source is required):
* [GLoBES, General Long Baseline Experiment Simulator](https://www.mpi-hd.mpg.de/personalhomes/globes/index.html), built-in version [globes-cmake](https://github.com/WeiMXi/globes-cmake)

```shell
cd the-project-path
mkdir build
cd build
cmake ..
make
cd src/5d_bulk
./nu5d_mat_osc
```

then you get the output, run script in `build/src/data/prob/T2Kplot.py` to get a fig