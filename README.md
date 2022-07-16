# rdbench
2D reaction-diffusion system benchmark using MPI and MPI-IO.

![](https://raw.githubusercontent.com/range3/rdbench/master/rdbench-viz/viz.gif)

# Build and install
## CMake
```bash
git clone https://github.com/range3/rdbench.git
cd rdbench
mkdir build
cd build
cmake ..
cmake --build .
cmake --install . --prefix /path/to/install
```

## Spack
```bash
git clone https://github.com/range3/rdbench.git
spack repo add rdbench/spack-repo
spack install rdbench
```

# Usage
```bash
mkdir -p /path/to/output_dir
mpirun -np 4 rdbench -o /path/to/output_dir/output_file_prefix_
```
