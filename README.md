# C⁵D: Sequential Continuous Convex Collision Detection Using Cone Casting

This repository contains an implementation of the C⁵D algorithm proposed in the SIGGRAPH 2025 paper *C⁵D: Sequential Continuous Convex Collision Detection Using Cone Casting*. 

## Build

```sh
mkdir build
cd build
cmake ..
make
```

CMake should download third-party dependencies (Eigen, CGAL, etc.) automatically. However, if you want to use your local installations, please modify the `CMakeLists.txt` file accordingly.

## Usage

The executable `datagen` generates synthetic data to test the C⁵D algorithm. 

```sh
cd ..
mkdir data
./build/bin/datagen configs/base.yaml data
```

The executable `compare` compares the variants of the C⁵D algorithm with ACCD, the baseline algorithm. The results will be saved in csv files. 

```sh
mkdir outputs
./build/bin/compare data/base outputs/base
```

You can then use the python script `scripts/plot_curve.py` to plot the results.

```sh
python scripts/plot_curve.py 
```
