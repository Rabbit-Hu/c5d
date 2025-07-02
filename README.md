# C⁵D: Sequential Continuous Convex Collision Detection Using Cone Casting

This repository contains an implementation of the C⁵D algorithm proposed in the SIGGRAPH 2025 paper *C⁵D: Sequential Continuous Convex Collision Detection Using Cone Casting*. 

## Build (Linux)

```sh
mkdir build
cd build
cmake ..
make
```

CMake should download third-party dependencies (Eigen, CGAL, etc.) automatically. However, if you want to use your local installations, please modify the `CMakeLists.txt` file accordingly. 

## Build (Windows)

For Windows, it is recommended to use the [vcpkg](https://github.com/microsoft/vcpkg) package manager to install Boost and CGAL. For other dependencies, you can also use vcpkg or install them manually shall there be any issues with the automatic download.

```sh
vcpkg install boost
vcpkg install cgal
```

You will need to modify the `CMakeLists.txt` file to include the vcpkg toolchain file. Put the following line before the `project()` command in `CMakeLists.txt`:

```cmake
include("C:/Users/xiaodiyuan/Repos/vcpkg/scripts/buildsystems/vcpkg.cmake")
```

Then, replace the CGAL and Boost include and link commands with the following:

```cmake
# include(boost)
# target_link_libraries(c5d PUBLIC Boost::timer Boost::filesystem Boost::regex)
find_package(Boost REQUIRED COMPONENTS timer filesystem regex)
target_link_libraries(c5d PUBLIC Boost::timer Boost::filesystem Boost::regex)

# ...

# include(cgal)
# target_link_libraries(c5d PUBLIC CGAL::CGAL CGAL::CGAL_Core)
find_package(CGAL CONFIG REQUIRED COMPONENTS Core)
target_link_libraries(c5d PRIVATE CGAL::CGAL CGAL::CGAL_Core)
```

I am still working on the Windows build, so please let me know if you encounter any issues.

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
