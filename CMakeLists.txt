cmake_minimum_required(VERSION 3.18.0)

project(
    C5D
    DESCRIPTION "Sequential Continuous Convex Collision Detection Using Cone Casting"
    LANGUAGES CXX
    VERSION 0.1.0
)

# Options
option(C5D_DOUBLE_PRECISION "Enable double precision for ABD" ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# place binaries and libraries according to GNU standards
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "DOUBLE PRECISION: ${C5D_DOUBLE_PRECISION}")

add_library(c5d)

# Use O0 optimization for Debug builds
set(C5D_OPTIMIZATION_LEVEL "O2" CACHE STRING "Optimization level for C5D")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Setting optimization level to O0 because build type is Debug")
    set(C5D_OPTIMIZATION_LEVEL "O0")
endif()
message(STATUS "Optimization level for C5D: ${C5D_OPTIMIZATION_LEVEL}")
target_compile_options(c5d PRIVATE "-${C5D_OPTIMIZATION_LEVEL}")

# Double precision
if(C5D_DOUBLE_PRECISION)
    message(STATUS "C5D is using double precision")
    target_compile_definitions(c5d PUBLIC C5D_DOUBLE_PRECISION)
else()
    message(STATUS "C5D is using single precision")
endif()

set(C5D_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
add_subdirectory("${C5D_SOURCE_DIR}")

# Add source to the library
set(C5D_SOURCES
    ${C5D_SOURCE_DIR}/c5d_exp.hpp
    ${C5D_SOURCE_DIR}/c5d_exp.cpp
    ${C5D_SOURCE_DIR}/ACCD.hpp
    ${C5D_SOURCE_DIR}/ACCD.cpp
    ${C5D_SOURCE_DIR}/gjk.hpp
    ${C5D_SOURCE_DIR}/gjk.cpp
)
source_group(TREE ${C5D_SOURCE_DIR} PREFIX "Source Files" FILES ${C5D_SOURCES})
target_sources(c5d PRIVATE ${C5D_SOURCES})

######## Libraries ########
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

include(boost)
target_link_libraries(c5d PUBLIC Boost::timer Boost::filesystem Boost::regex)

include(eigen)
target_link_libraries(c5d PUBLIC Eigen3::Eigen)

include(spdlog)
target_link_libraries(c5d PUBLIC spdlog::spdlog)

include(cgal)
target_link_libraries(c5d PUBLIC CGAL::CGAL CGAL::CGAL_Core)

include(libigl)
target_link_libraries(c5d PUBLIC igl::core igl::predicates)

include(json)
target_link_libraries(c5d PUBLIC nlohmann_json::nlohmann_json)

include(yaml-cpp)
target_link_libraries(c5d PUBLIC yaml-cpp::yaml-cpp)

# Set -Wall -Wextra
target_compile_options(c5d PUBLIC -Wall -Wextra)

