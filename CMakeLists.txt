cmake_minimum_required(VERSION 3.12)
project(LED)

set(CMAKE_CXX_STANDARD 20)

option(LED_BUILTIN_GLOBES "Use built-in GLoBES (network or pre-downloaded source is required)." OFF)

set(LED_BUILTIN_GLOBES_VERSION 3.2.18 CACHE STRING "Set built-in Mustard version (if required).")

set(LED_PROJECT_CMAKE_DIR ${PROJECT_SOURCE_DIR}/cmake)
set(LED_PROJECT_3RDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)

# GSL
include(${LED_PROJECT_CMAKE_DIR}/FindGSL.cmake)

# OpenBLAS
include(${LED_PROJECT_CMAKE_DIR}/FindOpenBLAS.cmake)

# GLoBES
include(${LED_PROJECT_CMAKE_DIR}/FindGLobes.cmake)

# CUDAToolkit
include(${LED_PROJECT_CMAKE_DIR}/FindCUDAToolkit.cmake)

add_subdirectory(experiments)
add_subdirectory(data)

add_subdirectory(src)