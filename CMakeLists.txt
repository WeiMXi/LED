cmake_minimum_required(VERSION 3.12)
project(LIV)

set(CMAKE_CXX_STANDARD 20)

option(LIV_BUILTIN_GLOBES "Use built-in GLoBES (network or pre-downloaded source is required)." OFF)

set(LIV_BUILTIN_GLOBES_VERSION 3.2.18 CACHE STRING "Set built-in Mustard version (if required).")

set(LIV_PROJECT_CMAKE_DIR ${PROJECT_SOURCE_DIR}/cmake)
set(LIV_PROJECT_3RDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)

# GSL
include(${LIV_PROJECT_CMAKE_DIR}/FindGSL.cmake)

# OpenBLAS
include(${LIV_PROJECT_CMAKE_DIR}/FindOpenBLAS.cmake)

# GLoBES
include(${LIV_PROJECT_CMAKE_DIR}/FindGLobes.cmake)

add_subdirectory(experiments)
add_subdirectory(data)

add_subdirectory(src)