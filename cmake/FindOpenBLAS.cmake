message(STATUS "Looking for OpenBLAS")

set(LIV_OPENBLAS_MINIMUM_REQUIRED 0.3.0)
find_package(OpenBLAS ${LIV_OPENBLAS_MINIMUM_REQUIRED} REQUIRED)


message(STATUS "Looking for OpenBLAS - found (version: ${OpenBLAS_VERSION})")
