message(STATUS "Looking for CUDAToolkit")

set(LED_CUDAToolkit_MINIMUM_REQUIRED 12.0)
find_package(CUDAToolkit ${LED_GSL_MINIMUM_REQUIRED} REQUIRED)

message(STATUS "Looking for CUDAToolkit - found (version: ${CUDAToolkit_VERSION})")
