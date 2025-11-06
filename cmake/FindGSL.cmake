message(STATUS "Looking for GSL")

set(LIV_GSL_MINIMUM_REQUIRED 2.6)
find_package(GSL ${LIV_GSL_MINIMUM_REQUIRED} REQUIRED)

message(STATUS "Looking for GSL - found (version: ${GSL_VERSION})")
