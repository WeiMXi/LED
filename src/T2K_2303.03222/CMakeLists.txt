find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
add_executable(T2K_SM_EventRates T2K_SM_Eventrates.c myio.c T2K_Setup.h)

target_link_libraries(T2K_SM_EventRates GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_SM_EventRates openblas)
endif()

add_executable(T2K_LED_EventRates T2K_LED_Eventrates.c++ T2K_Setup.h)

target_link_libraries(T2K_LED_EventRates GSL::gsl globeslib LED_Lib MPI::MPI_C)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_LED_EventRates openblas)
endif()

add_executable(T2K_SM_Fit T2K_SM_Fit.c myio.c T2K_Setup.h)

target_link_libraries(T2K_SM_Fit GSL::gsl globeslib MPI::MPI_C)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_SM_Fit openblas)
endif()

add_executable(T2K_LED_Scan T2K_LED_Scan.c++ T2K_Setup.h)

target_link_libraries(T2K_LED_Scan GSL::gsl globeslib LED_Lib MPI::MPI_C)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_LED_Scan openblas)
endif()