find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(NOvA_SM_EventRates NOvA_SM_Eventrates.c++ NOvA_Setup.h)
add_executable(NOvA_LED_EventRates NOvA_LED_Eventrates.c++ NOvA_Setup.h)
add_executable(NOvA_SM_Fit NOvA_SM_Fit.c++ NOvA_Setup.h)
add_executable(NOvA_LED_Fit NOvA_LED_Fit.c++ NOvA_Setup.h)
add_executable(NOvA_LED_Scan NOvA_LED_Scan.c++ NOvA_Setup.h)

target_link_libraries(NOvA_SM_EventRates GSL::gsl globeslib LED_Lib MPI::MPI_C)
target_link_libraries(NOvA_LED_EventRates GSL::gsl globeslib LED_Lib MPI::MPI_C)
target_link_libraries(NOvA_SM_Fit GSL::gsl globeslib LED_Lib MPI::MPI_C)
target_link_libraries(NOvA_LED_Fit GSL::gsl globeslib LED_Lib MPI::MPI_C)
target_link_libraries(NOvA_LED_Scan GSL::gsl globeslib LED_Lib MPI::MPI_C)

if(OpenBLAS_FOUND)
     target_link_libraries(NOvA_SM_EventRates openblas)
     target_link_libraries(NOvA_LED_EventRates openblas)
endif()