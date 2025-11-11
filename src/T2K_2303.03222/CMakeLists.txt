add_executable(T2K_SM_EventRates T2K_SM_eventrates.c myio.c T2K_setup.h)

target_link_libraries(T2K_SM_EventRates GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_SM_EventRates openblas)
endif()

add_executable(T2K_LED_EventRates T2K_LED_eventrates.c myio.c nu5d_mat_osc.h T2K_setup.h)

target_link_libraries(T2K_LED_EventRates GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_LED_EventRates openblas)
endif()

add_executable(T2K_SM_Fit T2K_SM_Fit.c myio.c T2K_setup.h)

target_link_libraries(T2K_SM_Fit GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(T2K_SM_Fit openblas)
endif()