# file(GLOB_RECURSE NOVA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

add_executable(NOvA_SM_EventRates NOvA_SM_eventrates.c myio.c NOvA_setup.h)
add_executable(NOvA_LED_EventRates NOvA_LED_eventrates.c myio.c NOvA_setup.h)

target_link_libraries(NOvA_SM_EventRates GSL::gsl globeslib)
target_link_libraries(NOvA_LED_EventRates GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(NOvA_SM_EventRates openblas)
     target_link_libraries(NOvA_LED_EventRates openblas)
endif()