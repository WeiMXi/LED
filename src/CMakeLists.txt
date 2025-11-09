file(GLOB_RECURSE LED_Lib_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/*.c++ ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/*.h++)
add_library(LED_Lib STATIC ${LED_Lib_SRC})
target_link_libraries(LED_Lib GSL::gsl globeslib)
if(OpenBLAS_FOUND)
     target_link_libraries(LED_Lib openblas)
endif()
target_include_directories(LED_Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


add_subdirectory(5d_bulk)
add_subdirectory(NOvA_2108.08219)
add_subdirectory(T2K_2303.03222)