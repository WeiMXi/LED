add_executable(nu5d_mat_osc nu5d_mat_osc.c++)

target_link_libraries(nu5d_mat_osc GSL::gsl globeslib LED_Lib)
if(OpenBLAS_FOUND)
     target_link_libraries(nu5d_mat_osc openblas)
endif()


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test.glb
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})