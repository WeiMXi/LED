add_executable(nu5d_mat_osc_NH nu5d_mat_osc_NH.c++)
add_executable(nu5d_mat_osc_IH nu5d_mat_osc_IH.c++)

target_link_libraries(nu5d_mat_osc_NH GSL::gsl globeslib LED_Lib)
target_link_libraries(nu5d_mat_osc_IH GSL::gsl globeslib LED_Lib)
if(OpenBLAS_FOUND)
     target_link_libraries(nu5d_mat_osc_NH openblas)
endif()

if(OpenBLAS_FOUND)
     target_link_libraries(nu5d_mat_osc_IH openblas)
endif()


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T2K_Prob.glb
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/NOvA_Prob.glb
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})