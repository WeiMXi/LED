add_executable(LEDProbability_NH LEDProbability_NH.c++)

target_link_libraries(LEDProbability_NH GSL::gsl globeslib LED_Lib)
if(OpenBLAS_FOUND)
     target_link_libraries(LEDProbability_NH openblas)
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T2K_Prob.glb
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/NOvA_Prob.glb
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})