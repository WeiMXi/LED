########################################
##### LED Library CMakeLists.txt #######
file(GLOB_RECURSE LED_Lib_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/*.c++
    ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/*.h++
)

list(REMOVE_ITEM LED_Lib_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/Engine/CUDA/CUDAHermDiag.h++
    ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/Engine/CUDA/CUDAHermDiag.cu
)

if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    file(GLOB LED_Lib_CUDA_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/Engine/CUDA/CUDAHermDiag.cu
        ${CMAKE_CURRENT_SOURCE_DIR}/ledlib/Engine/CUDA/CUDAHermDiag.h++
    )
    list(APPEND LED_Lib_SRC ${LED_Lib_CUDA_SRCS})
endif()

add_library(LED_Lib STATIC ${LED_Lib_SRC})

# include dirs
target_include_directories(LED_Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${GSL_INCLUDE_DIRS})

# link GSL and your other libs
target_link_libraries(LED_Lib PUBLIC GSL::gsl globeslib)

if(OpenBLAS_FOUND)
    target_link_libraries(LED_Lib PUBLIC openblas)
endif()


if(CUDAToolkit_FOUND)
    target_link_libraries(LED_Lib PUBLIC CUDA::cudart CUDA::cublas CUDA::cusolver)
    set_target_properties(LED_Lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_ARCHITECTURES "native")
    set_target_properties(LED_Lib PROPERTIES
        CUDA_ARCHITECTURES "${CUDA_ARCHITECTURES}"
    )
    target_compile_definitions(LED_Lib PUBLIC HAVE_CUDA)
endif()
########################################


add_subdirectory(probability)
add_subdirectory(NOvA_2108.08219)
add_subdirectory(T2K_2303.03222)