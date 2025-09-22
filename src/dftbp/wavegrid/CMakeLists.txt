set(curdir ${CMAKE_CURRENT_SOURCE_DIR})


set(sources-fpp
    ${curdir}/slater.F90
    ${curdir}/wavegrid.F90)

add_subdirectory(molorb)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)

# Cuda offloading support
if(WITH_CUDA)
    set(sources_lib_cuda
        molorb/cuda/kernel.cu
        molorb/cuda/host_logic.cu
    )
    add_library(wavegrid_cuda_objlib OBJECT ${sources_lib_cuda})

    if(WITH_OMP)
        target_compile_options(wavegrid_cuda_objlib PRIVATE "-Xcompiler=-fopenmp")
        target_link_libraries(wavegrid_cuda_objlib PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()




# Add the CUDA object files directly into dftbplus
