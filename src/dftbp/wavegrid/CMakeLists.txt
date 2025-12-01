set(curdir ${CMAKE_CURRENT_SOURCE_DIR})


set(sources-fpp ${curdir}/wavegrid.F90)

add_subdirectory(basis)
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
        # Tell nvcc to pass -fopenmp to the host compiler (e.g. gcc)
        target_compile_options(wavegrid_cuda_objlib PRIVATE "-Xcompiler=-fopenmp")
        # Tell Cmake that if anyone uses this objlib, to also link OpenMP_CXX
        # (Alternatively, linking in OpenMP_CXX elsewhere would also work)
        target_link_libraries(wavegrid_cuda_objlib PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

