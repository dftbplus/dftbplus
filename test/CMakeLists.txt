include(FetchContent)

if (TEST_AUTO_DOWNLOAD)
  FetchContent_Declare(
    slakos
    URL https://github.com/dftbplus/testparams/archive/0e1d95abc70339d017fd53d40c41c041057fd03f.tar.gz
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/slakos/origin"
  )
  FetchContent_MakeAvailable(slakos)

  FetchContent_Declare(
    gbsa
    URL https://github.com/grimme-lab/gbsa-parameters/archive/6836c4d997e4135e418cfbe273c96b1a3adb13e2.tar.gz
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/gbsa/origin"
  )
  FetchContent_MakeAvailable(gbsa)
endif()

add_subdirectory(src/dftbp)
add_subdirectory(app)
add_subdirectory(tools/dptools)
