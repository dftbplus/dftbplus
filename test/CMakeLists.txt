include(FetchContent)

if (TEST_AUTO_DOWNLOAD)
  FetchContent_Declare(
    slakos
    URL https://github.com/dftbplus/testparams/archive/1f5aefe45b7951f3fd1c0b63e40b80f80fdab612.tar.gz
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
