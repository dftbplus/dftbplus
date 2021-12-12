cmake_minimum_required(VERSION 3.16)

project(TestDftbPlusBuild LANGUAGES Fortran)

find_package(DftbPlus REQUIRED)

add_executable(test_build test_build.f90)
target_link_libraries(test_build DftbPlus::DftbPlus)


