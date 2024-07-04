cmake_minimum_required(VERSION 3.19...3.30)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Please use out of source build:
  cmake -B build")
endif()


project(MUMPS
        LANGUAGES C Fortran
        VERSION 5.7.2.0
        DESCRIPTION "Sparse direct parallel solver"
        HOMEPAGE_URL "http://mumps-solver.org/"
)