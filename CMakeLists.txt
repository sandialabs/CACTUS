cmake_minimum_required(VERSION 3.11)

project(CACTUS
        VERSION 0.0
        DESCRIPTION "CACTUS"
        HOMEPAGE_URL http://www.github.com/SNL-WaterPower/CACTUS
        LANGUAGES Fortran)

# set compiler-specific options
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-O2 -fdefault-real-8 -ffree-line-length-0 -cpp")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-r8 -mkl -fpp")
    set(bounds "-O2 -check bounds")
    set(BLA_VENDOR Intel10_64lp)  # force use of MKL BLAS
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-r8 -O -Mvect=sse -Mcache_align -Mpre -Mnoframe -Mlre -Mflushz -Mautoinline")
    set(bounds "-C -Mbounds")
endif()

# find LAPACK, OpenMP
enable_language(C)  # For compatibility, some systems complain about this
find_package(LAPACK)

option(OPENMP "Compile with OpenMP enabled?" ON)
if(OPENMP MATCHES ON)
    find_package(OPENMP)
    if (NOT OpenMP_Fortran_FOUND)
        message(FATAL_ERROR "OPENMP option was set to true (default), but no OpenMP was found. Install OpenMP, or try again with -DOPENMP=OFF. Exiting.")
    endif()
    string(APPEND dialect " ${OpenMP_Fortran_FLAGS}")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

# add VERSION variable
# see https://cmake.org/pipermail/cmake/2018-October/068389.html
find_package(Git QUIET REQUIRED)

execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --abbrev=4 --dirty --always --tags
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE CACTUS_VERSION
    # ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

set_property(GLOBAL APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/.git/index")

add_definitions(-DVERSION="${CACTUS_VERSION}")

# create source lists
file(GLOB_RECURSE mod_sources ${PROJECT_SOURCE_DIR}/mod/**.f95)
file(GLOB_RECURSE prog_sources ${PROJECT_SOURCE_DIR}/src/*.f95)

add_executable(cactus ${prog_sources} ${mod_sources})
target_compile_options(cactus PUBLIC ${LAPACK_LINKER_FLAGS})
target_link_libraries(cactus ${LAPACK_LIBRARIES})

# install
install(TARGETS cactus
        RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin/
        )