cmake_minimum_required(VERSION 3.16)

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

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

# find LAPACK
# (this sets LAPACK_LINKER_FLAGS, LAPACK_LIBARIES, among others)
find_package(LAPACK)

# add VERSION variable
add_definitions(-DVERSION="SNAPSHOT")

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