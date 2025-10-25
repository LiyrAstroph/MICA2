cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(mica2 LANGUAGES C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra -finline-functions -fcommon -ffast-math")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")

##if using clang compiler, add c++11/c++14/c++17... support (https://clang.llvm.org/cxx_status.html)
message("C Compiler: " ${CMAKE_C_COMPILER} "," ${CMAKE_C_COMPILER_ID} "," ${CMAKE_C_STANDARD})
if (${CMAKE_C_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang")
  set(CMAKE_C_FLAGS "-std=c11 ${CMAKE_C_FLAGS}")
endif()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# when installing directly in CMake, redirect its intall root to dist dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "default install path" FORCE )
endif()

set(SRC "src")
set(CDNEST "cdnest")
include_directories(${CDNEST})
include_directories(${SRC})

add_executable(${CMAKE_PROJECT_NAME}
  ${SRC}/main.c
  ${SRC}/allvars.c
  ${SRC}/run.c 
  ${SRC}/dnest_con.c
  ${SRC}/dnest_line.c 
  ${SRC}/dnest_pmap.c
  ${SRC}/dnest_vmap.c
  ${SRC}/dnest_mmap.c
  ${SRC}/read.c
  ${SRC}/mc_con.c
  ${SRC}/mc_line.c
  ${SRC}/mc_pmap.c
  ${SRC}/mc_vmap.c
  ${SRC}/mc_mmap.c
  ${SRC}/init.c
  ${SRC}/mathfun.c
  ${SRC}/line_decomp.c
  ${SRC}/line_decomp_mmap.c
  ${SRC}/system.c
  ${SRC}/version.c
  ${SRC}/error.c
  ${SRC}/mpifunc.c
  ${SRC}/mygetopt.c
  ${CDNEST}/dnest.c
  ${CDNEST}/dnestvars.c
  ${CDNEST}/dnestpostprocess.c
  )

# Math library
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC m)

#=====================================================
# MPI
find_package(MPI REQUIRED)
if(MPI_FOUND)
  include_directories(${MPI_C_INCLUDE_DIRS})
  target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${MPI_C_LIBRARIES})
else(MPI_FOUND)
message(FATAL_ERROR "MPI not found")
endif(MPI_FOUND)

#=====================================================
# GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL not found")
endif(GSL_FOUND)


find_package(MKL CONFIG)
if(MKL_FOUND)
  #message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
  target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE -DIntelMKL)
  target_compile_options(${CMAKE_PROJECT_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::MKL>)

else(MKL_FOUND)
  #=====================================================
  # LAPACK
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
    target_link_options(${CMAKE_PROJECT_NAME} PUBLIC ${LAPACK_LINKER_FLAGS})
    target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} lapacke)
  else(LAPACK_FOUND)
    message(FATAL_ERROR "LAPACKE not found")
  endif(LAPACK_FOUND)
endif(MKL_FOUND)

