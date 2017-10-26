#
# Minimum cmake version:
# Setting C++ standard seems to require at least 3.7
#
cmake_minimum_required(VERSION 3.7 FATAL_ERROR)

#
# Set project name and language
#
project(hbtplus)
enable_language(CXX)

#
# Set default build mode to Release
#
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (e.g. Debug, Release)" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

#
# Store library paths in executables
#
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#
# Put executables directly in build directory
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#
# Find OpenMP flag for this compiler, if necessary
#
option(HBT_USE_OPENMP "Whether to enable OpenMP support" ON)
if(HBT_USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_compile_options(${OpenMP_CXX_FLAGS})
    link_libraries(${OpenMP_CXX_FLAGS})
  else(OPENMP_FOUND)
    error("Unable to find OpenMP flag for this compiler")
  endif(OPENMP_FOUND)
endif()

#
# Find MPI libraries for C++
#
find_package(MPI REQUIRED)
add_definitions(${MPI_CXX_COMPILE_FLAGS})
include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${MPI_CXX_LIBRARIES})

#
# Find HDF5 library
#
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})

#
# Find GSL library (optional)
#
option(HBT_USE_GSL "Enable to use GSL for eigenvalue decomposition of inertial tensors" OFF)
if(HBT_USE_GSL)
  find_package(GSL)
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
    add_definitions(-DHAS_GSL)
  else(GSL_FOUND)
    error("Unable to find GSL library")
  endif(GSL_FOUND)
endif(HBT_USE_GSL)

#
# List of HBT source files.
# Headers are found automatically by cmake.
#
set(HBT_SRC 
  HBT.cpp
  src/config_parser.cpp
  src/halo.cpp
  src/linkedlist_parallel.cpp
  src/particle_exchanger.cpp
  src/subhalo.cpp
  src/subhalo_unbind.cpp
  src/geometric_tree.cpp
  src/hdf_wrapper.cpp
  src/mpi_wrapper.cpp
  src/snapshot.cpp
  src/subhalo_merge.cpp
  src/gravity_tree.cpp
  src/linkedlist.cpp
  src/mymath.cpp
  src/snapshot_exchanger.cpp
  src/subhalo_tracking.cpp
  src/io/apostle_io.cpp
  src/io/gadget_group_io.cpp
  src/io/gadget_io.cpp
  src/io/halo_io.cpp
  src/io/snapshot_io.cpp
  src/io/subhalo_io.cpp
  )

#
# Specify which C++ standard we need
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
# Function for defining optional features
#
function(hbt_option NAME DEFINE DESCRIPTION DEFAULT)
  option(${NAME} ${DESCRIPTION} ${DEFAULT})
  if(${NAME})
    add_definitions(-D${DEFINE})
  endif(${NAME})
endfunction()

#
# Optional features:
#
#          CMake option name           Preprocessor macro      Description                                       Default
hbt_option(HBT_DM_ONLY                 DM_ONLY                 "Enable to save memory in DMONLY runs"            OFF)
hbt_option(HBT_INPUT_INT8              INPUT_INT8              "Enable to read 8 byte integers from input files" ON)
hbt_option(HBT_INT8                    HBT_INT8                "Enable to use 8 byte integers"                   ON)
hbt_option(HBT_REAL8                   HBT_REAL8               "Enable to use 8 byte reals"                      OFF)
hbt_option(HBT_UNSIGNED_LONG_ID_OUTPUT UNSIGNED_LONG_ID_OUTPUT "Enable to output IDs as unsigned long"           ON)

#
# HBT executable
#
add_executable(HBT HBT.cpp ${HBT_SRC})
target_link_libraries(HBT ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES})

#
# Hide some variables we usually don't need from the default view in the cmake GUI
#
mark_as_advanced(HDF5_DIR)
mark_as_advanced(HDF5_CXX_LIBRARY_dl)
mark_as_advanced(HDF5_CXX_LIBRARY_hdf5)
mark_as_advanced(HDF5_CXX_LIBRARY_hdf5_cpp)
mark_as_advanced(HDF5_CXX_LIBRARY_hdf5_hl)
mark_as_advanced(HDF5_CXX_LIBRARY_hdf5_hl_cpp)
mark_as_advanced(HDF5_CXX_LIBRARY_m)
mark_as_advanced(HDF5_CXX_LIBRARY_rt)
mark_as_advanced(HDF5_CXX_LIBRARY_z)
mark_as_advanced(CMAKE_INSTALL_PREFIX)

#
# Show C++ compiler flags in default view
#
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_RELEASE)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_DEBUG)
