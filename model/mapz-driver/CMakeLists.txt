cmake_minimum_required(VERSION 3.14.1)
project(mapz LANGUAGES Fortran)

set(srcs
  ../fv_mapz.F90
  ../fv_fill.F90
  ../fv_arrays.F90
  ../fv_cmp.F90
  ../lin_cloud_microphys.F90

  ../../tools/fv_timing.F90
  ../../tools/fv_mp_mod.F90

  stubs/fv_grid_utils_stub.f90
  stubs/external_ic_stub.F90
  
  input_types.f90
  input.f90
  output.f90
  main.F90)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Dependency - MPI
find_package(MPI REQUIRED)
message(STATUS "MPI Fortran libraries: ${MPI_Fortran_LIBRARIES}")

# Baselibs
find_package(Baselibs REQUIRED)

# Dependency - GFDL FMS
include(${CMAKE_CURRENT_LIST_DIR}/cmake/gfdl-fms.cmake)

add_executable(mapz ${srcs})
target_link_libraries(
  mapz
  ${FMS_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})
target_compile_definitions(mapz PRIVATE MAPL_MODE SPMD TIMING MAPZ_DRIVER_DEBUG)
target_compile_options(mapz PRIVATE -traceback)
