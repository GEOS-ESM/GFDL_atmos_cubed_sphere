# Kokkos
if(NOT EXISTS ${KOKKOS_PATH})
  message(FATAL_ERROR "ERROR: Must specify KOKKOS_PATH with cmake ... -DKOKKOS_PATH=<path>.")
endif ()

add_subdirectory(${KOKKOS_PATH} ${CMAKE_BINARY_DIR}/kokkos)

include_directories(${Kokkos_INCLUDE_DIRS_RET})
