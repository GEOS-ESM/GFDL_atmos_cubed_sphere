# Kokkos
if(NOT EXISTS ${KOKKOS_SRC_DIR})
  message(FATAL_ERROR "ERROR: Must specify KOKKOS_SRC_DIR with cmake ... -DKOKKOS_SRC_DIR=<path>.")
endif ()

add_subdirectory(${KOKKOS_SRC_DIR} ${CMAKE_BINARY_DIR}/kokkos)

include_directories(${Kokkos_INCLUDE_DIRS_RET})
