# Kokkos
if(NOT DEFINED ENV{KOKKOS_SRC_DIR})
  message(FATAL_ERROR "ERROR: env var KOKKOS_SRC_DIR is not defined")
endif()

add_subdirectory($ENV{KOKKOS_SRC_DIR} ${CMAKE_CURRENT_BINARY_DIR}/kokkos)

include_directories(${Kokkos_INCLUDE_DIRS_RET})
