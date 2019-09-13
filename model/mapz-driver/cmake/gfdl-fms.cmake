# GFDL FMS

# Include directory
if(NOT DEFINED ENV{FMS_INCLUDE_DIR})
  message(FATAL_ERROR "ERROR: env var FMS_INCLUDE_DIR is not defined")
else()
  message(STATUS "GFDL FMS include dir: $ENV{FMS_DIR}")
endif()
include_directories($ENV{FMS_INCLUDE_DIR})

# Libraries
if(NOT DEFINED ENV{FMS_LIB_DIR})
  message(FATAL_ERROR "ERROR: env var FMS_LIB_DIR is not defined")
else()
  message(STATUS "GFDL FMS lib dir: $ENV{FMS_LIB_DIR}")
endif()

find_library(FMS_LIBRARIES NAMES libGFDL_fms_r4.a PATHS $ENV{FMS_LIB_DIR})
message(STATUS "GFDL FMS: ${FMS_LIBRARIES}")
