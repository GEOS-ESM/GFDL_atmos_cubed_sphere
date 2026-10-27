esma_set_this()

set(srcs
  ../calc_gas_specific_heat_MLT.F90
  ../fv_mapz.F90
  ../fv_fill.F90
  ../fv_arrays.F90
  ../fv_cmp.F90
  ../lin_cloud_microphys.F90

  ../../tools/fv_timing.F90
  ../../tools/fv_mp_mod.F90

  ../../NRL_MSIS/msis_constants.F90
  ../../NRL_MSIS/msis_utils.F90
  ../../NRL_MSIS/msis_init.F90
  ../../NRL_MSIS/msis_gfn.F90
  ../../NRL_MSIS/msis_tfn.F90
  ../../NRL_MSIS/msis_dfn.F90
  ../../NRL_MSIS/msis_calc.F90
  ../../NRL_MSIS/msis_gtd8d.F90
  ../../NRL_MSIS/msis_wrapper.F90

  stubs/fv_grid_utils_stub.f90
  stubs/external_ic_stub.F90
  
  input/input_array_extents.f90
  input/input_scalars.f90
  input/input_arrays.f90
  input/input_flagstruct.f90
  input/input_gridstruct.f90
  input/input_domain_layout.f90
  output/output.f90
  main.F90)

if (FV_PRECISION STREQUAL R4)
  set(FMS fms_r4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
  set(FMS fms_r8)
elseif (FV_PRECISION STREQUAL R8)
   set(FMS fms_r8)
endif ()

ecbuild_add_executable(
  TARGET mapz-driver
  SOURCES ${srcs}
  LIBS ${FMS} MAPL)
target_compile_definitions(mapz-driver PRIVATE MAPL_MODE SPMD TIMING)
target_compile_options(mapz-driver PRIVATE ${TRACEBACK})
set_target_properties(${this} PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
