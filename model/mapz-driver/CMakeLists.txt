esma_set_this()

set(srcs
  ../fv_mapz.F90
  ../fv_fill.F90
  ../fv_arrays.F90
  ../fv_cmp.F90

  ../../tools/fv_timing.F90
  ../../tools/fv_mp_mod.F90

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

ecbuild_add_executable(
  TARGET mapz-driver
  SOURCES ${srcs}
  LIBS FMS::fms MAPL)
target_compile_definitions(mapz-driver PRIVATE MAPL_MODE SPMD TIMING)
target_compile_options(mapz-driver PRIVATE ${TRACEBACK})
set_target_properties(${this} PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
