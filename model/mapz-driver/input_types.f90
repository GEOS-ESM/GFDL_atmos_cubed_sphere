module input_types_mod

  implicit none
  
  private

  public ArrayExtents_T, Scalars_T, Arrays_T

  type ArrayExtents_T
     integer :: is, ie, isd, ied
     integer :: js, je, jsd, jed
     integer :: km
     integer :: nq
  end type ArrayExtents_T

  type Scalars_T
     ! All of these are inputs to Lagrangian_to_Eulerian
     logical :: last_step
     real :: consv
     real :: mdt
     real :: pdt
     integer :: nwat
     integer :: sphum
     real :: r_vir
     real :: cp
     real :: akap
     integer :: kord_mt
     integer :: kord_wz
     integer :: kord_tm
     integer :: ng
     logical :: fill
     logical :: reproduce_sum
     logical :: out_dt
     real :: ptop
     logical :: do_sat_adj
     logical :: hydrostatic
     logical :: hybrid_z
     logical :: do_omega
     logical :: adiabatic
     logical :: do_adiabatic_init
     integer :: remap_option
  end type Scalars_T

  type Arrays_T
     real, allocatable :: ps(:, :)
     real, allocatable :: pe(:, :, :)
     real, allocatable :: delp(:, :, :)
     real, allocatable :: pkz(:, :, :)
     real, allocatable :: pk(:, :, :)
     real, allocatable :: q_con(:, :, :)
     real, allocatable :: u(:, :, :) , v(:, :, :), w(:, :, :)
     real, allocatable :: delz(:, :, :)
     real, allocatable :: pt(:, :, :)
     real, allocatable :: q(:, :, :, :) ! this is the assumed shape array
     real, allocatable :: hs(:, :)
     real, allocatable :: cappa(:, :, :)
     integer, allocatable :: kord_tr(:)
     real, allocatable :: peln(:, :, :)
     real, allocatable :: te0_2d(:, :)
     real, allocatable :: ua(:, :, :), va(:, :, :)
     real, allocatable :: omga(:, :, :)
     real, allocatable :: te(:, :, :)
     real, allocatable :: ws(:, :)     
     real, allocatable :: dtdt(:, :, :)
     real, allocatable :: ak(:)
     real, allocatable :: bk(:)
     real, allocatable :: pfull(:)
     real, allocatable :: mfx(:, :, :), mfy(:, :, :)
     real, allocatable :: cx(:, :, :), cy(:, :, :)
  end type Arrays_T

end module input_types_mod
