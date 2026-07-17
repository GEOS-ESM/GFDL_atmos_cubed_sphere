module msis_wrapper
  ! Simple MSIS wrapper:
  ! - reads F10.7/AP from file 'F107_ap_appended.txt'
  ! - initialize MSIS via msisinit('msis21.parm')
  ! - call msis_point(year,doy,ut_seconds,alt,lat,lon,slt, O,N2,O2,T)

  use msis_init, only : msisinit

  implicit none
  
  private
  
  public :: msis_wrapper_init, msis_prepare_time, msis_point

  ! --- Local Variables ---
  integer, allocatable :: r_iyd(:)    ! stored as year*1000 + doy
  integer, allocatable :: r_hour(:)
  real(4), allocatable :: r_ap(:), r_f107(:), r_f107a(:)
  integer :: nrec = 0
  logical :: loaded = .false.
  logical :: msis_inited = .false.

  ! --- Hard-coded values in case of issues reading the file in or no dates ---
  real(4), parameter :: DEFAULT_AP   = 10.0_4
  real(4), parameter :: DEFAULT_F107 = 150.0_4
  real(4), parameter :: DEFAULT_F107A= 150.0_4

  ! --- Cached space-weather state ---
  ! Only the hourly forcing indices are cached. MSIS itself is still called
  ! for every requested time, location, and altitude.
  logical, save :: index_cache_ready = .false.
  integer, save :: cached_iyd  = -1
  integer, save :: cached_hour = -1
  real(4), save :: cached_ap    = DEFAULT_AP
  real(4), save :: cached_f107  = DEFAULT_F107
  real(4), save :: cached_f107a = DEFAULT_F107A

  ! Enable only for a short verification run. Output is produced by each MPI rank.
  logical, parameter :: DEBUG_MSIS_INDICES = .false.

  ! --- Explicit interface for external MSIS routine ---
  interface
    subroutine gtd8d(iyd,ut,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
      integer, intent(in)      :: iyd
      real(4), intent(in)      :: ut          ! UT in seconds
      real(4), intent(in)      :: alt, glat, glong, stl, f107a, f107
      real(4), intent(in)      :: ap(7)
      integer, intent(in)      :: mass
      real(4), intent(out)     :: d(10)
      real(4), intent(out)     :: t(2)
    end subroutine gtd8d
  end interface

contains

  subroutine msis_wrapper_init()
    ! Initialize MSIS and load f107_ap_appended.txt
    if (msis_inited .and. loaded) return

    if (.not. msis_inited) then
      call msisinit(parmpath='./', parmfile='msis21.parm')
      msis_inited = .true.
    end if

    call load_f107_file()
    
    if (.not. loaded) then
      print *, 'GEOS_MLT_MSIS_ERROR: F107/AP table was not loaded.'
      error stop 'GEOS-MLT MSIS initialization failed'
    end if

    index_cache_ready = .false.
  end subroutine msis_wrapper_init


  subroutine msis_prepare_time(year, doy, ut_seconds)
    ! Select the time-varying space-weather indices once per model hour.
    ! Subsequent msis_point calls reuse these scalar indices but still run MSIS.
    integer, intent(in) :: year, doy, ut_seconds

    integer :: requested_iyd, requested_hour, idx

    requested_iyd  = year*1000 + doy
    requested_hour = max(0, min(23, ut_seconds / 3600))

    if (index_cache_ready) then
      if (requested_iyd  == cached_iyd .and. &
          requested_hour == cached_hour) return
    end if

    if (.not. loaded) then
      print *, 'GEOS_MLT_MSIS_ERROR: index table is not loaded.'
      print *, 'Requested year,doy,hour=', year, doy, requested_hour
      error stop 'GEOS-MLT MSIS table was not initialized'
    end if

    ! This linear table search now occurs only when the model hour changes.
    idx = find_record(requested_iyd, requested_hour)

    if (idx < 1) then
      print *, 'GEOS_MLT_MSIS_ERROR: no exact index record.'
      print *, 'Requested year,doy,hour=', year, doy, requested_hour
      error stop 'GEOS-MLT missing F107/AP record'
    end if

    cached_iyd   = requested_iyd
    cached_hour  = requested_hour
    cached_ap    = r_ap(idx)
    cached_f107  = r_f107(idx)
    cached_f107a = r_f107a(idx)
    index_cache_ready = .true.

    if (DEBUG_MSIS_INDICES) then
      write(*,'(A,3(I0,1X),A,3(F10.3,1X))') &
           'GEOS_MLT_MSIS_INDEX year,doy,hour=', &
           year, doy, requested_hour, ' Ap,F107,F107A=', &
           cached_ap, cached_f107, cached_f107a
    end if
  end subroutine msis_prepare_time

  subroutine msis_point(year, doy, ut_seconds, alt, glat, glong, stl, &
                        O_out, N2_out, O2_out, T_out)
    integer, intent(in) :: year, doy, ut_seconds
    real(4), intent(in) :: alt, glat, glong, stl
    real(4), intent(out):: O_out, N2_out, O2_out, T_out

    integer :: iyd, hour
    real(4) :: ut, ap(7), d(10), t(2)
    integer :: mass

    ! Compute day-of-year / iyd even if file not loaded 
    iyd = year*1000 + doy
    hour = max(0, min(23, ut_seconds / 3600))

    ! msis_prepare_time must be called before the grid-point/level MSIS loops.
    if (.not. index_cache_ready) then
      error stop 'Call msis_prepare_time before msis_point'
    end if

    ! Do not allow indices from a previous model hour to be used accidentally.
    if (iyd /= cached_iyd .or. hour /= cached_hour) then
      print *, 'GEOS_MLT_MSIS_ERROR: stale index cache.'
      print *, 'Requested iyd,hour=', iyd, hour
      print *, 'Cached    iyd,hour=', cached_iyd, cached_hour
      error stop 'GEOS-MLT stale MSIS index cache'
    end if

    ! Prepare inputs and call MSIS: UT passed in seconds
    ut = real(ut_seconds, kind=4)
    mass = 1
    ap(:) = cached_ap
    
! Logging
    !print *, 'F107, ap, stl, alt:', cached_f107, ap, stl, alt
    !print *, iyd, ut
    call gtd8d(iyd, ut, alt, glat, glong, stl, &
               cached_f107a, cached_f107, ap, mass, d, t)

    !print *, 'MSIS Temperature:', t(2)
    O_out  = d(2)
    N2_out = d(3)
    O2_out = d(4)
    T_out  = t(2)
  end subroutine msis_point


  subroutine load_f107_file()
    ! Read inpute file (F107_ap_appended.txt): columns: year doy hour ap f107 f107a
    character(len=*), parameter :: fname = 'F107_ap_appended.txt'
    integer :: unit, ios, count
    character(len=256) :: line
    integer :: yy, doy, hr
    real(4) :: apv, f107v, f107av

    nrec = 0

    !print *, 'Trying to open file: ', fname
    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    !print *, 'open iostat = ', ios
    if (ios /= 0) then
       loaded = .false.
       print *, 'Failed to open file.'
       return
    end if

    ! Count non-empty lines
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      nrec = nrec + 1
    end do
    close(unit)

    if (nrec <= 0) then
      loaded = .false.
      return
    end if

    allocate(r_iyd(nrec)); allocate(r_hour(nrec))
    allocate(r_ap(nrec)); allocate(r_f107(nrec)); allocate(r_f107a(nrec))

    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      loaded = .false.
      return
    end if

    count = 0
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read(line,*,iostat=ios) yy, doy, hr, apv, f107v, f107av
      if (ios /= 0) then
        cycle
      end if
      count = count + 1
      r_iyd(count)   = yy*1000 + doy
      r_hour(count)  = hr
      r_ap(count)    = apv
      r_f107(count)  = f107v
      r_f107a(count) = f107av
    end do
    close(unit)

    if (count < nrec) then
      nrec = count
      call shrink_arrays(nrec)
    end if

    if (nrec == 0) then
      loaded = .false.
      return
    end if

    loaded = .true.
  end subroutine load_f107_file


  ! I used AI to help write this subroutine. It should save on memory rather than hold entire file in memory (~6mb).
  subroutine shrink_arrays(new_n)
     integer, intent(in) :: new_n
     integer, allocatable :: tmpi(:)
     real(4), allocatable :: tmpf(:)

     allocate(tmpi(new_n))
     tmpi = r_iyd(1:new_n)
     deallocate(r_iyd)
     allocate(r_iyd(new_n))
     r_iyd = tmpi
     deallocate(tmpi)

     allocate(tmpi(new_n))
     tmpi = r_hour(1:new_n)
     deallocate(r_hour)
     allocate(r_hour(new_n))
     r_hour = tmpi
     deallocate(tmpi)

     allocate(tmpf(new_n))
     tmpf = r_ap(1:new_n)
     deallocate(r_ap)
     allocate(r_ap(new_n))
     r_ap = tmpf
     deallocate(tmpf)

     allocate(tmpf(new_n))
     tmpf = r_f107(1:new_n)
     deallocate(r_f107)
     allocate(r_f107(new_n))
     r_f107 = tmpf
     deallocate(tmpf)

     allocate(tmpf(new_n))
     tmpf = r_f107a(1:new_n)
     deallocate(r_f107a)
     allocate(r_f107a(new_n))
     r_f107a = tmpf
     deallocate(tmpf)
  end subroutine shrink_arrays

  integer function find_record(iyd, hour)
    integer, intent(in) :: iyd, hour
    integer :: i

    find_record = -1

    if (nrec <= 0) return

    ! Find exact day of year and record
    do i = 1, nrec
      if (r_iyd(i) == iyd .and. r_hour(i) == hour) then
        find_record = i
        return
      end if
    end do

  end function find_record

  function day_of_year(year, month, day) result(doy)
    
    integer :: doy
    integer, intent(in) :: year, month, day
    integer :: mdays(12), m
    logical :: leap

    ! month lengths
    mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    leap = .false.
    if (mod(year,400) == 0 .or. (mod(year,4) == 0 .and. mod(year,100) /= 0)) leap = .true.
    if (leap) mdays(2) = 29

    doy = 0
    do m = 1, month-1
      doy = doy + mdays(m)
    end do
    doy = doy + day
  end function day_of_year

end module msis_wrapper
