module msis_wrapper
  ! Simple MSIS wrapper:
  ! - reads F10.7/AP from file 'F107_ap_appended.txt'
  ! - initialize MSIS via msisinit('msis21.parm')
  ! - call msis_point(year,month,day,hour,alt,lat,lon,slt, O,N2,O2,T)

  use msis_init, only : msisinit
  implicit none
  private
  public :: msis_wrapper_init, msis_point

  integer, allocatable :: r_iyd(:)    ! stored as year*1000 + doy
  integer, allocatable :: r_hour(:)
  real(4), allocatable :: r_ap(:), r_f107(:), r_f107a(:)
  integer :: nrec = 0
  logical :: loaded = .false.
  logical :: msis_inited = .false.

  ! Hard-coded values in case of issues reading the file in or no dates
  real(4), parameter :: DEFAULT_AP   = 10.0_4
  real(4), parameter :: DEFAULT_F107 = 150.0_4
  real(4), parameter :: DEFAULT_F107A= 150.0_4

  ! Explicit interface for external MSIS routine (must be in module spec part)
  interface
    subroutine gtd8d(iyd,ut,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
      integer, intent(in)      :: iyd
      real(4), intent(in)      :: ut          ! UT in hours
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
    if (.not. msis_inited) then
      call msisinit(parmpath='/discover/nobackup/jmpettit/GEOS_MLT_v7/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/@fvdycore/NRL_MSIS', parmfile='msis21.parm')
      msis_inited = .true.
    end if
    call load_f107_file()
    if (.not. loaded) then
      print *, 'Warning: F107/AP not loaded; using hard-coded defaults (ap=10,f107=150,f107a=150).'
    end if
  end subroutine msis_wrapper_init

  subroutine msis_point(year, month, day, hour, alt, glat, glong, stl, &
                        O_out, N2_out, O2_out, T_out)
    integer, intent(in) :: year, month, day, hour
    real(4), intent(in) :: alt, glat, glong, stl
    real(4), intent(out):: O_out, N2_out, O2_out, T_out

    integer :: doy, iyd, idx
    real(4) :: ut, ap(7), d(10), t(2)
    integer :: mass, i
    real(4) :: apv, f107v, f107av

    ! Compute day-of-year / iyd even if file not loaded (avoid uninitialized iyd)
    doy = day_of_year(year, month, day)
    iyd = year*1000 + doy

    ! Determine F107/AP values using F107_ap_appended.txt otherwise use defaults
    if (loaded) then
      idx = find_record(iyd, hour)
      if (idx >= 1) then
        apv = r_ap(idx); f107v = r_f107(idx); f107av = r_f107a(idx)
      else
        ! No matching date/hour found in file -> use defaults
        apv = DEFAULT_AP
        f107v = DEFAULT_F107
        f107av = DEFAULT_F107A
        print *, 'Warning: no F107/AP record found for', year, doy, hour, &
                 '; using defaults ap=',apv,' f107=',f107v,' f107a=',f107av
      end if
    else
      ! File not loaded -> use defaults (no error)
      apv = DEFAULT_AP
      f107v = DEFAULT_F107
      f107av = DEFAULT_F107A
    end if

    ! Prepare inputs and call MSIS: UT passed in hours
    ut = real(hour, kind=4)
    mass = 1
    do i = 1, 7
      ap(i) = apv
    end do

    call gtd8d(iyd, ut, alt, glat, glong, stl, f107av, f107v, ap, mass, d, t)

    O_out  = d(2)
    N2_out = d(3)
    O2_out = d(4)
    T_out  = t(2)
  end subroutine msis_point

  subroutine load_f107_file()
    ! Read fixed file: columns expected: year doy hour ap f107 f107a
    character(len=*), parameter :: fname = 'F107_ap_appended.txt'
    integer :: unit, ios, count
    character(len=256) :: line
    integer :: yy, doy, hr
    real(4) :: apv, f107v, f107av

    nrec = 0
    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      loaded = .false.
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

  subroutine shrink_arrays(new_n)
    integer, intent(in) :: new_n
    integer, allocatable :: tmpi(:)
    real(4), allocatable :: tmpf(:)

    tmpi = r_iyd(1:new_n); deallocate(r_iyd); allocate(r_iyd(new_n)); r_iyd = tmpi
    tmpi = r_hour(1:new_n); deallocate(r_hour); allocate(r_hour(new_n)); r_hour = tmpi

    tmpf = r_ap(1:new_n); deallocate(r_ap); allocate(r_ap(new_n)); r_ap = tmpf
    tmpf = r_f107(1:new_n); deallocate(r_f107); allocate(r_f107(new_n)); r_f107 = tmpf
    tmpf = r_f107a(1:new_n); deallocate(r_f107a); allocate(r_f107a(new_n)); r_f107a = tmpf
  end subroutine shrink_arrays

  integer function find_record(iyd, hour)
    integer, intent(in) :: iyd, hour
    integer :: i
    integer(kind=8) :: target, rec_time, best_diff
    integer :: best_idx

    if (nrec <= 0) then
      find_record = -1
      return
    end if

    ! try exact match first
    do i = 1, nrec
      if (r_iyd(i) == iyd .and. r_hour(i) == hour) then
        find_record = i
        return
      end if
    end do

    ! otherwise nearest in time (comparing iyd*24 + hour)
    target = int(iyd, kind=8) * 24_8 + int(hour, kind=8)
    best_diff = huge(0_8)
    best_idx = -1
    do i = 1, nrec
      rec_time = int(r_iyd(i), kind=8) * 24_8 + int(r_hour(i), kind=8)
      if (abs(rec_time - target) < best_diff) then
        best_diff = abs(rec_time - target)
        best_idx = i
      end if
    end do
    find_record = best_idx
  end function find_record

  integer function day_of_year(year, month, day) result(doy)
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
