!#######################################################################
! MSIS® (NRL-SOF-014-1) SOFTWARE
! NRLMSIS® empirical atmospheric model software. Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! nrlmsis2.1_license.txt, which is part of this software package. BY
! USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
! CONDITIONS OF THE LICENSE.  
!#######################################################################

!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================

!==================================================================================================
! MSIS_GEOS: MSIS needed for GEOS-MLT
!==================================================================================================
program msis_geos

  use msis_init, only          : msisinit

  implicit none

  integer, parameter          :: nrec = 200

  integer                     :: iyd, mass, day, hour
  real(4)                     :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4)                     :: d(10),t(2)
  
  integer                     :: i
  character(128)              :: dummy

  !Initialize model
  call msisinit(parmpath='',parmfile='msis21.parm')

  !Open input and output files, loop through records, and call model
  open(77,file='F107_ap_appended.txt',status='old')
  open(78,file='msis2.1_geos.txt',status='replace')
  read(77,*) dummy
  write(78,'(10a7,3a13)') &
    'iyd','sec','alt','glat','glong','T','stl','f107a','f107','Ap','O','N2','O2'
  do i = 1,200
	  read(77,*) iyd,day,hour,apd,f107,f107a
    ap(1) = apd
    sec = hour*3600

!  Write call here to read in altitude, glat, glon, and stl prior to MSIS call
    
    alt=100
    glat=40.0
    glong=40.0
    stl=14.25

    call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)


!   Write in call here to return the O, O2, N2, and T back to GEOS-MLT

    write(78,'(2i7,4f7.1,f7.2,3f7.1,3e13.4)')  &
      iyd,int(sec),alt,glat,glong,t(2),stl,f107a,f107,ap(1),d(2:4)



  enddo
  close(77)
  close(78)

  stop

end program msis_geos
