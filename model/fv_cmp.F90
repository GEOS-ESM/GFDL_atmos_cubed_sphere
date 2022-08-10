!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the GFDL Cloud Microphysics.
!*
!* The GFDL Cloud Microphysics is free software: you can 
!8 redistribute it and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The GFDL Cloud Microphysics is distributed in the hope it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the GFDL Cloud Microphysics.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'fv_cmp' implements the fast procesesses in the GFDL
!! microphysics
!!>@author Shian-Jiann Lin, Linjiong Zhou
! Fast saturation adjustment is part of the gfdl cloud microphysics
! =======================================================================

module fv_cmp_mod
! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>rvgas, rdgas, grav, hlv, hlf, cp_air</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td> r_grid</td>
!   </tr>
!   <tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is_master</td>
!   </tr>
!   <tr>
!     <td>gfdl_lin_cloud_microphys_mod</td>
!     <td>ql_gen, qi_gen, qi0_max, ql_mlt, ql0_max, qi_lim, qs_mlt,
!         tau_r2g, tau_smlt, tau_i2s, tau_v2l, tau_l2v, tau_imlt, tau_l2r,
!         preciprad, dw_ocean, dw_land</td>
!   </tr>
! </table>
    
    use constants_mod, only: rvgas, rdgas, grav, hlv, hlf, cp_air
    use fv_mp_mod, only: is_master
    use fv_arrays_mod, only: r_grid
    use gfdl_lin_cloud_microphys_mod, only: ql_gen, qi_gen, qi0_crt, qi0_max, ql_mlt, ql0_max, qi_lim, qs_mlt
    use gfdl_lin_cloud_microphys_mod, only: icloud_f, sat_adj0, t_sub, cld_min, dt_fr
    use gfdl_lin_cloud_microphys_mod, only: tau_r2g, tau_smlt, tau_i2s, tau_v2l, tau_l2v, tau_i2v, tau_imlt, tau_l2r, tau_frz
    use gfdl_lin_cloud_microphys_mod, only: preciprad, dw_ocean, dw_land, do_qa
    
    implicit none
    
    private
    
    public fv_sat_adj, qs_init
    
    ! real, parameter :: cp_air = cp_air ! 1004.6, heat capacity of dry air at constant pressure, come from constants_mod
    real, parameter :: cp_vap = 4.0 * rvgas !< 1846.0, heat capacity of water vapor at constant pressure
    real, parameter :: cv_air = cp_air - rdgas !< 717.55, heat capacity of dry air at constant volume
    real, parameter :: cv_vap = 3.0 * rvgas !< 1384.5, heat capacity of water vapor at constant volume
    real, parameter :: pi = 3.1415926535897931 !< gfs: ratio of circle circumference to diameter
    
    ! http: / / www.engineeringtoolbox.com / ice - thermal - properties - d_576.html
    ! c_ice = 2050.0 at 0 deg c
    ! c_ice = 1972.0 at - 15 deg c
    ! c_ice = 1818.0 at - 40 deg c
    ! http: / / www.engineeringtoolbox.com / water - thermal - properties - d_162.html
    ! c_liq = 4205.0 at 4 deg c
    ! c_liq = 4185.5 at 15 deg c
    ! c_liq = 4178.0 at 30 deg c
    
    ! real, parameter :: c_ice = 2106.0 ! ifs: heat capacity of ice at 0 deg c
    ! real, parameter :: c_liq = 4218.0 ! ifs: heat capacity of liquid at 0 deg c
    real, parameter :: c_ice = 1972.0 !< gfdl: heat capacity of ice at - 15 deg c
    real, parameter :: c_liq = 4185.5 !< gfdl: heat capacity of liquid at 15 deg c
   

    real, parameter :: qpmin = 1.e-8 ! min value for suspended rain/snow/liquid/ice condensate
    real, parameter :: qvmin = 1.e-20 !< min value for water vapor (treated as zero)
    real, parameter :: qcmin = 1.e-12 !< min value for cloud condensates
 
    real, parameter :: dc_vap = cp_vap - c_liq !< - 2339.5, isobaric heating / cooling
    real, parameter :: dc_ice = c_liq - c_ice !< 2213.5, isobaric heating / colling
    
    real, parameter :: tice = 273.16 !< freezing temperature
    real, parameter :: t_wfr = tice - 40. !< homogeneous freezing temperature

    real, parameter :: t_evap = 273.16 + 10.0
 
    real, parameter :: lv0 = hlv - dc_vap * tice !< 3.13905782e6, evaporation latent heat coefficient at 0 deg k
    real, parameter :: li00 = hlf - dc_ice * tice !< - 2.7105966e5, fusion latent heat coefficient at 0 deg k
    
    ! real (kind = r_grid), parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c
    real (kind = r_grid), parameter :: e00 = 611.21 !< ifs: saturation vapor pressure at 0 deg c
    
    real (kind = r_grid), parameter :: d2ice = dc_vap + dc_ice !< - 126, isobaric heating / cooling
    real (kind = r_grid), parameter :: li2 = lv0 + li00 !< 2.86799816e6, sublimation latent heat coefficient at 0 deg k
    
    real, parameter :: lat2 = (hlv + hlf) ** 2 !< used in bigg mechanism
    
    real :: d0_vap !< the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
    real :: lv00 !< the same as lv0, except that cp_vap can be cp_vap or cv_vap
    
    real, allocatable :: table (:), table2 (:), tablew (:), des2 (:), desw (:)
    
    logical :: mp_initialized = .false.
    
contains


!>@brief The subroutine 'fv_sat_adj' performs the fast processes in the GFDL microphysics.
!>@details This is designed for single-moment 6-class cloud microphysics schemes.
!! It handles the heat release due to in situ phase changes.
subroutine fv_sat_adj (mdt, zvir, is, ie, js, je, ng, hydrostatic, consv_te, &
        te0, qv, ql, qi, qr, qs, qg, hs, dpln, delz, pt, dp, cappa, &
        area, dtdt, out_dt, last_step, qa)
 
    implicit none
    
    integer, intent (in) :: is, ie, js, je, ng
    
    logical, intent (in) :: hydrostatic, consv_te, out_dt, last_step
    
    real, intent (in) :: zvir, mdt ! remapping time step
    
    real, intent (in), dimension (is - ng:ie + ng, js - ng:je + ng) :: dp, delz, hs
    real, intent (in), dimension (is:ie, js:je) :: dpln
   
    real, intent (inout), dimension (is - ng:ie + ng, js - ng:je + ng) :: pt, qv, ql, qi, qr, qs, qg
    real, intent (inout), dimension (is - ng:, js - ng:) :: cappa
    real, intent (inout), dimension (is:ie, js:je) :: dtdt
    
    real, intent (out), dimension (is - ng:ie + ng, js - ng:je + ng) :: qa, te0
    
    real (kind = r_grid), intent (in), dimension (is - ng:ie + ng, js - ng:je + ng) :: area
    
    real, dimension (is:ie) :: wqsat, dq2dt, qpz, cvm, t0, pt1, qstar
    real, dimension (is:ie) :: icp2, lcp2, tcp2, tcp3
    real, dimension (is:ie) :: den, q_liq, q_sol, q_cond, evap, src, sink, hvar
    real, dimension (is:ie) :: mc_air, lhl, lhi
    
    real :: qsw, rh
    real :: tc, qsi, dqsdt, dq, dq0, pidep, qi_crt, tmp, dtmp
    real :: tin, rqi, q_plus, q_minus, dqh
    real :: sdt, dt_bigg
    real :: fac_smlt, fac_r2g, fac_i2s, fac_imlt, fac_l2r, fac_v2l, fac_l2v, fac_i2v, fac_frz
    real :: factor, qim, tice0, c_air, c_vap, dw
    real :: a1, newqi, newql

    integer :: i, j
    
    sdt = 0.5 * mdt ! half remapping time step
    dt_bigg = mdt ! bigg mechinism time step
    
    tice0 = tice - 0.01 ! 273.15, standard freezing temperature
    
    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------
    
    fac_i2s = sat_adj0 * (1. - exp (- mdt / tau_i2s))
    fac_v2l = sat_adj0 * (1. - exp (- sdt / tau_v2l))
    fac_r2g = sat_adj0 * (1. - exp (- mdt / tau_r2g))
    fac_l2r = sat_adj0 * (1. - exp (- mdt / tau_l2r))
    
    fac_l2v = sat_adj0 * (1. - exp (- sdt / tau_l2v))
    fac_i2v = sat_adj0 * (1. - exp (- sdt / tau_i2v))

    fac_imlt = sat_adj0 * (1. - exp (- sdt / tau_imlt))
    fac_smlt = sat_adj0 * (1. - exp (- mdt / tau_smlt))
    
    fac_frz = sat_adj0 * (1. - exp (- mdt / tau_frz))
    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------
    
    if (hydrostatic) then
        c_air = cp_air
        c_vap = cp_vap
    else
        c_air = cv_air
        c_vap = cv_vap
    endif
    d0_vap = c_vap - c_liq
    lv00 = hlv - d0_vap * tice
    ! dc_vap = cp_vap - c_liq ! - 2339.5
    ! d0_vap = cv_vap - c_liq ! - 2801.0
    
    do j = js, je ! start j loop

      ! Compute subgrid variance (hvar)
        do i = is, ie
            ! higher than 10 m is considered "land" and will have higher subgrid variability
            dw = dw_ocean + (dw_land - dw_ocean) * min (1., abs (hs (i, j)) / (10. * grav))
            ! "scale - aware" subgrid variability: 100 - km as the base
            hvar (i) = min (0.3, max (0.01, dw * sqrt (sqrt (area (i, j)) / 100.e3)))
        enddo
        
        do i = is, ie
            q_liq (i) = ql (i, j) + qr (i, j)
            q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
            qpz (i) = q_liq (i) + q_sol (i)
#ifdef USE_COND
            pt1 (i) = pt (i, j) / ((1 + zvir * qv (i, j)) * (1 - qpz (i)))
#else
            pt1 (i) = pt (i, j) / (1 + zvir * qv (i, j))
#endif
            t0 (i) = pt1 (i) ! true temperature
            qpz (i) = qpz (i) + qv (i, j) ! total_wat conserved in this routine
        enddo
        
        ! -----------------------------------------------------------------------
        ! define air density based on hydrostatical property
        ! -----------------------------------------------------------------------
        
        if (hydrostatic) then
            do i = is, ie
                den (i) = dp (i, j) / (dpln (i, j) * rdgas * pt (i, j))
            enddo
        else
            do i = is, ie
                den (i) = - dp (i, j) / (grav * delz (i, j)) ! moist_air density
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! define heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            mc_air (i) = (1. - qpz (i)) * c_air ! constant
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            if (hydrostatic) then
                do i = is, ie
                    te0 (i, j) = - c_air * t0 (i)
                enddo
            else
                do i = is, ie
#ifdef USE_COND
                    te0 (i, j) = - cvm (i) * t0 (i)
#else
                    te0 (i, j) = - c_air * t0 (i)
#endif
                enddo
            endif
        endif
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud ice with snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qi (i, j) < 0.) then
                qs (i, j) = qs (i, j) + qi (i, j)
                qi (i, j) = 0.
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of cloud ice to cloud water and rain
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qi (i, j) > qcmin .and. pt1 (i) > tice) then
                sink (i) = min (qi (i, j), fac_imlt * (pt1 (i) - tice) / icp2 (i))
                qi (i, j) = qi (i, j) - sink (i)
                ! sjl, may 17, 2017
                ! tmp = min (sink (i), dim (ql_mlt, ql (i, j))) ! max ql amount
                ! ql (i, j) = ql (i, j) + tmp
                ! qr (i, j) = qr (i, j) + sink (i) - tmp
                ! sjl, may 17, 2017
                ql (i, j) = ql (i, j) + sink (i)
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative snow with graupel or graupel with available snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qs (i, j) < 0.) then
                qg (i, j) = qg (i, j) + qs (i, j)
                qs (i, j) = 0.
            elseif (qg (i, j) < 0.) then
                tmp = min (- qg (i, j), max (0., qs (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qs (i, j) = qs (i, j) - tmp
            endif
        enddo
        
        ! after this point cloud ice & snow are positive definite
        
        ! -----------------------------------------------------------------------
        ! fix negative cloud water with rain or rain with available cloud water
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) < 0.) then
                tmp = min (- ql (i, j), max (0., qr (i, j)))
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) - tmp
            elseif (qr (i, j) < 0.) then
                tmp = min (- qr (i, j), max (0., ql (i, j)))
                ql (i, j) = ql (i, j) - tmp
                qr (i, j) = qr (i, j) + tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! enforce complete freezing of cloud water to cloud ice below - 48 c
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = tice - 48. - pt1 (i)
            if (ql (i, j) > 0. .and. dtmp > 0.) then
                sink (i) = min (ql (i, j), dtmp / icp2 (i))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
        enddo
        
        ! -----------------------------------------------------------------------
        ! condensation / evaporation between water vapor and cloud water
        ! -----------------------------------------------------------------------
        if (do_qa) then 
          call wqs2_vect (is, ie, pt1, den, wqsat, dq2dt)
          do i = is, ie
            dq0 = wqsat(i) - qv(i,j)
            factor = min (1., fac_l2v * (10. * dq0 / wqsat(i)))
            evap (i) = min (ql (i,j), factor * ql(i,j) / (1. + tcp3 (i) * dq2dt (i)))
            qv (i, j) = qv (i, j) + evap (i)
            ql (i, j) = ql (i, j) - evap (i)
            q_liq (i) = q_liq (i) - evap (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) - evap (i) * lhl (i) / cvm (i)
          enddo
          ! -----------------------------------------------------------------------
          ! update latend heat coefficient
          ! -----------------------------------------------------------------------
          do i = is, ie
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lhi (i) = li00 + dc_ice * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp3 (i) = lcp2 (i) + icp2 (i) * min (1., dim (tice, pt1 (i)) / 48.)
          enddo
        endif

        ! -----------------------------------------------------------------------
        ! homogeneous freezing of cloud water to cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = (t_wfr - dt_fr) - pt1 (i) ! [ - 40, - 48]
            if (ql (i, j) > 0. .and. dtmp > 0.) then
                newqi = new_ice_condensate(pt1 (i), ql (i,j), qi (i,j))
                sink (i) = max(0.0,min (newqi, fac_frz * dtmp / icp2 (i)))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! bigg mechanism (heterogeneous freezing of cloud water to cloud ice)
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = tice0 - pt1 (i)
            if (ql (i, j) > qcmin .and. dtmp > 0.) then
                newqi = new_ice_condensate(pt1 (i), ql (i,j), qi (i,j))
                sink (i) = 3.3333e-10 * dt_bigg * (exp (0.66 * dtmp) - 1.) * den (i) * ql (i, j) ** 2
                sink (i) = max(0.0,min (newqi, fac_frz * dtmp / icp2 (i), sink (i)))
                ql (i, j) = ql (i, j) - sink (i)
                qi (i, j) = qi (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! freezing of rain to graupel
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = tice0 - pt1 (i)
            if (qr (i, j) > qpmin .and. dtmp > 0.) then
                tmp = min (1., (dtmp * 0.025) ** 2) * qr (i, j) ! no limit on freezing below - 40 deg c
                sink (i) = min (tmp, fac_r2g * dtmp / icp2 (i))
                qr (i, j) = qr (i, j) - sink (i)
                qg (i, j) = qg (i, j) + sink (i)
                q_liq (i) = q_liq (i) - sink (i)
                q_sol (i) = q_sol (i) + sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) + sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! melting of snow to rain or cloud water
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            dtmp = pt1 (i) - tice0
            if (qs (i, j) > qpmin .and. dtmp > 0.) then
                tmp = min (1., (dtmp * 0.1) ** 2) * qs (i, j) ! no limter on melting above 10 deg c
                sink (i) = min (tmp, fac_smlt * dtmp / icp2 (i))
                tmp = min (sink (i), dim (qs_mlt, ql (i, j))) ! max ql due to snow melt
                qs (i, j) = qs (i, j) - sink (i)
                ql (i, j) = ql (i, j) + tmp
                qr (i, j) = qr (i, j) + sink (i) - tmp
                ! qr (i, j) = qr (i, j) + sink (i)
                q_liq (i) = q_liq (i) + sink (i)
                q_sol (i) = q_sol (i) - sink (i)
                cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
                pt1 (i) = pt1 (i) - sink (i) * lhi (i) / cvm (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud water to rain
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (ql (i, j) > ql0_max) then
                sink (i) = fac_l2r * (ql (i, j) - ql0_max)
                qr (i, j) = qr (i, j) + sink (i)
                ql (i, j) = ql (i, j) - sink (i)
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            lhl (i) = lv00 + d0_vap * pt1 (i)
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
            tcp2 (i) = lcp2 (i) + icp2 (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! sublimation / deposition between water vapor and cloud ice
        ! -----------------------------------------------------------------------
        if (do_qa) then
          do i = is, ie
            src (i) = 0.
            if (pt1 (i) < t_sub) then ! too cold to be accurate; freeze qv as a fix
                src (i) = dim (qv (i, j), qvmin)
            elseif (pt1 (i) < tice0) then
                qsi = iqs2 (pt1 (i), den (i), dqsdt)
                dq = (qv (i, j) - qsi)
                sink (i) = fac_i2v * dq / (1. + tcp2 (i) * dqsdt)
                if (qi (i, j) > qcmin) then
                    pidep = sdt * dq * 349138.78 * exp (0.875 * log (qi (i, j) * den (i))) &
                         / (qsi * den (i) * lat2 / (0.0243 * rvgas * pt1 (i) ** 2) + 4.42478e4)
                else
                    pidep = 0.
                endif
                if (dq > 0.) then ! vapor - > ice
                    tmp = tice - pt1 (i)
           ! WRF    qi_crt = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.15 * tmp)))
           ! GFDL   qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den (i)
           ! GEOS impose CALIPSO ice polynomial from 0 C to -40 C on qi_crt  
                    qi_crt = calipso_ice_polynomial(pt1(i)) * qi_gen / den (i)
                    src (i) = min (sink (i), max (qi_crt - qi (i, j), pidep), tmp / tcp2 (i))
                else
                    pidep = pidep * min (1., dim (pt1 (i), t_sub) * 0.2)
                    src (i) = max (pidep, sink (i), - qi (i, j))
                endif
            endif
            qv (i, j) = qv (i, j) - src (i)
            qi (i, j) = qi (i, j) + src (i)
            q_sol (i) = q_sol (i) + src (i)
            cvm (i) = mc_air (i) + qv (i, j) * c_vap + q_liq (i) * c_liq + q_sol (i) * c_ice
            pt1 (i) = pt1 (i) + src (i) * (lhl (i) + lhi (i)) / cvm (i)
          enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! virtual temp updated
        ! -----------------------------------------------------------------------
        
        do i = is, ie
#ifdef USE_COND
            q_cond (i) = q_liq (i) + q_sol (i)
            tmp = 1. + zvir * qv (i, j)
            pt (i, j) = pt1 (i) * tmp * (1. - q_cond (i))
            tmp = rdgas * tmp
            cappa (i, j) = tmp / (tmp + cvm (i))
#else
            pt (i, j) = pt1 (i) * (1. + zvir * qv (i, j))
#endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! fix negative graupel with available cloud ice
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            if (qg (i, j) < 0.) then
                tmp = min (- qg (i, j), max (0., qi (i, j)))
                qg (i, j) = qg (i, j) + tmp
                qi (i, j) = qi (i, j) - tmp
            endif
        enddo
        
        ! -----------------------------------------------------------------------
        ! autoconversion from cloud ice to snow
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            qim = qi0_crt * (1.0 - min(qi_lim,exp(0.05 * (pt1(i)-tice)))) / den (i)
            if (qi (i, j) > qim) then
                sink (i) = fac_i2s * (qi (i, j) - qim)
                qi (i, j) = qi (i, j) - sink (i)
                qs (i, j) = qs (i, j) + sink (i)
            endif
        enddo
        
        if (out_dt) then
            do i = is, ie
                dtdt (i, j) = dtdt (i, j) + pt1 (i) - t0 (i)
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! fix energy conservation
        ! -----------------------------------------------------------------------
        
        if (consv_te) then
            do i = is, ie
                if (hydrostatic) then
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + c_air * pt1 (i))
                else
#ifdef USE_COND
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + cvm (i) * pt1 (i))
#else
                    te0 (i, j) = dp (i, j) * (te0 (i, j) + c_air * pt1 (i))
#endif
                endif
            enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! update latend heat coefficient
        ! -----------------------------------------------------------------------
        
        do i = is, ie
            lhi (i) = li00 + dc_ice * pt1 (i)
            lhl (i) = lv00 + d0_vap * pt1 (i)
            cvm (i) = mc_air (i) + (qv (i, j) + q_liq (i) + q_sol (i)) * c_vap
            lcp2 (i) = lhl (i) / cvm (i)
            icp2 (i) = lhi (i) / cvm (i)
        enddo
        
        ! -----------------------------------------------------------------------
        ! compute cloud fraction
        ! -----------------------------------------------------------------------
        
        if (do_qa .and. last_step) then
            ! FV3 will do the cloud PDF

            ! -----------------------------------------------------------------------
            ! combine water species
            ! -----------------------------------------------------------------------
            
            if (preciprad) then
                do i = is, ie
                   q_sol (i) = qi (i, j) + qs (i, j) + qg (i, j)
                   q_liq (i) = ql (i, j) + qr (i, j)
                enddo
            else
                do i = is, ie
                   q_sol (i) = qi (i, j)
                   q_liq (i) = ql (i, j)
                enddo
            endif
            do i = is, ie
                q_cond (i) = q_sol (i) + q_liq (i)
            enddo
            
            ! -----------------------------------------------------------------------
            ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
            ! -----------------------------------------------------------------------
            
            do i = is, ie
                
                tin = pt1 (i) - (lcp2 (i) * q_cond (i) + icp2 (i) * q_sol (i)) ! minimum temperature
                ! tin = pt1 (i) - ((lv00 + d0_vap * pt1 (i)) * q_cond (i) + &
                ! (li00 + dc_ice * pt1 (i)) * q_sol (i)) / (mc_air (i) + qpz (i) * c_vap)
                
                ! -----------------------------------------------------------------------
                ! determine saturated specific humidity
                ! -----------------------------------------------------------------------
                
                if (tin <= t_wfr) then
                    ! ice phase:
                    qstar (i) = iqs1 (tin, den (i))
                elseif (tin >= tice) then
                    ! liquid phase:
                    qstar (i) = wqs1 (tin, den (i))
                else
                    ! mixed phase:
                    qsi = iqs1 (tin, den (i))
                    qsw = wqs1 (tin, den (i))
                    rqi = calipso_ice_polynomial(tin)
                    qstar (i) = rqi * qsi + (1. - rqi) * qsw
                endif
                
                ! -----------------------------------------------------------------------
                ! partial cloudiness by pdf:
                ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
                ! binary cloud scheme; qa = 0.5 if qstar (i) == qpz
                ! -----------------------------------------------------------------------
                
                rh = qpz (i) / qstar (i)
                
                ! -----------------------------------------------------------------------
                ! icloud_f = 0: bug - fixed
                ! icloud_f = 1: old fvgfs gfdl) mp implementation
                ! icloud_f = 2: top-hat pdf !!! used to be binary cloud scheme (0 / 1)
                ! icloud_f = 3: triangular pdf
                ! -----------------------------------------------------------------------
                if (rh > 0.75 .and. qpz (i) > qcmin) then
                    dq = hvar (i) * qpz (i)
                    q_plus = qpz (i) + dq
                    q_minus = qpz (i) - dq
                    if (icloud_f == 3) then
                     ! triangular
                     if(q_plus.le.qstar(i)) then
                       qa (i, j) = qcmin  ! little/no cloud cover 
                     elseif ( (qpz(i).le.qstar(i)).and.(qstar(i).lt.q_plus) ) then ! partial cloud cover
                       qa (i, j) = max(qcmin,min(1., qa (i, j) + (q_plus-qstar(i))*(q_plus-qstar(i)) / ( (q_plus-q_minus)*(q_plus-qpz(i)) )))
                     elseif ( (q_minus.le.qstar(i)).and.(qstar(i).lt.qpz(i)) ) then ! partial cloud cover
                       qa (i, j) = max(qcmin,min(1., 1. - ( (qstar(i)-q_minus)*(qstar(i)-q_minus) / ( (q_plus-q_minus)*(qpz(i)-q_minus) ))))
                      elseif ( qstar(i).le.q_minus ) then
                       qa (i, j) = 1.0  ! air fully saturated; 100 % cloud cover
                      endif
                    else
                      if (icloud_f == 2) then
                       ! top-hat
                       if(q_plus.le.qstar(i)) then
                         qa (i, j) = qcmin  ! little/no cloud cover 
                       elseif (qstar(i) < q_plus .and. q_cond (i) > qcmin) then
                         qa (i, j) = max(qcmin, min(1., (q_plus - qstar(i)) / (dq + dq) )) ! partial cloud cover
                       elseif (qstar(i) .le. q_minus) then
                         qa (i, j) = 1.0  ! air fully saturated; 100 % cloud cover
                       endif
                      else
                        if (qstar (i) < q_minus) then
                            qa (i, j) = 1.
                        else
                            if (qstar (i) < q_plus) then
                                if (icloud_f == 0) then
                                    qa (i, j) = (q_plus - qstar (i)) / (dq + dq)
                                else
                                    qa (i, j) = (q_plus - qstar (i)) / (2. * dq * (1. - q_cond (i)))
                                endif
                            else
                                qa (i, j) = 0.
                            endif
                            ! impose minimum cloudiness if substantial q_cond (i) exist
                            if (q_cond (i) > qcmin) then
                                qa (i, j) = max (cld_min, qa (i, j))
                            endif
                            qa (i, j) = min (1., qa (i, j))
                        endif
                      endif
                    endif
                else
                    qa (i, j) = 0.
                endif
                
            enddo
            
        endif
        
    enddo ! end j loop
    
end subroutine fv_sat_adj

! =======================================================================
!>@brief the function 'wqs1' computes the 
!! saturated specific humidity for table ii
! =======================================================================
real function wqs1 (ta, den)
    
    implicit none
    
    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density
    
    real, intent (in) :: ta, den
    
    real :: es, ap1, tmin
    
    integer :: it
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs1 = es / (rvgas * ta * den)
    
end function wqs1

! =======================================================================
!>@brief the function 'wqs1' computes the  saturated specific humidity 
!! for table iii
! =======================================================================
real function iqs1 (ta, den)
    
    implicit none
    
    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density
    
    real, intent (in) :: ta, den
    
    real :: es, ap1, tmin
    
    integer :: it
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs1 = es / (rvgas * ta * den)
    
end function iqs1

! =======================================================================
!>@brief The function 'wqs2'computes the gradient of saturated specific 
!! humidity for table ii
! =======================================================================
real function wqs2 (ta, den, dqdt)
    
    implicit none
    
    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density
    
    real, intent (in) :: ta, den
    
    real, intent (out) :: dqdt
    
    real :: es, ap1, tmin
    
    integer :: it
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)
    
end function wqs2

! =======================================================================
!>@brief The function wqs2_vect computes the gradient of saturated 
!! specific humidity for table ii.
!! It is the same as "wqs2", but written as vector function.
! =======================================================================
subroutine wqs2_vect (is, ie, ta, den, wqsat, dqdt)
    
    implicit none
    
    ! pure water phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density
    
    integer, intent (in) :: is, ie
    
    real, intent (in), dimension (is:ie) :: ta, den
    
    real, intent (out), dimension (is:ie) :: wqsat, dqdt
    
    real :: es, ap1, tmin
    
    integer :: i, it
    
    tmin = tice - 160.
    
    do i = is, ie
        ap1 = 10. * dim (ta (i), tmin) + 1.
        ap1 = min (2621., ap1)
        it = ap1
        es = tablew (it) + (ap1 - it) * desw (it)
        wqsat (i) = es / (rvgas * ta (i) * den (i))
        it = ap1 - 0.5
        ! finite diff, del_t = 0.1:
        dqdt (i) = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta (i) * den (i))
    enddo
    
end subroutine wqs2_vect

! =======================================================================
!>@brief The function 'iqs2' computes the gradient of saturated specific 
!! humidity for table iii.
! =======================================================================
real function iqs2 (ta, den, dqdt)
    
    implicit none
    
    ! water - ice phase; universal dry / moist formular using air density
    ! input "den" can be either dry or moist air density
    
    real, intent (in) :: ta, den
    
    real, intent (out) :: dqdt
    
    real :: es, ap1, tmin
    
    integer :: it
    
    tmin = tice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)
    
end function iqs2

! =======================================================================
! initialization
! prepare saturation water vapor pressure tables
! =======================================================================
!>@brief The subroutine 'qs_init' initializes lookup tables for the saturation mixing ratio.
subroutine qs_init (kmp)
    
    implicit none
    
    integer, intent (in) :: kmp
    
    integer, parameter :: length = 2621
    
    integer :: i
    
    if (mp_initialized) return
    
    if (is_master ()) write (*, *) 'top layer for gfdl_mp = ', kmp
    
    ! generate es table (dt = 0.1 deg c)
    
    allocate (table (length))
    allocate (table2 (length))
    allocate (tablew (length))
    allocate (des2 (length))
    allocate (desw (length))
    
    call qs_table (length)
    call qs_table2 (length)
    call qs_tablew (length)
    
    do i = 1, length - 1
        des2 (i) = max (0., table2 (i + 1) - table2 (i))
        desw (i) = max (0., tablew (i + 1) - tablew (i))
    enddo
    des2 (length) = des2 (length - 1)
    desw (length) = desw (length - 1)
    
    mp_initialized = .true.
    
end subroutine qs_init

! =======================================================================
! saturation water vapor pressure table i
! 3 - phase table
! =======================================================================

subroutine qs_table (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, esh20
    real (kind = r_grid) :: wice, wh2o, fac0, fac1, fac2
    real (kind = r_grid) :: esupc (400)
    
    integer :: i
    real :: tk
 
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------
    
    do i = 1, 1600
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / tice) + fac1) / rvgas
        table (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! compute es over water between - 40 deg c and 102 deg c.
    ! -----------------------------------------------------------------------
    
    do i = 1, 1421
        tem = 233.16 + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        esh20 = e00 * exp (fac2)
        if (i <= 400) then
            esupc (i) = esh20
        else
            table (i + 1200) = esh20
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 40 deg c and 0 deg c
    ! -----------------------------------------------------------------------
    
    do i = 1, 400
        tem = 233.16 + delt * real (i - 1)
     !  wice = 0.05 * (tice - tem)
     !  wh2o = 0.05 * (tem - 253.16)
        tk = tem
        wice = calipso_ice_polynomial(tk)
        wh2o = 1.0 - wice
        table (i + 1200) = wice * table (i + 1200) + wh2o * esupc (i)
    enddo
    
end subroutine qs_table

! =======================================================================
! saturation water vapor pressure table ii
! 1 - phase table
! =======================================================================

subroutine qs_tablew (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem, fac0, fac1, fac2
    
    integer :: i
    
    tmin = tice - 160.
    
    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------
    
    do i = 1, n
        tem = tmin + delt * real (i - 1)
        fac0 = (tem - tice) / (tem * tice)
        fac1 = fac0 * lv0
        fac2 = (dc_vap * log (tem / tice) + fac1) / rvgas
        tablew (i) = e00 * exp (fac2)
    enddo
    
end subroutine qs_tablew

! =======================================================================
! saturation water vapor pressure table iii
! 2 - phase table
! =======================================================================

subroutine qs_table2 (n)
    
    implicit none
    
    integer, intent (in) :: n
    
    real (kind = r_grid) :: delt = 0.1
    real (kind = r_grid) :: tmin, tem0, tem1, fac0, fac1, fac2
    
    integer :: i, i0, i1
    
    tmin = tice - 160.
    
    do i = 1, n
        tem0 = tmin + delt * real (i - 1)
        fac0 = (tem0 - tice) / (tem0 * tice)
        if (i <= 1600) then
            ! -----------------------------------------------------------------------
            ! compute es over ice between - 160 deg c and 0 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * li2
            fac2 = (d2ice * log (tem0 / tice) + fac1) / rvgas
        else
            ! -----------------------------------------------------------------------
            ! compute es over water between 0 deg c and 102 deg c.
            ! -----------------------------------------------------------------------
            fac1 = fac0 * lv0
            fac2 = (dc_vap * log (tem0 / tice) + fac1) / rvgas
        endif
        table2 (i) = e00 * exp (fac2)
    enddo
    
    ! -----------------------------------------------------------------------
    ! smoother around 0 deg c
    ! -----------------------------------------------------------------------
    
    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1
    
end subroutine qs_table2

real function new_ice_condensate(tk, qlk, qik)

     real, intent(in) :: tk, qlk, qik
     real :: ptc, ifrac

     ifrac = calipso_ice_polynomial(tk)
     new_ice_condensate = max(0.0,ifrac*(qlk+qik) - qik)

end function new_ice_condensate

real function new_liq_condensate(tk, qlk, qik)

     real, intent(in) :: tk, qlk, qik
     real :: lfrac

     lfrac = 1.0 - calipso_ice_polynomial(tk)
     if (qlk+qik > qcmin) then
        new_liq_condensate = max(0.0,lfrac*(qlk+qik) - qlk)
     else
        new_liq_condensate = 0.0
     endif

end function new_liq_condensate

real function calipso_ice_polynomial(tk)
  ! Citation: Hu, Y., S. Rodier, K. Xu, W. Sun, J. Huang, B. Lin, P. Zhai, and D. Josset (2010), 
  !           Occurrence, liquid water content, and fraction of supercooled water clouds from 
  !           combined CALIOP/IIR/MODIS measurements, J. Geophys. Res., 115, D00H34, 
  !           doi:10.1029/2009JD012384. 
     real, intent(in) :: tk ! temperature in K
     real :: tc, ptc

     tc = min(0.0,max(t_wfr-tice, tk-tice)) ! convert to celcius
     ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
     calipso_ice_polynomial = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
   ! Returning the fraction of ice for given T(K)

end function calipso_ice_polynomial

end module fv_cmp_mod
