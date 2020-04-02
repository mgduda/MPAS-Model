! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

module mpas_atm_dissipation_models

  use mpas_kind_types, only : RKIND
  use mpas_atmphys_constants
  use mpas_constants
  use mpas_log
  use mpas_derived_types, only : MPAS_LOG_CRIT

  logical, parameter :: debug_dissipation = .false.
  logical, parameter :: les_test = .true.
  !! real (kind=RKIND), parameter :: tke_heat_flux = 0.03  ! shear case from Moeng et al., first hour
  real (kind=RKIND), parameter :: tke_heat_flux = 0.0
  !! real (kind=RKIND), parameter :: tke_drag_coefficient = 0.0013 ! ocean roughness length
  real (kind=RKIND), parameter :: tke_drag_coefficient = 0.00935
  real (kind=RKIND), parameter :: epsilon_bv = 1.e-06
  real (kind=RKIND), parameter :: c_k = 0.1


    contains

      subroutine smagorinsky_2d( kdiff, u, v, c_s, config_len_disp,                                 &
                                 deformation_coef_c2, deformation_coef_s2, deformation_coef_cs,     &
                                 invDt, h_mom_eddy_visc4, config_visc4_2dsmag, h_theta_eddy_visc4,  &
                                 cellStart, cellEnd, nEdgesOnCell, edgesOnCell,                     &
                                 nCells, nEdges                                                    )

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here

      implicit none
      
      integer, intent(in) :: cellStart, cellEnd, nCells, nEdges
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: v
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_c2
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_s2
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_cs
      real (kind=RKIND), intent(in) :: c_s, config_len_disp, invDt, config_visc4_2dsmag
      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(out) :: kdiff
      real (kind=RKIND), intent(out) :: h_mom_eddy_visc4, h_theta_eddy_visc4

      ! local variables

      integer :: iCell, iEdge, k
      real (kind=RKIND), dimension(nVertLevels) :: d_11, d_22, d_12, dudx, dudy, dvdx, dvdy

        if(debug_dissipation)  call mpas_log_write(' begin smagorinsky_2d ')


            do iCell = cellStart,cellEnd
               dudx(1:nVertLevels) = 0.0
               dudy(1:nVertLevels) = 0.0
               dvdx(1:nVertLevels) = 0.0
               dvdy(1:nVertLevels) = 0.0
               do iEdge=1,nEdgesOnCell(iCell)
                  do k=1,nVertLevels
                     dudx(k) = dudx(k) + deformation_coef_c2(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell)) &
                                       - deformation_coef_cs(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                     dudy(k) = dudy(k) + deformation_coef_cs(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell)) &
                                       - deformation_coef_s2(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                     dvdx(k) = dvdx(k) + deformation_coef_cs(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell)) &
                                       + deformation_coef_c2(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                     dvdy(k) = dvdy(k) + deformation_coef_s2(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell)) &
                                       + deformation_coef_cs(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                  end do
               end do
!DIR$ IVDEP
               do k=1, nVertLevels
                  ! here is the Smagorinsky formulation, 
                  ! followed by imposition of an upper bound on the eddy viscosity
                  d_11(k) = 2*dudx(k)
                  d_22(k) = 2*dvdy(k)
                  d_12(k) = dudy(k) + dvdx(k)
                  kdiff(k,iCell) = (c_s * config_len_disp)**2 * sqrt(0.25*(d_11(k)-d_22(k))**2 + d_12(k)**2)
                  kdiff(k,iCell) = min(kdiff(k,iCell),(0.01*config_len_disp**2) * invDt)
               end do
            end do

            h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
            h_theta_eddy_visc4 = h_mom_eddy_visc4

        if(debug_dissipation)  call mpas_log_write(' exiting smagorinsky_2d ')

      end subroutine smagorinsky_2d

!---------------------------------------

      subroutine les_models( config_les_model, eddy_visc_horz, eddy_visc_vert,                  &
                             u, v, uCell, vCell,                                                &
                             w, c_s, bv_freq2, zgrid, config_len_disp,                          &
                             deformation_coef_c2, deformation_coef_s2, deformation_coef_cs,     &
                             deformation_coef_c, deformation_coef_s,                            &
                             invDt, h_mom_eddy_visc4, config_visc4_2dsmag, h_theta_eddy_visc4,  &
                             scalars, tend_scalars, index_tke, rho_zz, meshScalingDel2,         &
                             cellStart, cellEnd, nEdgesOnCell, edgesOnCell, cellsOnEdge,        &
                             nCells, nEdges, nVertLevels, maxEdges, num_scalars                )

      implicit none
      
      character (len=StrKIND), intent(in) :: config_les_model

      integer, intent(in) :: cellStart, cellEnd, nCells, nEdges, nVertLevels, maxEdges, index_tke, num_scalars
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: v
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: uCell
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: vCell
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: bv_freq2
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: scalars
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1) :: tend_scalars
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: w
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: zgrid
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_c2
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_s2
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_cs
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_c
      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: deformation_coef_s
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), intent(in) :: c_s, config_len_disp, invDt, config_visc4_2dsmag
      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell
      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(out) :: eddy_visc_horz, eddy_visc_vert
      real (kind=RKIND), intent(out) :: h_mom_eddy_visc4, h_theta_eddy_visc4

      ! local variables

      integer :: iCell, iEdge, k, ie, cell1, cell2
      real (kind=RKIND), dimension(nVertLevels) :: d_11, d_22, d_33, d_12, d_13, d_23 
      real (kind=RKIND), dimension(nVertLevels) :: dudx, dudy, dvdx, dvdy
      real (kind=RKIND), dimension(nVertLevels+1) :: dwdx, dwdy
      real (kind=RKIND), dimension(nVertLevels) :: dudz, dvdz, dwdz
      real (kind=RKIND) :: rdz, def2, pr_inv, wk
      real (kind=RKIND) :: shear_production, buoyancy, dissipation, delta_z, delta_s, bv, tke_length, diss_length
      real (kind=RKIND) :: l_horizontal, l_vertical, c_dissipation
      real (kind=RKIND) :: prandtl_horizontal_inv, prandtl_vertical_inv
      real (kind=RKIND) :: eddy_visc_h, eddy_visc_v

      logical, parameter :: test_tke=.true.
      ! real (kind=RKIND), parameter :: epsilon_bv = 1.e-06

      if(debug_dissipation)  call mpas_log_write(' begin les_models ')
      if(debug_dissipation)  call mpas_log_write(' les scheme is '//trim(config_les_model))

      pr_inv = 1./prandtl

      ! set up coefficients for 4th-order horizontal background filter

      h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
      h_theta_eddy_visc4 = h_mom_eddy_visc4

      do iCell = cellStart,cellEnd

        dudx(1:nVertLevels) = 0.0
        dudy(1:nVertLevels) = 0.0
        dvdx(1:nVertLevels) = 0.0
        dvdy(1:nVertLevels) = 0.0
        dwdx(1:nVertLevels+1) = 0.0
        dwdy(1:nVertLevels+1) = 0.0

        dudz(1:nVertLevels) = 0.0
        dvdz(1:nVertLevels) = 0.0
        dwdz(1:nVertLevels) = 0.0
               
        do iEdge=1,nEdgesOnCell(iCell)

          ie = EdgesOnCell(iEdge,iCell)
          cell1 = cellsOnEdge(1,ie)
          cell2 = cellsOnEdge(2,ie)

          do k=1,nVertLevels
            dudx(k) = dudx(k) + deformation_coef_c2(iEdge,iCell)*u(k,ie) &
                              - deformation_coef_cs(iEdge,iCell)*v(k,ie)
            dudy(k) = dudy(k) + deformation_coef_cs(iEdge,iCell)*u(k,ie) &
                              - deformation_coef_s2(iEdge,iCell)*v(k,ie)
            dvdx(k) = dvdx(k) + deformation_coef_cs(iEdge,iCell)*u(k,ie) &
                              + deformation_coef_c2(iEdge,iCell)*v(k,ie)
            dvdy(k) = dvdy(k) + deformation_coef_s2(iEdge,iCell)*u(k,ie) &
                              + deformation_coef_cs(iEdge,iCell)*v(k,ie)
          end do

          do k=1,nVertLevels+1
            wk = 0.5*(w(k,cell1)+w(k,cell2))
             dwdx(k) = dwdx(k) + deformation_coef_c(iEdge,iCell)*wk
             dwdy(k) = dwdy(k) + deformation_coef_s(iEdge,iCell)*wk
          end do

        end do

        do k=1,nVertLevels
          rdz = 1./(zgrid(k+1,iCell)-zgrid(k,iCell))
          dwdz(k) = (w(k+1,iCell)-w(k,iCell))*rdz
        end do

        do k=2,nVertLevels-1
          rdz = 1./(zgrid(k+2,iCell)+zgrid(k+1,iCell)-zgrid(k,iCell)-zgrid(k-1,iCell))
          dudz(k) = (uCell(k+1,iCell)-uCell(k-1,iCell))*rdz
          dvdz(k) = (vCell(k+1,iCell)-vCell(k-1,iCell))*rdz
        end do

        k = 1
        rdz = 1./(zgrid(k+1,iCell)-zgrid(k,iCell))
        dudz(k) = (uCell(k+1,iCell)-uCell(k,iCell))*rdz
        dvdz(k) = (vCell(k+1,iCell)-vCell(k,iCell))*rdz

        k = nVertLevels-1
        rdz = 1./(zgrid(k+1,iCell)-zgrid(k,iCell))
        dudz(k+1) = (uCell(k+1,iCell)-uCell(k,iCell))*rdz
        dvdz(k+1) = (vCell(k+1,iCell)-vCell(k,iCell))*rdz

        do k=1, nVertLevels
          d_11(k) = 2.*dudx(k)
          d_22(k) = 2.*dvdy(k)
          d_33(k) = 2.*dwdz(k)
          d_12(k) = dudy(k) + dvdx(k)
          d_13(k) = dwdx(k) + dudz(k)
          d_23(k) = dwdy(k) + dvdz(k)
        end do

        if (config_les_model == "3d_smagorinsky") then

          do k=1, nVertLevels
            def2 = 0.5*(d_11(k)**2 + d_22(k)**2 + d_33(k)**2) + d_12(k)**2 + d_13(k)**2 + d_23(k)**2
            eddy_visc_horz(k,iCell) = (c_s * config_len_disp)**2 * sqrt(max(0.,def2 - pr_inv*bv_freq2(k,iCell)))
            eddy_visc_horz(k,iCell) = min(eddy_visc_horz(k,iCell),(0.01*config_len_disp**2) * invDt)
            eddy_visc_vert(k,iCell) = eddy_visc_horz(k,iCell)
          end do

        else if (config_les_model == "prognostic_1.5_order") then

          do k=1,nVertLevels  ! bound the tke here, currently hardwired
            scalars(index_tke,k,iCell) = max(0.,min(100.,scalars(index_tke,k,iCell)))
          end do

          do k=1,nVertLevels

            delta_z = zgrid(k+1,iCell)-zgrid(k,iCell)
            delta_s = ((config_len_disp**2)*delta_z)**(1./3.)
            bv = max( sqrt(abs(bv_freq2(k,iCell))), epsilon_bv )
            tke_length = delta_s
            if(bv_freq2(k,iCell) .gt. 1.e-06) &
                tke_length = 0.76*sqrt(scalars(index_tke,k,iCell))/bv
            diss_length = min(delta_s,max(tke_length,0.01*delta_s))
            if(bv_freq2(k,iCell) <= 0) diss_length = delta_s

            ! non-isotropic mixing                                                                            

            l_horizontal = config_len_disp
            l_vertical = min(delta_z,tke_length)
            if(bv_freq2(k,iCell) <= 0) diss_length = delta_z

            ! isotropic mixing

            ! l_horizontal = min(delta_s,tke_length)
            ! if(bv_freq2(k,iCell) <= 0) diss_length = delta_s
            ! l_vertical = l_horizontal

            ! eddy viscocities set here if we are running the 1.5 order prognostic tke scheme
            eddy_visc_h = c_k*l_horizontal*sqrt(scalars(index_tke,k,iCell))
            eddy_visc_v = c_k*l_vertical*sqrt(scalars(index_tke,k,iCell))

            ! testing
            eddy_visc_horz(k,iCell) = eddy_visc_h
            eddy_visc_vert(k,iCell) = eddy_visc_v

            ! terms for the prognostic tke integration

            shear_production = eddy_visc_h*(d_11(k)**2 + d_22(k)**2 + d_12(k)**2) &
                              +eddy_visc_v*(d_33(k)**2 + d_13(k)**2 + d_13(k)**2)

            buoyancy = -eddy_visc_v*bv_freq2(k,iCell)

            ! dissipation                                                                                     

            c_dissipation = 1.9*c_k + max( 0.0, 0.93 - 1.9*c_k )*diss_length/delta_s
            !  if( (k.eq. 1) .or. (k.eq.nVertLevels) ) c_dissipation = 3.9

            dissipation = -c_dissipation*(scalars(index_tke,k,iCell)**(1.5))/diss_length

            ! computing eddy viscosities *********                                                            

            prandtl_horizontal_inv = 3.
            prandtl_vertical_inv = 1.0+(2.0*l_vertical/delta_z) ! going to need to store this off somewhere   

            ! RHS term for the subgrid ke.                                                                    

            tend_scalars(index_tke,k,iCell) = rho_zz(k,iCell)*( shear_production + buoyancy + dissipation )

          end do

        else

          call mpas_log_write(' in les_models, no les scheme for '//trim(config_les_model), messageType=MPAS_LOG_CRIT)
          
        end if ! end of config_les_model test

      end do ! loop over all owned cells (columns)

      if(debug_dissipation)  call mpas_log_write(' les_models ')

      end subroutine les_models

!---------------------------------------

      subroutine calculate_n2( bn2, theta_m, exner, pressure_b, pp, zgrid, scalars, index_qv, index_qc, qtot, &
                               cellStart, cellEnd, nCells)

      use mpas_atm_dimensions  !  pull nVertLevels and num_scalars from here

      integer, intent(in) :: cellStart, cellEnd, nCells
      integer, intent(in) :: index_qv, index_qc
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(out) :: bn2
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m, exner, pressure_b, pp, qtot
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: zgrid
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(in) :: scalars
! local
      real (kind=RKIND) :: dz, rdz, esw
      real (kind=RKIND), parameter :: qc_cr = 0.00001  ! in kg/kg
      real (kind=RKIND), dimension(nVertLevels) :: theta, qvsw, temp, coefa
      logical :: dry_bv_frequency

         if(debug_dissipation)  call mpas_log_write(' begin BV frequency calculations ')

            do iCell = cellStart,cellEnd
!DIR$ IVDEP
               do k=1, nVertLevels

                  theta(k) = theta_m(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))

                  temp(k) = exner(k,iCell) * theta(k)

                  p = pressure_b(k,iCell) + pp(k,iCell)
                  esw = 1000. * svp1 * exp(svp2 * (temp(k) - svpt0) / (temp(k) - svp3))
                  if (p < esw) esw = p * 0.99     ! fix for pressure < esw
                  qvsw(k) = ep_2 * esw / (p - esw)

                  coefa(k) = ( 1.0 + xlv * qvsw(k)/ R_d / temp(k) ) / &
                     ( 1.0 + xlv * xlv *qvsw(k) / Cp / R_v / temp(k) / temp(k) ) 

               end do

               do k=2, nVertLevels-1
                  dz = 0.5 * (zgrid(k+2,iCell)+zgrid(k+1,iCell)) - 0.5 * (zgrid(k,iCell)+zgrid(k-1,iCell))
                  rdz = 1.0/dz

                  ! if ( scalars(index_qc,k,iCell) < qc_cr ) then
                  !   ! Dry Brunt-Vaisala frequency
                  !   bn2(k,iCell) = gravity * ((theta(k+1) - theta(k-1) ) / theta(k)  / dz  &
                  !                + rvord * (scalars(index_qv,k+1,iCell) - scalars(index_qv,k-1,iCell)) / dz &
                  !               - ( qtot(k+1, iCell) - qtot(k-1, iCell) ) / dz )
                  ! else
                  !   ! Moist Brunt-Vaisala frequency according to Durran and Klemp (1982) Eq. 36
                  !   bn2(k,iCell) = gravity * ( coefa(k) * ((theta(k+1) - theta(k-1) ) / theta(k) / dz & 
                  !               + xlv / cp / temp(k) * ( qvsw(k+1) - qvsw(k-1)) / dz ) &
                  !               - ( qtot(k+1, iCell) - qtot(k-1, iCell) ) / dz )
                  ! endif

                  dry_bv_frequency = .true.
                  if(index_qc .gt. 0) then  ! if moist simulation, qc exists
                    if ( scalars(index_qc,k,iCell) .ge. qc_cr ) dry_bv_frequency = .false.
                  end if

                  if (dry_bv_frequency) then
                     ! Dry Brunt-Vaisala frequency
                     bn2(k,iCell) = gravity * ((theta(k+1) - theta(k-1) ) / theta(k)  * rdz  &
                                  + rvord * (scalars(index_qv,k+1,iCell) - scalars(index_qv,k-1,iCell)) * rdz &
                                 - ( qtot(k+1, iCell) - qtot(k-1, iCell) ) * rdz )
                  else
                     ! Moist Brunt-Vaisala frequency according to Durran and Klemp (1982) Eq. 36
                     bn2(k,iCell) = gravity * ( coefa(k) * ((theta(k+1) - theta(k-1) ) / theta(k) * rdz & 
                                 + xlv / cp / temp(k) * ( qvsw(k+1) - qvsw(k-1)) * rdz ) &
                                 - ( qtot(k+1, iCell) - qtot(k-1, iCell) ) * rdz )
                  endif

               end do
               bn2(1,iCell) = bn2(2,iCell)
               bn2(nVertLevels,iCell) = bn2(nVertLevels-1,iCell)

            end do

         if(debug_dissipation)  call mpas_log_write(' exiting BV frequency calculations ')

      end subroutine calculate_n2

!---------------------------------------

      subroutine u_dissipation( edgeStart, edgeEnd, edgeSolveStart, edgeSolveEnd, vertexStart, vertexEnd,    &
                                cellStart, cellEnd, nCells, nEdges, nVertices, vertexDegree,                 &
                                cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnVertex,                     &
                                nEdgesOnCell, edgesOnCell_sign, edgesOnVertex_sign,                          &
                                invAreaCell, invAreaTriangle, invDvEdge, invDcEdge,                          &
                                angleEdge, dcEdge, dvEdge, meshScalingDel2, meshScalingDel4,                 &
                                config_mix_full, h_mom_eddy_visc4, v_mom_eddy_visc2,                         &
                                config_del4u_div_factor, zgrid, kdiff,                                       &
                                delsq_u, delsq_vorticity, delsq_divergence,                                  &
                                u, divergence, vorticity, rho_edge, u_init, v_init, tend_u_euler            )

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here

      implicit none

      integer, intent(in) :: edgeStart, edgeEnd, edgeSolveStart, edgeSolveEnd
      integer, intent(in) :: vertexStart, vertexEnd, vertexDegree
      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: nCells, nEdges, nVertices
      logical, intent(in) :: config_mix_full

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge
      integer, dimension(2,nEdges+1), intent(in) :: verticesOnEdge
      integer, dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell
      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(vertexDegree,nVertices+1), intent(in) :: edgesOnVertex

      real (kind=RKIND), intent(in) :: h_mom_eddy_visc4
      real (kind=RKIND), intent(in) :: v_mom_eddy_visc2
      real (kind=RKIND), intent(in) :: config_del4u_div_factor

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(vertexDegree,nVertices+1), intent(in) :: edgesOnVertex_sign
      real (kind=RKIND), dimension(nVertices+1), intent(in) :: invAreaTriangle
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: angleEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: zgrid

      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: u
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: divergence
      real (kind=RKIND), dimension(nVertLevels,nVertices+1), intent(in) :: vorticity
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: kdiff


      !  scratch space from calling routine
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: delsq_u
      real (kind=RKIND), dimension(nVertLevels,nVertices+1) :: delsq_vorticity
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_divergence

      real (kind=RKIND), dimension(nVertLevels), intent(in) :: u_init, v_init

      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(out) :: tend_u_euler

      ! local variables

      integer :: iEdge, cell1, cell2, vertex1, vertex2, iVertex, iCell, i, k
      real (kind=RKIND) :: r_dc, r_dv, u_diffusion, kdiffu, r, edge_sign, u_mix_scale
      real (kind=RKIND) :: z1, z2, z3, z4, zm, z0, zp
      real (kind=RKIND), dimension(nVertLevels) :: u_mix

!$OMP BARRIER 

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         ! First, storage to hold the result from the first del^2 computation.

         delsq_u(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            r_dc = invDcEdge(iEdge)
            r_dv = min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP  
            do k=1,nVertLevels

               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               ! only valid for h_mom_eddy_visc4 == constant
              u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) * r_dc  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) * r_dv

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion

               kdiffu = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))

               ! include 2nd-orer diffusion here

               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) &
                                       + rho_edge(k,iEdge)* kdiffu * u_diffusion * meshScalingDel2(iEdge)

            end do
         end do

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active

!$OMP BARRIER
            do iVertex=vertexStart,vertexEnd
               delsq_vorticity(1:nVertLevels,iVertex) = 0.0
               do i=1,vertexDegree
                  iEdge = edgesOnVertex(i,iVertex)
                  edge_sign = invAreaTriangle(iVertex) * dcEdge(iEdge) * edgesOnVertex_sign(i,iVertex)
                  do k=1,nVertLevels
                     delsq_vorticity(k,iVertex) = delsq_vorticity(k,iVertex) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

            do iCell=cellStart,cellEnd
               delsq_divergence(1:nVertLevels,iCell) = 0.0
               r = invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = r * dvEdge(iEdge) * edgesOnCell_sign(i,iCell)
                  do k=1,nVertLevels
                     delsq_divergence(k,iCell) = delsq_divergence(k,iCell) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

!$OMP BARRIER
            do iEdge=edgeSolveStart,edgeSolveEnd
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               vertex1 = verticesOnEdge(1,iEdge)
               vertex2 = verticesOnEdge(2,iEdge)

               u_mix_scale = meshScalingDel4(iEdge)*h_mom_eddy_visc4
               r_dc = u_mix_scale * config_del4u_div_factor * invDcEdge(iEdge)
               r_dv = u_mix_scale * min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP
               do k=1,nVertLevels

                  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                  !                    only valid for h_mom_eddy_visc4 == constant
                  !
                  ! Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor
                  !    relative to the rotational part.  The stability constraint on the divergence component is much less
                  !    stringent than the rotational part, and this flexibility may be useful.
                  !
                  u_diffusion =  rho_edge(k,iEdge) *  ( ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) * r_dc  &
                                                       -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) * r_dv )
                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - u_diffusion

               end do
            end do

         end if ! 4th order mixing is active
      !
      !  vertical mixing for u - 2nd order filter in physical (z) space              
      !
         if ( v_mom_eddy_visc2 > 0.0 ) then

            if (config_mix_full) then  ! mix full state 

               do iEdge=edgeSolveStart,edgeSolveEnd

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=2,nVertLevels-1

                     z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                     z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                     z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                     z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                        (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      &
                                       -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            else  ! idealized cases where we mix on the perturbation from the initial 1-D state

               do iEdge=edgeSolveStart,edgeSolveEnd

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     u_mix(k) = u(k,iEdge) - u_init(k) * cos( angleEdge(iEdge) ) &
                                           + v_init(k) * sin( angleEdge(iEdge) )
                  end do

                  do k=2,nVertLevels-1

                     z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                     z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                     z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                     z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                        (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      &
                                       -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            end if  ! mix perturbation state

         end if  ! vertical mixing of horizontal momentum

      end subroutine u_dissipation

!------------------------

      subroutine u_dissipation_3d( edgeStart, edgeEnd, edgeSolveStart, edgeSolveEnd, vertexStart, vertexEnd,    &
                                   cellStart, cellEnd, nCells, nEdges, nVertices, vertexDegree,                 &
                                   cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnVertex,                     &
                                   nEdgesOnCell, edgesOnCell_sign, edgesOnVertex_sign,                          &
                                   invAreaCell, invAreaTriangle, invDvEdge, invDcEdge,                          &
                                   angleEdge, dcEdge, dvEdge, meshScalingDel2, meshScalingDel4,                 &
                                   config_mix_full, h_mom_eddy_visc4, v_mom_eddy_visc2,                         &
                                   config_del4u_div_factor, zgrid,                                              &
                                   eddy_visc_horz, eddy_visc_vert, zz, rdzu, rdzw,                              &
                                   fzm, fzp, config_les_model,                                                  &
                                   delsq_u, delsq_vorticity, delsq_divergence,                                  &
                                   u, v, divergence, vorticity, rho_edge, rho_zz, u_init, v_init, tend_u_euler )

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             

      implicit none

      integer, intent(in) :: edgeStart, edgeEnd, edgeSolveStart, edgeSolveEnd
      integer, intent(in) :: vertexStart, vertexEnd, vertexDegree
      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: nCells, nEdges, nVertices
      logical, intent(in) :: config_mix_full

      character (len=StrKIND) :: config_les_model

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge
      integer, dimension(2,nEdges+1), intent(in) :: verticesOnEdge
      integer, dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell
      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(vertexDegree,nVertices+1), intent(in) :: edgesOnVertex

      real (kind=RKIND), intent(in) :: h_mom_eddy_visc4
      real (kind=RKIND), intent(in) :: v_mom_eddy_visc2
      real (kind=RKIND), intent(in) :: config_del4u_div_factor

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(vertexDegree,nVertices+1), intent(in) :: edgesOnVertex_sign
      real (kind=RKIND), dimension(nVertices+1), intent(in) :: invAreaTriangle
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: angleEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: zgrid

      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: u
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: v
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: divergence
      real (kind=RKIND), dimension(nVertLevels,nVertices+1), intent(in) :: vorticity
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_horz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_vert
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzm
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzp


      !  scratch space from calling routine
      real (kind=RKIND), dimension(nVertLevels,nEdges+1) :: delsq_u
      real (kind=RKIND), dimension(nVertLevels,nVertices+1) :: delsq_vorticity
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_divergence

      real (kind=RKIND), dimension(nVertLevels), intent(in) :: u_init, v_init

      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(out) :: tend_u_euler

      ! local variables                                                                                                                               

      integer :: iEdge, cell1, cell2, vertex1, vertex2, iVertex, iCell, i, k
      real (kind=RKIND) :: r_dc, r_dv, u_diffusion, kdiffu, r, edge_sign, u_mix_scale
      real (kind=RKIND) :: z1, z2, z3, z4, zm, z0, zp
      real (kind=RKIND), dimension(nVertLevels) :: u_mix

      real (kind=RKIND), dimension(nVertLevels+1) :: turb_vflux
      real (kind=RKIND) :: rho_k_cell1, rho_k_cell2, rho_k_at_w
      real (kind=RKIND) :: zz_cell1, zz_cell2, zz_at_w

      real (kind=RKIND) :: velocity_magnitude

        if(debug_dissipation) then
            call mpas_log_write(' begin u_dissipation_3d ')
            call mpas_log_write(' les model is '//trim(config_les_model))
            call mpas_log_write(' 4th order hyperviscosity is $r ', realArgs=(/h_mom_eddy_visc4/))
            call mpas_log_write(' 4th order divergence factor is $r ', realArgs=(/config_del4u_div_factor/))
        end if

!$OMP BARRIER                                                                                                                                         

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).                                                                          
         ! First, storage to hold the result from the first del^2 computation.                                                                        

         delsq_u(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            r_dc = invDcEdge(iEdge)
            r_dv = min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP                                                                                                                                           
            do k=1,nVertLevels

               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity                                                         
               !                    only valid for h_mom_eddy_visc4 == constant                                                                       
              u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) * r_dc  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) * r_dv

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion

               kdiffu = 0.5*(eddy_visc_horz(k,cell1)+eddy_visc_horz(k,cell2))

               ! include 2nd-orer diffusion here                                                                                                      
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) &
                                       + rho_edge(k,iEdge)* kdiffu * u_diffusion * meshScalingDel2(iEdge)

            end do
         end do

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active                                                                               

!$OMP BARRIER                                                                                                                                         

            do iVertex=vertexStart,vertexEnd
               delsq_vorticity(1:nVertLevels,iVertex) = 0.0
               do i=1,vertexDegree
                  iEdge = edgesOnVertex(i,iVertex)
                  edge_sign = invAreaTriangle(iVertex) * dcEdge(iEdge) * edgesOnVertex_sign(i,iVertex)
                  do k=1,nVertLevels
                     delsq_vorticity(k,iVertex) = delsq_vorticity(k,iVertex) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

            do iCell=cellStart,cellEnd
               delsq_divergence(1:nVertLevels,iCell) = 0.0
               r = invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = r * dvEdge(iEdge) * edgesOnCell_sign(i,iCell)
                  do k=1,nVertLevels
                     delsq_divergence(k,iCell) = delsq_divergence(k,iCell) + edge_sign * delsq_u(k,iEdge)
                  end do
               end do
            end do

!$OMP BARRIER                                                                                                                                         

            do iEdge=edgeSolveStart,edgeSolveEnd
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               vertex1 = verticesOnEdge(1,iEdge)
               vertex2 = verticesOnEdge(2,iEdge)

               u_mix_scale = meshScalingDel4(iEdge)*h_mom_eddy_visc4
               r_dc = u_mix_scale * config_del4u_div_factor * invDcEdge(iEdge)
               r_dv = u_mix_scale * min(invDvEdge(iEdge), 4*invDcEdge(iEdge))

!DIR$ IVDEP                                                                                                                                           
               do k=1,nVertLevels

                  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity                                                      
                  !                    only valid for h_mom_eddy_visc4 == constant                                                                    
                  !                                                                                                                                   
                  ! Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor                                           
                  !    relative to the rotational part.  The stability constraint on the divergence component is much less                            
                  !    stringent than the rotational part, and this flexibility may be useful.                                                        
                  !                                                                                                                                   
                  u_diffusion =  rho_edge(k,iEdge) *  ( ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) * r_dc  &
                                                       -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) * r_dv )
                  tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - u_diffusion

               end do
            end do

         end if ! 4th order mixing is active                                                                                                          

      !                                                                                                                                               
      !  vertical mixing for u - 2nd order filter in physical (z) space                                                                               
      !                                                                                                                                               
         if ( v_mom_eddy_visc2 > 0.0 ) then

            if (config_mix_full) then  ! mix full state                                                                                               

               do iEdge=edgeSolveStart,edgeSolveEnd

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=2,nVertLevels-1

                     z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                     z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                     z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                     z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                        (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      &
                                       -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            else  ! idealized cases where we mix on the perturbation from the initial 1-D state                                                       

               do iEdge=edgeSolveStart,edgeSolveEnd

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     u_mix(k) = u(k,iEdge) - u_init(k) * cos( angleEdge(iEdge) ) &
                                           + v_init(k) * sin( angleEdge(iEdge) )
                  end do

                  do k=2,nVertLevels-1

                     z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
                     z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
                     z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
                     z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                        (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      &
                                       -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            end if  ! mix perturbation state                                                                                             

         end if  ! vertical mixing of horizontal momentum for les formulation

        if ( config_les_model /= "none") then

               do iEdge=edgeSolveStart,edgeSolveEnd

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)
                  turb_vflux(nVertlevels+1) = 0.  !  no turbulent flux out of the domain
                  turb_vflux(1) = 0.  !  lower bc flux handled where ???

                  do k=2,nVertLevels
                     rho_k_cell1 =  fzm(k)*rho_zz(k  ,cell1)*zz(k  ,cell1)*eddy_visc_vert(k  ,cell1) &
                                   +fzp(k)*rho_zz(k-1,cell1)*zz(k-1,cell1)*eddy_visc_vert(k-1,cell1)
                     rho_k_cell2 =  fzm(k)*rho_zz(k  ,cell2)*zz(k  ,cell2)*eddy_visc_vert(k  ,cell2) &
                                   +fzp(k)*rho_zz(k-1,cell2)*zz(k-1,cell2)*eddy_visc_vert(k-1,cell2)
                     rho_k_at_w = 0.5*(rho_k_cell1+rho_k_cell2)

                     zz_cell1 = fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1)
                     zz_cell2 = fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2)
                     zz_at_w = 0.5*(zz_cell1+zz_cell2)
                     turb_vflux(k) = - rho_k_at_w*zz_at_w*rdzu(k)*(u(k,iEdge)-u(k-1,iEdge))
                  end do

                if( les_test ) then
                  velocity_magnitude = sqrt(u(1,iEdge)**2 + v(1,iEdge)**2)
                  turb_vflux(1) = -rho_edge(1,iEdge)*tke_drag_coefficient*u(1,iEdge)*velocity_magnitude
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                  ! end test conditions
                else
                  ! test conditions for supercell case
                  turb_vflux(1) = turb_vflux(2)
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                  ! end test conditions
                end if

                  do k=1,nVertLevels
                    tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - rdzw(k)*(turb_vflux(k+1)-turb_vflux(k))
                  end do

               end do

         end if 

        if(debug_dissipation)  call mpas_log_write(' exiting u_dissipation_3d ')

      end subroutine u_dissipation_3d

!------------------------

      subroutine w_dissipation( cellStart, cellEnd, cellSolveStart, cellSolveEnd,         &
                                nCells, nEdges,                                           &
                                nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, &
                                invAreaCell, invDcEdge, dvEdge,                           &
                                meshScalingDel2, meshScalingDel4,                         &
                                rdzw, rdzu,                                               &
                                v_mom_eddy_visc2, h_mom_eddy_visc4,                       &
                                delsq_w,                                                  &
                                w, rho_edge, kdiff, rho_zz,                               &
                                tend_w_euler                                             )


      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             

      implicit none

      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: nCells, nEdges

      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: EdgesOnCell

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      real (kind=RKIND), intent(in) :: h_mom_eddy_visc4
      real (kind=RKIND), intent(in) :: v_mom_eddy_visc2

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: w
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: kdiff
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(inout) :: tend_w_euler

      ! storage passed in from calling routine                                                                                                        
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_w

      ! local variables                                                                                                                               

      integer :: cell1, cell2, iEdge, iCell, i, k
      real (kind=RKIND) :: r_areaCell, edge_sign, w_turb_flux

!  !OMP BARRIER  why is this openmp barrier here???                                                                                                   

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (w) ).                                                                          
       !                                                                                                                                            
       ! First, storage to hold the result from the first del^2 computation.                                                                        
       !  we copied code from the theta mixing, hence the theta* names.                                                                             


         delsq_w(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_w_euler(1:nVertLevels+1,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               edge_sign = 0.5 * r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

!DIR$ IVDEP                                                                                                                                           
              do k=2,nVertLevels

                  w_turb_flux =  edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))
                  delsq_w(k,iCell) = delsq_w(k,iCell) + w_turb_flux
                  w_turb_flux = w_turb_flux * meshScalingDel2(iEdge) * 0.25 * &
                                  (kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2))
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + w_turb_flux
               end do
            end do
         end do

!$OMP BARRIER                                                                                                                                         

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active                                                                               

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...                                                 
               r_areaCell = h_mom_eddy_visc4 * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell) * invDcEdge(iEdge)

                  do k=2,nVertLevels
                     tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - edge_sign * (delsq_w(k,cell2) - delsq_w(k,cell1))
                  end do

               end do
            end do

         end if ! 4th order mixing is active                                                                                                          

      if ( v_mom_eddy_visc2 > 0.0 ) then  ! vertical mixing

         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP                                                                                                                                           
            do k=2,nVertLevels
               tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  &
                                           (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              &
                                          -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k)
            end do
         end do

      end if

      end subroutine w_dissipation

!------------------------

      subroutine w_dissipation_3d( cellStart, cellEnd, cellSolveStart, cellSolveEnd,         &
                                   nCells, nEdges,                                           &
                                   nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, &
                                   invAreaCell, invDcEdge, dvEdge,                           &
                                   meshScalingDel2, meshScalingDel4,                         &
                                   rdzw, rdzu,                                               &
                                   v_mom_eddy_visc2, h_mom_eddy_visc4,                       &
                                   delsq_w,                                                  &
                                   w, rho_edge, rho_zz, zz,                                  &
                                   eddy_visc_horz, eddy_visc_vert,                           &
                                   config_les_model,                                         &
                                   tend_w_euler                                             )


      ! 3D w dissipation using the 3D smagorinsky eddy viscosities.  
      ! This routine also includes the simpler mixing models, and the 4th-order horizontal filter

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             

      implicit none

      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: nCells, nEdges

      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: EdgesOnCell

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      character (len=StrKIND) :: config_les_model

      real (kind=RKIND), intent(in) :: h_mom_eddy_visc4
      real (kind=RKIND), intent(in) :: v_mom_eddy_visc2

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(in) :: w
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_horz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_vert
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge

      real (kind=RKIND), dimension(nVertLevels+1,nCells+1), intent(inout) :: tend_w_euler

      ! storage passed in from calling routine                                                                                                        
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_w
      real (kind=RKIND), dimension(nVertLevels+1) :: turb_vflux

      ! local variables                                                                                                                               

      integer :: cell1, cell2, iEdge, iCell, i, k
      real (kind=RKIND) :: r_areaCell, edge_sign, w_turb_flux

!  !OMP BARRIER  why is this openmp barrier here???                                                                                                   

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (w) ).                                                                          
       !                                                                                                                                            
       ! First, storage to hold the result from the first del^2 computation.                                                                        
       !  we copied code from the theta mixing, hence the theta* names.                                                                             

        if(debug_dissipation) then
            call mpas_log_write(' begin w_dissipation_3d ')
            call mpas_log_write(' les model is '//trim(config_les_model))
            call mpas_log_write(' 4th order hyperviscosity is $r ', realArgs=(/h_mom_eddy_visc4/))
        end if

         delsq_w(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_w_euler(1:nVertLevels+1,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               edge_sign = 0.5 * r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

!DIR$ IVDEP                                                                                                                                           
              do k=2,nVertLevels

                  w_turb_flux =  edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))
                  delsq_w(k,iCell) = delsq_w(k,iCell) + w_turb_flux
                  w_turb_flux = w_turb_flux * meshScalingDel2(iEdge) * 0.25 * &
                                  ( eddy_visc_horz(k  ,cell1)+eddy_visc_horz(k  ,cell2)      &
                                   +eddy_visc_horz(k-1,cell1)+eddy_visc_horz(k-1,cell2) )
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + w_turb_flux
               end do
            end do
         end do

!$OMP BARRIER                                                                                                                                         

         if (h_mom_eddy_visc4 > 0.0) then  ! 4th order mixing is active                                                                               

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...                                                 
               r_areaCell = h_mom_eddy_visc4 * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell) * invDcEdge(iEdge)

                  do k=2,nVertLevels
                     tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - edge_sign * (delsq_w(k,cell2) - delsq_w(k,cell1))
                  end do

               end do
            end do

         end if ! 4th order mixing is active                                                                                                          

      if ( v_mom_eddy_visc2 > 0.0 ) then  ! vertical mixing

         do iCell=cellSolveStart,cellSolveEnd
!DIR$ IVDEP                                                                                                                                           
            do k=2,nVertLevels
               tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  &
                                           (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              &
                                          -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k)
            end do
         end do

      end if

         if ( config_les_model /= "none") then

             do iCell = cellSolveStart,cellSolveEnd ! vertical mixing for each column
                ! compute turbulent fluxes
                do k=1,nVertLevels
                  turb_vflux(k) = - rho_zz(k,iCell)*eddy_visc_vert(k,iCell)*zz(k,iCell)*zz(k,iCell)  &
                                    *rdzu(k)*(w(k+1,iCell)-w(k,iCell))
                end do
                do k=2,nVertLevels
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) &
                                              - rdzw(k)*(turb_vflux(k+1)-turb_vflux(k))
                end do
             end do

          end if

        if(debug_dissipation)  call mpas_log_write(' exiting w_dissipation_3d ')

      end subroutine w_dissipation_3d

!-----------------------------------------------------

      subroutine theta_dissipation( cellStart, cellEnd, cellSolveStart, cellSolveEnd,         &
                                    nCells, nEdges,                                           &
                                    nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, &
                                    invAreaCell, invDcEdge, dvEdge,                           &
                                    meshScalingDel2, meshScalingDel4,                         &
                                    config_mix_full, t_init, zgrid,                           &
                                    rdzw, rdzu,                                               &
                                    v_theta_eddy_visc2, h_theta_eddy_visc4, prandtl_inv,      &
                                    delsq_theta,                                              &
                                    theta_m, rho_edge, kdiff, rho_zz,                         &
                                    tend_theta_euler                                         )


      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             
      implicit none

      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: nCells, nEdges

      logical, intent(in) :: config_mix_full

      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: EdgesOnCell

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      real (kind=RKIND), intent(in) :: h_theta_eddy_visc4
      real (kind=RKIND), intent(in) :: v_theta_eddy_visc2
      real (kind=RKIND), intent(in) :: prandtl_inv

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu
      real (kind=RKIND), dimension(nVertLevels+1, nCells+1), intent(in) :: zgrid

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: t_init

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: kdiff
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta_euler

      ! storage passed in from calling routine                                                                                                        
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_theta

      ! local variables                                                                                                                               
      integer :: cell1, cell2, iEdge, iCell, i, k
      real (kind=RKIND) :: r_areaCell, edge_sign, theta_turb_flux, pr_scale
      real (kind=RKIND) :: z1, z2, z3, z4, zm, z0, zp

         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels

!  we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below

                  theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * pr_scale
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux

               end do
            end do
          end do

!$OMP BARRIER
            
         if (h_theta_eddy_visc4 > 0.0) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge)

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
                  end do
               end do
            end do

         end if ! 4th order mixing is active 

         if ( v_theta_eddy_visc2 > 0.0 ) then  ! vertical mixing for theta_m

            if (config_mix_full) then

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
                                             -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
                                             -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            end if

         end if

      end subroutine theta_dissipation

!-----------------------------------------------------

      subroutine theta_dissipation_3d( cellStart, cellEnd, cellSolveStart, cellSolveEnd,         &
                                       nCells, nEdges,                                           &
                                       nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, &
                                       invAreaCell, invDcEdge, dvEdge,                           &
                                       meshScalingDel2, meshScalingDel4,                         &
                                       config_mix_full, t_init, zgrid,                           &
                                       rdzw, rdzu, fzm, fzp,                                     &
                                       v_theta_eddy_visc2, h_theta_eddy_visc4, prandtl_inv,      &
                                       delsq_theta,                                              &
                                       theta_m, rho_edge, rho_zz, zz,                            &
                                       eddy_visc_horz, eddy_visc_vert,                           &
                                       config_les_model,                                         &
                                       tend_theta_euler                                         )


      ! 3D theta_m dissipation using the 3D smagorinsky eddy viscosities.  
      ! This routine also includes the simpler mixing models, and the 4th-order horizontal filter

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             
      implicit none

      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: nCells, nEdges

      logical, intent(in) :: config_mix_full

      character (len=StrKIND) :: config_les_model

      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: EdgesOnCell

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      real (kind=RKIND), intent(in) :: h_theta_eddy_visc4
      real (kind=RKIND), intent(in) :: v_theta_eddy_visc2
      real (kind=RKIND), intent(in) :: prandtl_inv

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzm
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzp
      real (kind=RKIND), dimension(nVertLevels+1, nCells+1), intent(in) :: zgrid
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: zz

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: t_init

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_horz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_vert
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta_euler

      ! storage passed in from calling routine                                                                                                        
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_theta

      ! local variables                                                                                                                               
      integer :: cell1, cell2, iEdge, iCell, i, k
      real (kind=RKIND) :: r_areaCell, edge_sign, theta_turb_flux, pr_scale
      real (kind=RKIND) :: z1, z2, z3, z4, zm, z0, zp
      real (kind=RKIND), dimension(nVertLevels+1) :: turb_vflux
      real (kind=RKIND) :: rho_k_at_w, zz_at_w

        if(debug_dissipation) then
            call mpas_log_write(' begin theta_dissipation_3d ')
            call mpas_log_write(' les model is '//trim(config_les_model))
            call mpas_log_write(' 4th order hyperviscosity is $r ', realArgs=(/h_theta_eddy_visc4/))
        end if

         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels

!  we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below.
!  This is in conservative form.

                  theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(eddy_visc_horz(k,cell1)+eddy_visc_horz(k,cell2)) * pr_scale
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux

               end do
            end do
          end do

!$OMP BARRIER
            
         if (h_theta_eddy_visc4 > 0.0) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge)

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
                  end do
               end do
            end do

         end if ! 4th order mixing is active 

         !  idealized case vertical mixing

         if ( v_theta_eddy_visc2 > 0.0 ) then  ! vertical mixing for theta_m

            if (config_mix_full) then

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
                                             -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
                                             -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            end if

         end if

         if ( config_les_model /= "none") then

             do iCell = cellSolveStart,cellSolveEnd ! vertical mixing for each column
                ! compute turbulent fluxes
                turb_vflux(nVertlevels+1) = 0.  !  no turbulent flux out of the domain
                turb_vflux(1) = 0.  !  lower bc flux handled where ???
                do k=2,nVertLevels
                  rho_k_at_w =  fzm(k)*rho_zz(k  ,iCell)*zz(k  ,iCell)*eddy_visc_vert(k  ,iCell) &
                               +fzp(k)*rho_zz(k-1,iCell)*zz(k-1,iCell)*eddy_visc_vert(k-1,iCell)
                  zz_at_w = fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell)
                  turb_vflux(k) = - rho_k_at_w*zz_at_w*rdzu(k)*(theta_m(k,iCell)-theta_m(k-1,iCell))
                end do

                ! test boundary conditions for supercell and les test cases

                if( les_test ) then
                  turb_vflux(1) = tke_heat_flux*rho_zz(1,iCell)  ! this is correct for DRY CASE ONLY
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                else
                  turb_vflux(1) = turb_vflux(2)
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                end if


                do k=1,nVertLevels
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) &
                                              - rdzw(k)*(turb_vflux(k+1)-turb_vflux(k))
                end do
             end do

         end if

        if(debug_dissipation)  call mpas_log_write(' exiting theta_dissipation_3d ')

      end subroutine theta_dissipation_3d

!-----------------------------------------------------

      subroutine scalar_dissipation_3d_les( cellStart, cellEnd, cellSolveStart, cellSolveEnd,         &
                                            nCells, nEdges,                                           &
                                            nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, &
                                            invAreaCell, invDcEdge, dvEdge,                           &
                                            meshScalingDel2, meshScalingDel4,                         &
                                            config_mix_full, t_init, zgrid,                           &
                                            rdzw, rdzu, fzm, fzp,                                     &
                                            v_theta_eddy_visc2, h_theta_eddy_visc4, prandtl_inv,      &
                                            delsq_theta,                                              &
                                            theta_m, rho_edge, rho_zz, zz,                            &
                                            eddy_visc_horz, eddy_visc_vert,                           &
                                            bv_freq2, config_len_disp, scalars, tend_scalars,         &
                                            index_tke, index_qv, num_scalars_dummy, mix_scalars,      &
                                            config_les_model,                                         &
                                            uReconstructZonal, uReconstructMeridional,                &
                                            tend_theta_euler                                         )


      ! 3D theta_m dissipation using the 3D smagorinsky eddy viscosities.  
      ! This routine also includes the simpler mixing models, and the 4th-order horizontal filter

      use mpas_atm_dimensions  !  pull nVertLevels and maxEdges from here                                                                             
      implicit none

      integer, intent(in) :: cellStart, cellEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd
      integer, intent(in) :: nCells, nEdges
      integer, intent(in) :: num_scalars_dummy
      integer, intent(in) :: index_tke, index_qv

      logical, intent(in) :: config_mix_full, mix_scalars

      character (len=StrKIND) :: config_les_model

      integer, dimension(nCells+1), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells+1), intent(in) :: EdgesOnCell

      integer, dimension(2,nEdges+1), intent(in) :: cellsOnEdge

      real (kind=RKIND), intent(in) :: h_theta_eddy_visc4
      real (kind=RKIND), intent(in) :: v_theta_eddy_visc2
      real (kind=RKIND), intent(in) :: prandtl_inv
      real (kind=RKIND), intent(in) :: config_len_disp

      real (kind=RKIND), dimension(maxEdges,nCells+1), intent(in) :: edgesOnCell_sign
      real (kind=RKIND), dimension(nCells+1), intent(in) :: invAreaCell
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: dvEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: invDcEdge
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel2
      real (kind=RKIND), dimension(nEdges+1), intent(in) :: meshScalingDel4
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzw
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: rdzu
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzm
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: fzp
      real (kind=RKIND), dimension(nVertLevels+1, nCells+1), intent(in) :: zgrid
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: zz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: uReconstructZonal, uReconstructMeridional

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: t_init

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: bv_freq2
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(in) :: scalars
      real (kind=RKIND), dimension(num_scalars,nVertLevels,nCells+1), intent(inout) :: tend_scalars
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: theta_m
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_horz
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: eddy_visc_vert
      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(in) :: rho_zz
      real (kind=RKIND), dimension(nVertLevels,nEdges+1), intent(in) :: rho_edge

      real (kind=RKIND), dimension(nVertLevels,nCells+1), intent(inout) :: tend_theta_euler

      ! storage passed in from calling routine                                                                                                        
      real (kind=RKIND), dimension(nVertLevels,nCells+1) :: delsq_theta

      ! local variables                                                                                                                               
      integer :: cell1, cell2, iEdge, iCell, i, k, iScalar
      real (kind=RKIND) :: r_areaCell, edge_sign, theta_turb_flux, pr_scale
      real (kind=RKIND) :: z1, z2, z3, z4, zm, z0, zp
      real (kind=RKIND), dimension(nVertLevels+1) :: turb_vflux, prandtl_inverse
      real (kind=RKIND), dimension(num_scalars,nVertLevels+1) :: turb_vflux_scalars
      real (kind=RKIND) :: rho_k_at_w, zz_at_w
      real (kind=RKIND) :: delta_z, delta_s, tke_length, bv_frequency2

      logical, parameter :: vmix_scalars = .true.
      logical, parameter :: hmix_scalars = .true.
      logical, parameter :: hmix4_scalars = .false.
      real (kind=RKIND) :: moisture_flux, heat_flux, theta_m_flux, c_h, c_q

        if(debug_dissipation) then
            call mpas_log_write(' begin scalar_dissipation_3d ')
            call mpas_log_write(' les model is '//trim(config_les_model))
            call mpas_log_write(' 4th order hyperviscosity is $r ', realArgs=(/h_theta_eddy_visc4/))
        end if

         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels

!  we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below.
!  This is in conservative form.

                  theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(eddy_visc_horz(k,cell1)+eddy_visc_horz(k,cell2)) * pr_scale
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux

               end do
            end do
          end do

!$OMP BARRIER
            
         if (h_theta_eddy_visc4 > 0.0) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge)

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
                  end do
               end do
            end do

         end if ! 4th order mixing is active 

       if(mix_scalars .and. hmix_scalars) then  ! dissipation for scalars, including 4th-order filter.  Likely needs optimization
         
         do iScalar=1,num_scalars
         
         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            ! tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = invAreaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * invDcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
               do k=1,nVertLevels

!  we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below.
!  This is in conservative form.

                  theta_turb_flux = edge_sign*(scalars(iScalar,k,cell2) - scalars(iScalar,k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(eddy_visc_horz(k,cell1)+eddy_visc_horz(k,cell2)) * pr_scale
                  tend_scalars(iScalar,k,iCell) = tend_scalars(iScalar,k,iCell) + theta_turb_flux

               end do
            end do
          end do

!$OMP BARRIER
            
         if ((h_theta_eddy_visc4 > 0.0).and. hmix4_scalars) then  ! 4th order mixing is active

            do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
               r_areaCell = h_theta_eddy_visc4 * prandtl_inv * invAreaCell(iCell)
               do i=1,nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i,iCell)
                  edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)*invDcEdge(iEdge)

                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  do k=1,nVertLevels
                     tend_scalars(iScalar,k,iCell) = tend_scalars(iScalar,k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
                  end do
               end do
            end do

         end if ! 4th order mixing is active 
         
        end do ! loop over scalars for horizontal mixing

      end if ! horizontal scalar mixing


         !  idealized case vertical mixing.  No scalar mixing here.

         if ( v_theta_eddy_visc2 > 0.0 ) then  ! vertical mixing for theta_m

            if (config_mix_full) then

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
                                             -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

               do iCell = cellSolveStart,cellSolveEnd
                  do k=2,nVertLevels-1
                     z1 = zgrid(k-1,iCell)
                     z2 = zgrid(k  ,iCell)
                     z3 = zgrid(k+1,iCell)
                     z4 = zgrid(k+2,iCell)

                     zm = 0.5*(z1+z2)
                     z0 = 0.5*(z2+z3)
                     zp = 0.5*(z3+z4)

                     tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                              ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
                                             -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
                  end do
               end do

            end if

         end if

         if ( config_les_model /= "none") then

             do iCell = cellSolveStart,cellSolveEnd ! vertical mixing for each column
                ! compute turbulent fluxes


                turb_vflux(nVertlevels+1) = 0.  !  no turbulent flux out of the domain
                turb_vflux(1) = 0.  !  lower bc flux handled where ???

                if ( config_les_model == "3d_smagorinsky") then
                  do k=2,nVertLevels
                    prandtl_inverse(k) = prandtl_inv
                  end do
                else   ! prognostic_1.5_order
                  do k=2,nVertLevels
                    delta_z = 0.5*(zgrid(k+1,iCell)-zgrid(k-1,iCell))
                    delta_s = ((config_len_disp**2)*delta_z)**(1./3.)
                    bv_frequency2 = 0.5*(bv_freq2(k,iCell)+bv_freq2(k-1,iCell))
                    tke_length = delta_s
                    if(bv_frequency2 .gt. 1.e-06) &
                        tke_length = 0.76*sqrt(scalars(index_tke,k,iCell))/sqrt(bv_frequency2)
                    prandtl_inverse(k) = 1. + 2.*tke_length/delta_z
                  end do
                end if

                do k=2,nVertLevels

                  ! delta_z = 0.5*(zgrid(k+1,iCell)-zgrid(k-1,iCell))
                  ! delta_s = ((config_len_disp**2)*delta_z)**(1./3.)
                  ! bv_frequency2 = 0.5*(bv_freq2(k)+bv_freq(k-1))
                  ! bv = max( sqrt(abs(bv_frequency2)), epsilon_bv )

                  rho_k_at_w =  fzm(k)*rho_zz(k  ,iCell)*zz(k  ,iCell)*eddy_visc_vert(k  ,iCell) &
                               +fzp(k)*rho_zz(k-1,iCell)*zz(k-1,iCell)*eddy_visc_vert(k-1,iCell)
                  zz_at_w = fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell)
                  turb_vflux(k) = - prandtl_inverse(k)*rho_k_at_w*zz_at_w*rdzu(k)*(theta_m(k,iCell)-theta_m(k-1,iCell))
                end do

                ! test boundary conditions for supercell and les test cases

                if( les_test ) then
                  moisture_flux = 0.
                  heat_flux = tke_heat_flux

                  ! bulk formulation will go here.

                  theta_m_flux = heat_flux*(1.0+(rv/rgas)*scalars(index_qv,1,iCell))  &
                                    +(rv/rgas)*theta_m(1,iCell)*moisture_flux/rho_zz(k,iCell)
                  ! turb_vflux(1) = tke_heat_flux*rho_zz(1,iCell)  ! this is correct for DRY CASE ONLY
                  turb_vflux(1) = theta_m_flux*rho_zz(1,iCell)
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                else
                  turb_vflux(1) = turb_vflux(2)
                  turb_vflux(nVertLevels+1) = turb_vflux(nVertLevels)
                end if


                do k=1,nVertLevels
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) &
                                              - rdzw(k)*(turb_vflux(k+1)-turb_vflux(k))
                end do

                if (mix_scalars .and. vmix_scalars) then

                ! compute turbulent fluxes
                turb_vflux_scalars(:,nVertlevels+1) = 0.  !  no turbulent flux out of the domain
                turb_vflux_scalars(:,1) = 0.  !  lower bc flux handled where ???
                do k=2,nVertLevels
                  rho_k_at_w =  fzm(k)*rho_zz(k  ,iCell)*zz(k  ,iCell)*eddy_visc_vert(k  ,iCell) &
                               +fzp(k)*rho_zz(k-1,iCell)*zz(k-1,iCell)*eddy_visc_vert(k-1,iCell)
                  zz_at_w = fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell)
                  do iScalar=1,num_scalars
                     turb_vflux_scalars(iScalar,k) = - prandtl_inverse(k)*rho_k_at_w*zz_at_w*rdzu(k)*  &
                                                      (scalars(iScalar,k,iCell)-scalars(iScalar,k-1,iCell))
                  end do
                end do

                if( les_test ) turb_vflux_scalars(index_qv,1) = moisture_flux ! lower b.c. for qv

                do k=1,nVertLevels
                  do iScalar=1,num_scalars
                    tend_scalars(iScalar,k,iCell) = tend_scalars(iScalar,k,iCell) &
                         - rdzw(k)*(turb_vflux_scalars(iScalar,k+1)-turb_vflux_scalars(iScalar,k))
                  end do
                end do

             end if ! mix scalars

             end do  ! loop over cells (columns)
             
         end if

        if(debug_dissipation)  call mpas_log_write(' exiting scalar_dissipation_3d ')

      end subroutine scalar_dissipation_3d_les

end module mpas_atm_dissipation_models
