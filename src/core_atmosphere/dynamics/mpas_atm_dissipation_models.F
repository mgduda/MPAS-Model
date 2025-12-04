! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

module mpas_atm_dissipation_models

  use mpas_kind_types, only : RKIND

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

      end subroutine smagorinsky_2d

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
               !                    only valid for h_mom_eddy_visc4 == constant                                                                       
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
                                           - v_init(k) * sin( angleEdge(iEdge) )
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

end module mpas_atm_dissipation_models
