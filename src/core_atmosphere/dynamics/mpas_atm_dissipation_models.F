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

end module mpas_atm_dissipation_models
