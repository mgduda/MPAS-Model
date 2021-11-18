! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_halo
!
!> \brief Routines for updating field halos
!> \author Michael Duda
!> \date   29 September 2021
!> \details
!>  This module provides routines for updating the halos of fields.
!
!-----------------------------------------------------------------------
module mpas_halo

    implicit none

    private

    public :: mpas_halo_exch_group_init, &
              mpas_halo_exch_group_finalize, &
              mpas_halo_exch_group_create, &
              mpas_halo_exch_group_complete, &
              mpas_halo_exch_group_destroy, &
              mpas_halo_exch_group_add_field, &
              mpas_halo_exch_group_full_halo_exch


    contains


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_init
    !
    !> \brief Initialize group halo exchanges
    !> \author Michael Duda
    !> \date   17 November 2021
    !> \details
    !>  This routine initializes group halo exchanges, and must be called
    !>  before any other routine for building or exchanging halos.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_init(domain, iErr)

        use mpas_derived_types, only : domain_type
        use mpas_pool_routines, only : mpas_pool_create_pool

        type (domain_type), intent(inout) :: domain
        integer, optional, intent(out) :: iErr


        if (present(iErr)) then
           iErr = 0
        end if

        call mpas_pool_create_pool(domain % haloGroups)

    end subroutine mpas_halo_exch_group_init


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_finalize
    !
    !> \brief Finalize group halo exchanges
    !> \author Michael Duda
    !> \date   17 November 2021
    !> \details
    !>  This routine finalize group halo exchanges, and must be called
    !>  after all other calls for building or exchanging halos.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_finalize(domain, iErr)

        use mpas_derived_types, only : domain_type
        use mpas_pool_routines, only : mpas_pool_destroy_pool

        type (domain_type), intent(inout) :: domain
        integer, optional, intent(out) :: iErr


        if (present(iErr)) then
           iErr = 0
        end if

        call mpas_pool_destroy_pool(domain % haloGroups)

    end subroutine mpas_halo_exch_group_finalize


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_create
    !
    !> \brief MPAS halo exchange group creation routine
    !> \author Doug Jacobsen
    !> \date   01/05/2016
    !> \details
    !>  This routine creates a new exchange group named 'groupName' within a
    !>  domain's list of exchange groups.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_create(domain, groupName, iErr)

        use mpas_dmpar, only : mpas_dmpar_exch_group_create
        use mpas_derived_types, only : domain_type

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr


        call mpas_dmpar_exch_group_create(domain, groupName, iErr)

    end subroutine mpas_halo_exch_group_create


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_complete
    !
    !> \brief Complete the creation of an exchange group
    !> \author Michael Duda
    !> \date   29 September 2021
    !> \details
    !>  Complete the creation of an exchange group that was defined via a call
    !>  to the mpas_halo_exch_group_create routine, and to which field were
    !>  added through calls to mpas_halo_exch_group_add_field. This routine
    !>  must be called for an exchange group before the group can be used in
    !>  calls to the mpas_halo_exch_group_full_halo_exch routine.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_complete(domain, groupName, iErr)

        use mpas_derived_types, only : domain_type

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr


    end subroutine mpas_halo_exch_group_complete


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_destroy
    !
    !> \brief MPAS halo exchange group destruction routine
    !> \author Doug Jacobsen
    !> \date   01/05/2016
    !> \details
    !>  This routine destroys an exchange group named 'groupName' within a
    !>  domain's list of exchange groups.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_destroy(domain, groupName, iErr)

        use mpas_dmpar, only : mpas_dmpar_exch_group_destroy
        use mpas_derived_types, only : domain_type

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr


        call mpas_dmpar_exch_group_destroy(domain, groupName, iErr)

    end subroutine mpas_halo_exch_group_destroy


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_add_field
    !
    !> \brief MPAS halo exchange group add field routine
    !> \author Doug Jacobsen
    !> \date   01/05/2016
    !> \details
    !>  This routine adds a field named 'fieldName' to an exchange group named
    !>  'groupName' within a domain's list of exchange groups.
    !>  The timeLevel argument allows control over which timeLevel(s) will be
    !>    exchanged as part of this group. If the timeLevel argument is omitted or if it
    !>    has a value of -1, all time levels will be exchanged.
    !>  The haloLayers argument allows an input array to define the halo layers
    !>    that should be exchanged as part of this exchange group. If it is
    !>    omitted, all halo layers will be exchanged.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_add_field(domain, groupName, fieldName, timeLevel, haloLayers, iErr)

        use mpas_dmpar, only : mpas_dmpar_exch_group_add_field
        use mpas_derived_types, only : domain_type

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        character (len=*), intent(in) :: fieldName
        integer, optional, intent(in) :: timeLevel
        integer, dimension(:), optional, intent(in) :: haloLayers
        integer, optional, intent(out) :: iErr


        call mpas_dmpar_exch_group_add_field(domain, groupName, fieldName, timeLevel, haloLayers, iErr)

    end subroutine mpas_halo_exch_group_add_field


    !-----------------------------------------------------------------------
    !  routine mpas_halo_exch_group_full_halo_exch
    !
    !> \brief MPAS halo exchange group full halo exchange routine
    !> \author Doug Jacobsen
    !> \date   01/11/2016
    !> \details
    !>  This routine performs a full halo exchange on an exchange group.
    !>  It is blocking, in that the routine doesn't return until the full
    !>  exchange is complete.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_exch_group_full_halo_exch(domain, groupName, iErr)

        use mpas_dmpar, only : mpas_dmpar_exch_group_full_halo_exch
        use mpas_derived_types, only : domain_type

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr


        call mpas_dmpar_exch_group_full_halo_exch(domain, groupName, iErr)

    end subroutine mpas_halo_exch_group_full_halo_exch

end module mpas_halo
