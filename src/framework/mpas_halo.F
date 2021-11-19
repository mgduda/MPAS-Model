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

        call mpas_pool_create_pool(domain % haloGroupPool)

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

        call mpas_pool_destroy_pool(domain % haloGroupPool)

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
        use mpas_derived_types, only : domain_type, mpas_pool_type
        use mpas_pool_routines, only : mpas_pool_create_pool, mpas_pool_add_subpool

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_pool_type), pointer :: newGroup


        call mpas_pool_create_pool(newGroup)
	call mpas_pool_add_subpool(domain % haloGroupPool, groupName, newGroup)

        !
        ! TEMPORARY: pass through call to existing dmpar routine
        !
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

        use mpas_derived_types, only : domain_type, mpas_pool_type, mpas_pool_iterator_type, MPAS_POOL_CONFIG, &
                                       mpas_halo_group
        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_remove_subpool, mpas_pool_destroy_pool, &
                                       mpas_pool_begin_iteration, mpas_pool_get_next_member, mpas_pool_get_dimension, &
                                       mpas_pool_remove_field
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_pool_type), pointer :: completedGroup
        type (mpas_pool_iterator_type) :: itr
        integer, dimension(:), pointer :: fieldHaloInfo
        type (mpas_halo_group), pointer :: newGroup

        call mpas_pool_get_subpool(domain % haloGroupPool, groupName, completedGroup)

        call mpas_log_write('--- Creating halo exchange group '//trim(groupName)//' ---')

        !
        ! Add new mpas_halo_group to list of haloGroups in domain
        !
        allocate(newGroup)
        newGroup % groupName = groupName
        newGroup % next => domain % haloGroups
        domain % haloGroups => newGroup

        !
        ! Figure out how many fields are in this group
        !
        newGroup % nFields = 0
        call mpas_pool_begin_iteration(completedGroup)
        do while (mpas_pool_get_next_member(completedGroup, itr))
            if (itr % memberType == MPAS_POOL_CONFIG) then
                newGroup % nFields = newGroup % nFields + 1
            end if
        end do

        allocate(newGroup % fields(newGroup % nFields))

        !
        ! Fill in field entries for this group
        !
        call mpas_pool_begin_iteration(completedGroup)
        do while (mpas_pool_get_next_member(completedGroup, itr))
            if (itr % memberType == MPAS_POOL_CONFIG) then
                call mpas_log_write('    '//trim(itr % memberName))
                call mpas_pool_get_dimension(completedGroup, trim(itr % memberName)//'.info', fieldHaloInfo)
                call mpas_log_write('     - fieldType = $i',   intArgs=[fieldHaloInfo(1)])
                call mpas_log_write('     - nDims = $i',       intArgs=[fieldHaloInfo(2)])
                call mpas_log_write('     - nTimeLevels = $i', intArgs=[fieldHaloInfo(3)])
                call mpas_log_write('     - nHaloLayers = $i', intArgs=[fieldHaloInfo(4)])
                call mpas_pool_remove_field(completedGroup, trim(itr % memberName)//'.field')
            end if
        end do

        call mpas_log_write('')

        call mpas_pool_destroy_pool(completedGroup)
        call mpas_pool_remove_subpool(domain % haloGroupPool, groupName)

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
        use mpas_derived_types, only : domain_type, mpas_halo_group, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_halo_group), pointer :: cursor, prev


        !
        ! Find this halo exhange group in the list of groups
        !
        nullify(prev)
        cursor => domain % haloGroups
        do while (associated(cursor))
            if (trim(cursor % groupName) == trim(groupName)) then
                exit
            end if
           
            prev => cursor
            cursor => cursor % next
        end do

        if (.not. associated(cursor)) then
            call mpas_log_write('Halo exchange group '//trim(groupName)//' not found in destroy routine.', &
                                messageType=MPAS_LOG_CRIT)
        end if

        !
        ! Unlink this exchange group
        !
        if (.not. associated(prev)) then
            domain % haloGroups => cursor % next
        else
            prev % next => cursor % next
        end if

        call mpas_log_write('--- Destroying halo exchange group '//trim(groupName)//' ---')

        !
        ! Deallocate this exchange group
        !
        deallocate(cursor % fields)


        !
        ! TEMPORARY: pass through call to existing dmpar routine
        !
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
        use mpas_derived_types, only : domain_type, mpas_pool_type, mpas_pool_field_info_type, MPAS_POOL_REAL, &
                                       field2DReal, field3DReal, MPAS_LOG_CRIT
        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_add_config, mpas_pool_get_field_info, &
                                       mpas_pool_add_dimension, mpas_pool_get_field, mpas_pool_add_field
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        character (len=*), intent(in) :: fieldName
        integer, optional, intent(in) :: timeLevel
        integer, dimension(:), optional, intent(in) :: haloLayers
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_pool_type), pointer :: group
        type (mpas_pool_field_info_type) :: info
        integer, allocatable, dimension(:) :: fieldHaloInfo
        type (field2DReal), pointer :: r2d
        type (field3DReal), pointer :: r3d


        call mpas_pool_get_subpool(domain % haloGroupPool, groupName, group)

        !
        ! Store an item in the pool to signal the field whose halo is to be exchanged
        !
        call mpas_pool_add_config(group, fieldName, 1)

        call mpas_pool_get_field_info(domain % blocklist % allFields, fieldName, info)

        !
        ! Store an item in the pool with basic info about this field
        !
        allocate(fieldHaloInfo(4))
        fieldHaloInfo(1) = info % fieldType
        fieldHaloInfo(2) = info % nDims
        fieldHaloInfo(3) = info % nTimeLevels
        fieldHaloInfo(4) = info % nHaloLayers
        call mpas_pool_add_dimension(group, fieldName//'.info', fieldHaloInfo)
        deallocate(fieldHaloInfo)

        !
        ! Store an item in the pool with list of halo layers to exchange, or (/-1/) if all layers
        !
        if (present(haloLayers)) then
            call mpas_pool_add_dimension(group, fieldName//'.layers', haloLayers)
        else
            call mpas_pool_add_dimension(group, fieldName//'.layers', (/ -1 /))
        end if

        !
        ! Store an item in the pool indicating which time level to exchange
        !
        if (present(timeLevel)) then
            call mpas_pool_add_dimension(group, fieldName//'.timelevel', timeLevel)
        else
            call mpas_pool_add_dimension(group, fieldName//'.timelevel', 1)
        end if

        !
        ! Store a reference to the field itself in the pool
        !
        if (info % fieldType == MPAS_POOL_REAL) then
            if (info % nDims == 2) then
                call mpas_pool_get_field(domain % blocklist % allFields, fieldName, r2d)
                call mpas_pool_add_field(group, fieldName//'.field', r2d)
            else if (info % nDims == 3) then
                call mpas_pool_get_field(domain % blocklist % allFields, fieldName, r3d)
                call mpas_pool_add_field(group, fieldName//'.field', r3d)
            else
                call mpas_log_write('Unsupported dimensionality for real field '//trim(fieldName), messageType=MPAS_LOG_CRIT)
            end if
        else
            call mpas_log_write('Unsupported field type for field '//trim(fieldName), messageType=MPAS_LOG_CRIT)
        end if


        !
        ! TEMPORARY: pass through call to existing dmpar routine
        !
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
        use mpas_derived_types, only : domain_type, mpas_halo_group, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_halo_group), pointer :: cursor


        !
        ! Find this halo exhange group in the list of groups
        !
        cursor => domain % haloGroups
        do while (associated(cursor))
            if (trim(cursor % groupName) == trim(groupName)) then
                exit
            end if
           
            cursor => cursor % next
        end do

        if (.not. associated(cursor)) then
            call mpas_log_write('Halo exchange group '//trim(groupName)//' not found in full_exch routine.', &
                                messageType=MPAS_LOG_CRIT)
        end if


        !
        ! TEMPORARY: pass through call to existing dmpar routine
        !
        call mpas_dmpar_exch_group_full_halo_exch(domain, groupName, iErr)

    end subroutine mpas_halo_exch_group_full_halo_exch

end module mpas_halo
