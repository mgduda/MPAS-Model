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

        use mpas_derived_types, only : domain_type, mpas_pool_type
        use mpas_pool_routines, only : mpas_pool_create_pool, mpas_pool_add_subpool

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        type (mpas_pool_type), pointer :: newGroup


        call mpas_pool_create_pool(newGroup)
	call mpas_pool_add_subpool(domain % haloGroupPool, groupName, newGroup)

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
                                       MPAS_POOL_REAL, MPAS_HALO_REAL, mpas_halo_group, MPAS_LOG_CRIT, &
                                       field2DReal, field3DReal
        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_remove_subpool, mpas_pool_destroy_pool, &
                                       mpas_pool_begin_iteration, mpas_pool_get_next_member, mpas_pool_get_dimension, &
                                       mpas_pool_remove_field, mpas_pool_get_field
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        integer :: i
        type (mpas_pool_type), pointer :: completedGroup
        type (mpas_pool_iterator_type) :: itr
        integer, dimension(:), pointer :: fieldHaloInfo
        integer, dimension(:), pointer :: haloLayers
        type (mpas_halo_group), pointer :: newGroup
        type (field2DReal), pointer :: r2d
        type (field3DReal), pointer :: r3d
        integer, pointer :: timeLevel


        call mpas_pool_get_subpool(domain % haloGroupPool, groupName, completedGroup)

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
        i = 1
        call mpas_pool_begin_iteration(completedGroup)
        do while (mpas_pool_get_next_member(completedGroup, itr))
            if (itr % memberType == MPAS_POOL_CONFIG) then
                newGroup % fields(i) % fieldName = trim(itr % memberName)

                call mpas_pool_get_dimension(completedGroup, trim(itr % memberName)//'.info', fieldHaloInfo)

                call mpas_pool_get_dimension(completedGroup, trim(itr % memberName)//'.timelevel', timeLevel)

                call mpas_pool_get_dimension(completedGroup, trim(itr % memberName)//'.layers', haloLayers)

                newGroup % fields(i) % nDims = fieldHaloInfo(2)
                newGroup % fields(i) % timeLevel = timeLevel

                select case (fieldHaloInfo(1))
                case (MPAS_POOL_REAL)
                    newGroup % fields(i) % fieldType = MPAS_HALO_REAL
                case default
                    call mpas_log_write('Only real-valued fields are supported in mpas_halo_exch_group_complete', &
                                        messageType=MPAS_LOG_CRIT)
                end select

                if (fieldHaloInfo(1) == MPAS_POOL_REAL) then
                    if (fieldHaloInfo(2) == 2) then
                        call mpas_pool_get_field(completedGroup, trim(itr % memberName)//'.field', r2d)
                        call mpas_halo_compact_halo_info(domain, r2d % sendList, r2d % recvList, r2d % dimSizes, &
                                                         newGroup % fields(i) % compactHaloInfo, haloLayers)
                    else if (fieldHaloInfo(2) == 3) then
                        call mpas_pool_get_field(completedGroup, trim(itr % memberName)//'.field', r3d)
                        call mpas_halo_compact_halo_info(domain, r3d % sendList, r3d % recvList, r3d % dimSizes, &
                                                         newGroup % fields(i) % compactHaloInfo, haloLayers)
                    else
                        call mpas_log_write('Unsupported dimensionality for real field in mpas_halo_exch_group_complete.', &
                                            messageType=MPAS_LOG_CRIT)
                    end if
                else
                    call mpas_log_write('Unsupported field type in mpas_halo_exch_group_complete.', &
                                        messageType=MPAS_LOG_CRIT)
                end if

                call mpas_pool_remove_field(completedGroup, trim(itr % memberName)//'.field')
                i = i + 1
            end if
        end do

        call mpas_halo_aggregate_group_info(newGroup)

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

        use mpas_derived_types, only : domain_type, mpas_halo_group, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        integer :: i
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

        !
        ! Deallocate this exchange group
        !
        do i = 1, cursor % nFields
            deallocate(cursor % fields(i) % compactHaloInfo)
        end do
        deallocate(cursor % fields)
        deallocate(cursor % groupPackOffsets)
        deallocate(cursor % groupSendNeighbors)
        deallocate(cursor % groupSendOffsets)
        deallocate(cursor % groupSendCounts)
        deallocate(cursor % groupUnpackOffsets)
        deallocate(cursor % groupRecvNeighbors)
        deallocate(cursor % groupRecvOffsets)
        deallocate(cursor % groupRecvCounts)
        deallocate(cursor)

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
        integer :: local_timeLevel
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
            local_timeLevel = timeLevel
        else
            local_timeLevel = 1
        end if
        call mpas_pool_add_dimension(group, fieldName//'.timelevel', local_timeLevel)

        !
        ! Store a reference to the field itself in the pool
        !
        if (info % fieldType == MPAS_POOL_REAL) then
            if (info % nDims == 2) then
                call mpas_pool_get_field(domain % blocklist % allFields, fieldName, r2d, timeLevel=local_timeLevel)
                call mpas_pool_add_field(group, fieldName//'.field', r2d)
            else if (info % nDims == 3) then
                call mpas_pool_get_field(domain % blocklist % allFields, fieldName, r3d, timeLevel=local_timeLevel)
                call mpas_pool_add_field(group, fieldName//'.field', r3d)
            else
                call mpas_log_write('Unsupported dimensionality for real field '//trim(fieldName), messageType=MPAS_LOG_CRIT)
            end if
        else
            call mpas_log_write('Unsupported field type for field '//trim(fieldName), messageType=MPAS_LOG_CRIT)
        end if

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

        use mpas_kind_types, only : RKIND
        use mpas_derived_types, only : domain_type, mpas_halo_group, MPAS_HALO_REAL, MPAS_LOG_CRIT
        use mpas_pool_routines, only : mpas_pool_get_array
        use mpas_log, only : mpas_log_write
        use mpi, only : MPI_Irecv, MPI_Wait, MPI_REAL, MPI_STATUS_IGNORE, MPI_SUCCESS

        type (domain_type), intent(inout) :: domain
        character (len=*), intent(in) :: groupName
        integer, optional, intent(out) :: iErr

        ! Local variables
        integer :: i, bufstart, bufend
        integer :: dim1, dim2, dim3
        integer :: i1, i2, i3, iNeighbor, j
        integer :: iHalo, iEndp, iOffset
        integer :: nSendEndpts, nRecvEndpts
        integer :: idx, nList
        integer :: mpi_ierr
        type (mpas_halo_group), pointer :: group
        integer, dimension(:), pointer :: compactHaloInfo


        !
        ! Find this halo exhange group in the list of groups
        !
        group => domain % haloGroups
        do while (associated(group))
            if (trim(group % groupName) == trim(groupName)) then
                exit
            end if
           
            group => group % next
        end do

        if (.not. associated(group)) then
            call mpas_log_write('Halo exchange group '//trim(groupName)//' not found in full_exch routine.', &
                                messageType=MPAS_LOG_CRIT)
        end if

        allocate(group % sendBuf(group % groupSendBufSize))
        allocate(group % recvBuf(group % groupRecvBufSize))
        allocate(group % sendRequests(group % nGroupSendNeighbors))
        allocate(group % recvRequests(group % nGroupRecvNeighbors))


        do i = 1, group % nGroupRecvNeighbors
            if (group % groupRecvCounts(i) > 0) then
                bufstart = group % groupRecvOffsets(i)
                bufend = group % groupRecvOffsets(i) + group % groupRecvCounts(i) - 1
!TO DO: how do we determine appropriate type here?
                call MPI_Irecv(group % recvBuf(bufstart:bufend), group % groupRecvCounts(i), MPI_REAL, &
                               group % groupRecvNeighbors(i), group % groupRecvNeighbors(i), domain % dminfo % comm, &
                               group % recvRequests(i), mpi_ierr)
            end if
        end do

        do i = 1, group % nFields

            compactHaloInfo => group % fields(i) % compactHaloInfo

            if (group % fields(i) % fieldType == MPAS_HALO_REAL) then
                dim1 = compactHaloInfo(2)
                dim2 = compactHaloInfo(3)
                dim3 = compactHaloInfo(4)

                nSendEndpts = compactHaloInfo(10)
                idx = 11

                select case (group % fields(i) % nDims)
                case (2)
                    call mpas_pool_get_array(domain % blocklist % allFields, trim(group % fields(i) % fieldName), &
                                             group % fields(i) % r2arr, timeLevel=group % fields(i) % timeLevel)

                    do iEndp = 1, nSendEndpts
                        ! Find this endpoint in the list of neighbors
                        do iNeighbor = 1, group % nGroupSendNeighbors
                            if (group % groupSendNeighbors(iNeighbor) == compactHaloInfo(idx)) exit
                        end do
                        idx = idx + 1

                        iOffset = group % groupPackOffsets(iNeighbor, i)

                        do iHalo = 1, compactHaloInfo(9)
                            nList = compactHaloInfo(idx)
                            idx = idx + 1

                            do j = 1, nList
                                do i1 = 1, dim1
                                    group % sendBuf(iOffset + dim1 * (compactHaloInfo(idx + nList) - 1) + i1) = &
                                        group % fields(i) % r2arr(i1, compactHaloInfo(idx))
                                end do
                                idx = idx + 1
                            end do
                            idx = idx + nList
                        end do
                    end do

                case (3)

                    call mpas_pool_get_array(domain % blocklist % allFields, trim(group % fields(i) % fieldName), &
                                             group % fields(i) % r3arr, group % fields(i) % timeLevel)

                    do iEndp = 1, compactHaloInfo(10)
                        ! Find this endpoint in the list of neighbors
                        do iNeighbor = 1, group % nGroupSendNeighbors
                            if (group % groupSendNeighbors(iNeighbor) == compactHaloInfo(idx)) exit
                        end do
                        idx = idx + 1

                        iOffset = group % groupPackOffsets(iNeighbor, i)

                        do iHalo = 1, compactHaloInfo(9)
                            nList = compactHaloInfo(idx)
                            idx = idx + 1

                            do j = 1, nList
                                do i2 = 1, dim2
                                do i1 = 1, dim1
                                    group % sendBuf(iOffset + dim1 * dim2 * (compactHaloInfo(idx + nList) - 1) + dim1 * (i2 - 1) + i1) = &
                                        group % fields(i) % r3arr(i1, i2, compactHaloInfo(idx))
                                end do
                                end do
                                idx = idx + 1
                            end do
                            idx = idx + nList
                        end do
                    end do
                case default
                    call mpas_log_write('Only 2-d and 3-d fields are supported by mpas_halo_exch_group_full_halo_exch', &
                                    messageType=MPAS_LOG_CRIT)
                end select
            else
                call mpas_log_write('Only real-valued fields are supported by mpas_halo_exch_group_full_halo_exch', &
                                    messageType=MPAS_LOG_CRIT)
            end if
        end do


        do i = 1, group % nGroupSendNeighbors
            if (group % groupSendCounts(i) > 0) then
                bufstart = group % groupSendOffsets(i)
                bufend = group % groupSendOffsets(i) + group % groupSendCounts(i) - 1
!TO DO: how do we determine appropriate type here?
                call MPI_Isend(group % sendBuf(bufstart:bufend), group % groupSendCounts(i), MPI_REAL, &
                               group % groupSendNeighbors(i), domain % dminfo % my_proc_id, domain % dminfo % comm, &
                               group % sendRequests(i), mpi_ierr)
            end if
        end do

        do i = 1, group % nGroupRecvNeighbors
            if (group % groupRecvCounts(i) > 0) then
                call MPI_Wait(group % recvRequests(i), MPI_STATUS_IGNORE, mpi_ierr)
            end if
        end do

        do i = 1, group % nFields

            compactHaloInfo => group % fields(i) % compactHaloInfo

            if (group % fields(i) % fieldType == MPAS_HALO_REAL) then
                dim1 = compactHaloInfo(2)
                dim2 = compactHaloInfo(3)
                dim3 = compactHaloInfo(4)

                nSendEndpts = compactHaloInfo(10)
                idx = 11
                do iendp = 1, nSendEndpts
                    idx = idx + 1                ! skip over endPointID
                    do ihalo = 1, compactHaloInfo(9)
                        nlist = compactHaloInfo(idx)
                        idx = idx + 1            ! skip over nList
                        idx = idx + 2 * nlist    ! skip over srcList and destList
                    end do
                end do

                nRecvEndpts = compactHaloInfo(idx)
                idx = idx + 1

                select case (group % fields(i) % nDims)
                case (2)
                    do iEndp = 1, nRecvEndpts
                        ! Find this endpoint in the list of neighbors
                        do iNeighbor = 1, group % nGroupRecvNeighbors
                            if (group % groupRecvNeighbors(iNeighbor) == compactHaloInfo(idx)) exit
                        end do
                        idx = idx + 1

                        iOffset = group % groupUnpackOffsets(iNeighbor, i)

                        do iHalo = 1, compactHaloInfo(9)
                            nList = compactHaloInfo(idx)
                            idx = idx + 1

                            do j = 1, nList
                                do i1 = 1, dim1
                                    group % fields(i) % r2arr(i1, compactHaloInfo(idx + nList)) = &
                                        group % recvBuf(iOffset + dim1 * (compactHaloInfo(idx) - 1) + i1)
                                end do
                                idx = idx + 1
                            end do
                            idx = idx + nList
                        end do
                    end do

                case (3)
                    do iEndp = 1, nRecvEndpts
                        ! Find this endpoint in the list of neighbors
                        do iNeighbor = 1, group % nGroupRecvNeighbors
                            if (group % groupRecvNeighbors(iNeighbor) == compactHaloInfo(idx)) exit
                        end do
                        idx = idx + 1

                        iOffset = group % groupUnpackOffsets(iNeighbor, i)

                        do iHalo = 1, compactHaloInfo(9)
                            nList = compactHaloInfo(idx)
                            idx = idx + 1

                            do j = 1, nList
                                do i2 = 1, dim2
                                do i1 = 1, dim1
                                    group % fields(i) % r3arr(i1, i2, compactHaloInfo(idx + nList)) = &
                                        group % recvBuf(iOffset + dim1 * dim2 * (compactHaloInfo(idx) - 1) + dim1 * (i2 - 1) + i1)
                                end do
                                end do
                                idx = idx + 1
                            end do
                            idx = idx + nList
                        end do
                    end do
                case default
                    call mpas_log_write('Only 2-d and 3-d fields are supported by mpas_halo_exch_group_full_halo_exch', &
                                    messageType=MPAS_LOG_CRIT)
                end select
            else
                call mpas_log_write('Only real-valued fields are supported by mpas_halo_exch_group_full_halo_exch', &
                                    messageType=MPAS_LOG_CRIT)
            end if
        end do

        do i = 1, group % nGroupSendNeighbors
            if (group % groupSendCounts(i) > 0) then
                call MPI_Wait(group % sendRequests(i), MPI_STATUS_IGNORE, mpi_ierr)
            end if
        end do

        deallocate(group % sendBuf)
        deallocate(group % recvBuf)
        deallocate(group % sendRequests)
        deallocate(group % recvRequests)

    end subroutine mpas_halo_exch_group_full_halo_exch


    !-----------------------------------------------------------------------
    !  routine mpas_halo_compact_halo_info
    !
    !> \brief Compacts information needed for halo exchanges
    !> \author Michael Duda
    !> \date 7 December 2017
    !> \details
    !>  This routine extracts all information needed to perform a halo exchange
    !>  from dynamic data types and places it into a single, contiguous array
    !>  for use by the mpas_halo_exch_halo_acc routines.
    !>  The resulting compactHaloInfo array has the following elements:
    !>    1 - The dimensionality of the field
    !>    2 - Dimension 1 of the field (i.e., the left-most dimension)
    !>    3 - Dimension 2 of the field
    !>    4 - Dimension 3 of the field
    !>    5 - Dimension 4 of the field
    !>    6 - Dimension 5 of the field
    !>    7 - The MPI communicator
    !>    8 - The MPI rank of the current process
    !>    9 - The number of halo layers for the field
    !>   10 - The number of endpoints to send to
    !>   11 - foreach (send endpoint) { endPointID foreach (halolayer) {nList srcList(1:nList) destList(1:nList)} }
    !>   12 - The number of endpoints to recv from
    !>   13 - foreach (recv endpoint) { endPointID foreach (halolayer) {nList srcList(1:nList) destList(1:nList)} }
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_compact_halo_info(domain, sendList, recvList, dimSizes, compactHaloInfo, haloLayers)

        use mpas_derived_types, only : domain_type, mpas_multihalo_exchange_list, mpas_exchange_list, MPAS_LOG_CRIT
        use mpas_log, only : mpas_log_write

        implicit none

        type (domain_type), intent(in) :: domain
        type (mpas_multihalo_exchange_list), pointer :: sendList
        type (mpas_multihalo_exchange_list), pointer :: recvList
        integer, dimension(:), intent(in) :: dimsizes
        integer, dimension(:), pointer :: compactHaloInfo
        integer, dimension(:), intent(in) :: haloLayers

        ! Local variables
        integer :: infosize
        integer :: i, iendpt, j, ioffset
        integer :: idx
        integer :: nSendEndpoints
        integer :: maxSendEndpoints
        integer :: totSendListSize
        integer :: nRecvEndpoints
        integer :: maxRecvEndpoints
        integer :: totRecvListSize
        integer :: nHaloLayers
        integer :: nEndpoints
        logical :: found
        integer, allocatable, dimension(:) :: sendEndpoints
        integer, allocatable, dimension(:) :: recvEndpoints
        type (mpas_multihalo_exchange_list), pointer :: sendListCursor, recvListCursor
        type (mpas_exchange_list), pointer :: exchListPtr
        integer :: activeHaloLayers
        logical, dimension(:), allocatable :: useHalo


        !
        ! Find number of halo layers
        !
        nHaloLayers = size(sendList % halos)
        if (nHaloLayers /= size(recvList % halos)) then
            call mpas_log_write('The number of halo layers in the recv list does not match the number in the send list', &
                                messageType=MPAS_LOG_CRIT)
        end if

        !
        ! Create logical array indicating, for each halo layer, whether that halo layer should be
        ! used in a halo exchange
        !
        allocate(useHalo(nHaloLayers))

        if (haloLayers(1) == -1) then     ! Use all halo layers
            useHalo(:) = .true.
            activeHaloLayers = nHaloLayers
        else
            useHalo(:) = .false.
            activeHaloLayers = size(haloLayers)
            do i = 1, activeHaloLayers
                useHalo(haloLayers(i)) = .true.
            end do
        end if

        !
        ! Find the maximum number of "endpoints" that we will need to send to for any halo layer,
        ! as well as the total size of the send lists for all halo layers
        !
        maxSendEndpoints = 0
        totSendListSize = 0
        sendListCursor => sendList
        do while (associated(sendListCursor))
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    nSendEndpoints = 0
                    exchListPtr => sendListCursor % halos(i) % exchList
                    do while (associated(exchListPtr))
                        nSendEndpoints = nSendEndpoints + 1
                        totSendListSize = totSendListSize + 2 * exchListPtr % nList    ! We have srcList and destList
                        exchListPtr => exchListPtr % next
                    end do
                    maxSendEndpoints = max(nSendEndpoints, maxSendEndpoints)
                end if
            end do
            sendListCursor => sendListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
        end do

        allocate(sendEndpoints(maxSendEndpoints * activeHaloLayers))

        !
        ! Gather a list of the unique endpoints that we will need to send to
        !
        nSendEndpoints = 0
        sendListCursor => sendList
        do while (associated(sendListCursor))
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    exchListPtr => sendListCursor % halos(i) % exchList
                    do while (associated(exchListPtr))
 
                        ! If the current endpoint is not already in the list, add it
                        do j=1,nSendEndpoints
                            if (exchListPtr % endPointID == sendEndpoints(j)) exit
                        end do
                        if (j > nSendEndpoints) then
                            nSendEndpoints = nSendEndpoints + 1
                            sendEndpoints(nSendEndpoints) = exchListPtr % endPointID
                        end if

                        exchListPtr => exchListPtr % next
                    end do
                end if
            end do
            sendListCursor => sendListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
        end do

        !
        ! Find the maximum number of "endpoints" that we will need to receive from for any halo layer,
        ! as well as the total size of the recv lists for all halo layers
        !
        maxRecvEndpoints = 0
        totRecvListSize = 0
        recvListCursor => recvList
        do while (associated(recvListCursor))
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    nRecvEndpoints = 0
                    exchListPtr => recvListCursor % halos(i) % exchList
                    do while (associated(exchListPtr))
                        nRecvEndpoints = nRecvEndpoints + 1
                        totRecvListSize = totRecvListSize + 2 * exchListPtr % nList    ! We have srcList and destList
                        exchListPtr => exchListPtr % next
                    end do
                    maxRecvEndpoints = max(nRecvEndpoints, maxRecvEndpoints)
                end if
            end do
            recvListCursor => recvListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
        end do

        allocate(recvEndpoints(maxRecvEndpoints * activeHaloLayers))

        !
        ! Gather a list of the unique endpoints that we will need to receive from
        !
        nRecvEndpoints = 0
        recvListCursor => recvList
        do while (associated(recvListCursor))
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    exchListPtr => recvListCursor % halos(i) % exchList
                    do while (associated(exchListPtr))

                        ! If the current endpoint is not already in the list, add it
                        do j=1,nRecvEndpoints
                            if (exchListPtr % endPointID == recvEndpoints(j)) exit
                        end do
                        if (j > nRecvEndpoints) then
                            nRecvEndpoints = nRecvEndpoints + 1
                            recvEndpoints(nRecvEndpoints) = exchListPtr % endPointID
                        end if

                        exchListPtr => exchListPtr % next
                    end do
                end if
            end do
            recvListCursor => recvListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
        end do


        !
        ! Compute the number of elements we will need in compactHaloInfo
        !
        infosize = 11 + nSendEndpoints + activeHaloLayers * nSendEndpoints + totSendListSize &
                      + nRecvEndpoints + activeHaloLayers * nRecvEndpoints + totRecvListSize

        allocate(compactHaloInfo(infosize))
        compactHaloInfo(:) = 0

        !
        ! 1-6: Add field dimensionality and dimensions
        !
        compactHaloInfo(1) = size(dimSizes)    ! Dimensionality of the field
        idx = 2
        do i=1,compactHaloInfo(1)
            compactHaloInfo(idx) = dimSizes(i)
            idx = idx + 1
        end do

        !
        ! 7-8: Add MPI info
        !
        idx = 7
        compactHaloInfo(idx) = domain % dminfo % comm
        idx = idx + 1
        compactHaloInfo(idx) = domain % dminfo % my_proc_id
        idx = idx + 1

        !
        ! 9: Add number of halo layers
        !
        compactHaloInfo(idx) = activeHaloLayers
        idx = idx + 1

        !
        ! 10: Add number of send endpoints
        !
        compactHaloInfo(idx) = nSendEndpoints
        idx = idx + 1

        !
        ! 11: foreach (endpoint) { endPointID foreach (halolayer) {nList srcList(1:nList) destList(1:nList)} }
        !
        do iendpt=1,nSendEndpoints
            compactHaloInfo(idx) = sendEndpoints(iendpt)
            idx = idx + 1
            ioffset = 0
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    sendListCursor => sendList
                    do while (associated(sendListCursor))
                        found = .false.
                        exchListPtr => sendListCursor % halos(i) % exchList
                        do while (associated(exchListPtr))
                            if (exchListPtr % endPointID == sendEndpoints(iendpt)) then
                                found = .true.
                                compactHaloInfo(idx) = exchListPtr % nList
                                idx = idx + 1
                                compactHaloInfo(idx:idx+exchListPtr % nList-1) = exchListPtr % srcList(1:exchListPtr % nList)
                                idx = idx + exchListPtr % nList
                                compactHaloInfo(idx:idx+exchListPtr % nList-1) = exchListPtr % destList(1:exchListPtr % nList) + ioffset
                                idx = idx + exchListPtr % nList
                                ioffset = ioffset + exchListPtr % nList
                                exit
                            end if
                            exchListPtr => exchListPtr % next
                        end do
                        if (.not. found) then
                            compactHaloInfo(idx) = 0
                            idx = idx + 1
                        end if
                        sendListCursor => sendListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
                    end do
                end if
            end do
        end do

        deallocate(sendEndpoints)

        !
        ! 12: Add number of receive endpoints
        !
        compactHaloInfo(idx) = nRecvEndpoints
        idx = idx + 1

        !
        ! 13: foreach (endpoint) { endPointID foreach (halolayer) {nList srcList(1:nList) destList(1:nList)} }
        !
        do iendpt=1,nRecvEndpoints
            compactHaloInfo(idx) = recvEndpoints(iendpt)
            idx = idx + 1
            ioffset = 0
            do i=1,nHaloLayers
                if (useHalo(i)) then
                    recvListCursor => recvList
                    do while (associated(recvListCursor))
                        found = .false.
                        exchListPtr => recvListCursor % halos(i) % exchList
                        do while (associated(exchListPtr))
                            if (exchListPtr % endPointID == recvEndpoints(iendpt)) then
                                found = .true.
                                compactHaloInfo(idx) = exchListPtr % nList
                                idx = idx + 1
                                compactHaloInfo(idx:idx+exchListPtr % nList-1) = exchListPtr % srcList(1:exchListPtr % nList) + ioffset
                                idx = idx + exchListPtr % nList
                                compactHaloInfo(idx:idx+exchListPtr % nList-1) = exchListPtr % destList(1:exchListPtr % nList)
                                idx = idx + exchListPtr % nList
                                ioffset = ioffset + exchListPtr % nList
                                exit
                            end if
                            exchListPtr => exchListPtr % next
                        end do
                        if (.not. found) then
                            compactHaloInfo(idx) = 0
                            idx = idx + 1
                        end if
                        recvListCursor => recvListCursor % next     ! Expected to iterate just once for MPAS-Atmosphere
                    end do
                end if
            end do
        end do

        deallocate(recvEndpoints)

        deallocate(useHalo)

    end subroutine mpas_halo_compact_halo_info


    !-----------------------------------------------------------------------
    !  routine mpas_halo_aggregate_group_info
    !
    !> \brief Aggregate exchange info from all fields in a halo exchange group
    !> \author Michael Duda
    !> \date   23 November 2021
    !> \details
    !>  Given an mpas_halo_group, this routine aggregates information from across
    !>  all mpas_halo_field members of the group.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_halo_aggregate_group_info(group, ierr)

        use mpas_derived_types, only : mpas_halo_group

        type (mpas_halo_group), intent(inout) :: group
        integer, intent(out), optional :: ierr

        ! Local variables
        integer :: i, j, idx, ihalo, iendp, nlist, idxsend, idxrecv
        integer :: ndims, ninnerelems
        integer :: maxGroupSendNeighbors, maxGroupRecvNeighbors
        integer, allocatable, dimension(:) :: sendNeighbors, recvNeighbors
        integer, allocatable, dimension(:,:) :: sendCounts, recvCounts


        if (present(ierr)) then
            ierr = 0
        end if

        !
        ! Compute an upper bound on the number of send and recv neighbors for this group
        !
        maxGroupSendNeighbors = 0
        maxGroupRecvNeighbors = 0
        do i = 1, group % nFields
            maxGroupSendNeighbors = maxGroupSendNeighbors + group % fields(i) % compactHaloInfo(10)

            ! Skip past send info to get to recv neighbor count
            idx = 11
            do iendp = 1, group % fields(i) % compactHaloInfo(10)
                idx = idx + 1                ! skip over endPointID
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

            maxGroupRecvNeighbors = maxGroupRecvNeighbors + group % fields(i) % compactHaloInfo(idx)
        end do


        !
        ! Create a list of unique send and recv neighbors for this group
        !
        allocate(sendNeighbors(maxGroupSendNeighbors))
        allocate(recvNeighbors(maxGroupRecvNeighbors))

        sendNeighbors(:) = -1
        recvNeighbors(:) = -1

        group % nGroupSendNeighbors = 0
        group % nGroupRecvNeighbors = 0

        do i = 1, group % nFields
            idxsend = 10
            idx = idxsend + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxsend)
               
                ! Try to locate this endPointID in the list
                do j = 1, group % nGroupSendNeighbors
                    if (sendNeighbors(j) == group % fields(i) % compactHaloInfo(idx)) then
                        exit
                    end if
                end do

                ! If endPointID was not found, add it to the list
                if (j > group % nGroupSendNeighbors) then
                    group % nGroupSendNeighbors = group % nGroupSendNeighbors + 1
                    sendNeighbors(group % nGroupSendNeighbors) = group % fields(i) % compactHaloInfo(idx)
                end if

                ! Skip over remaining info for this endpoint
                idx = idx + 1                ! skip over endPointID
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

            idxrecv = idx
            idx = idxrecv + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxrecv)
               
                ! Try to locate this endPointID in the list
                do j = 1, group % nGroupRecvNeighbors
                    if (recvNeighbors(j) == group % fields(i) % compactHaloInfo(idx)) then
                        exit
                    end if
                end do

                ! If endPointID was not found, add it to the list
                if (j > group % nGroupRecvNeighbors) then
                    group % nGroupRecvNeighbors = group % nGroupRecvNeighbors + 1
                    recvNeighbors(group % nGroupRecvNeighbors) = group % fields(i) % compactHaloInfo(idx)
                end if

                ! Skip over remaining info for this endpoint
                idx = idx + 1                ! skip over endPointID
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

        end do

        allocate(group % groupSendNeighbors(group % nGroupSendNeighbors))
        group % groupSendNeighbors(:) = sendNeighbors(1:group % nGroupSendNeighbors)

        allocate(group % groupRecvNeighbors(group % nGroupRecvNeighbors))
        group % groupRecvNeighbors(:) = recvNeighbors(1:group % nGroupRecvNeighbors)

        deallocate(sendNeighbors)
        deallocate(recvNeighbors)


        !
        ! Compute total size of send and receive buffers for this group
        !
        group % groupSendBufSize = 0
        group % groupRecvBufSize = 0

        do i = 1, group % nFields

            ndims = group % fields(i) % compactHaloInfo(1)
            ninnerelems = 1
            do j = 1, ndims - 1    ! Do not include right-most dimension (nCells, nEdges, or nVertices)
                ninnerelems = ninnerelems * group % fields(i) % compactHaloInfo(j + 1)    ! First dim is at (2), etc.
            end do

            idxsend = 10
            idx = idxsend + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxsend)
               
                idx = idx + 1                ! skip over endPointID
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    group % groupSendBufSize = group % groupSendBufSize + ninnerelems * nlist

                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

            idxrecv = idx
            idx = idxrecv + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxrecv)
               
                idx = idx + 1                ! skip over endPointID
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    group % groupRecvBufSize = group % groupRecvBufSize + ninnerelems * nlist

                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

        end do


        !
        ! Compute sizes and offsets in group send and recv buffers for all fields in group
        !
        allocate(group % groupPackOffsets(group % nGroupSendNeighbors, group % nFields))
        allocate(group % groupSendOffsets(group % nGroupSendNeighbors))
        allocate(group % groupSendCounts(group % nGroupSendNeighbors))
        group % groupPackOffsets(:,:) = -1
        group % groupSendOffsets(:) = -1
        group % groupSendCounts(:) = -1

        allocate(group % groupUnpackOffsets(group % nGroupRecvNeighbors, group % nFields))
        allocate(group % groupRecvOffsets(group % nGroupRecvNeighbors))
        allocate(group % groupRecvCounts(group % nGroupRecvNeighbors))
        group % groupUnpackOffsets(:,:) = -1
        group % groupRecvOffsets(:) = -1
        group % groupRecvCounts(:) = -1

        allocate(sendCounts(group % nGroupSendNeighbors, group % nFields))
        allocate(recvCounts(group % nGroupRecvNeighbors, group % nFields))
        sendCounts(:,:) = 0
        recvCounts(:,:) = 0

        do i = 1, group % nFields

            ndims = group % fields(i) % compactHaloInfo(1)
            ninnerelems = 1
            do j = 1, ndims - 1    ! Do not include right-most dimension (nCells, nEdges, or nVertices)
                ninnerelems = ninnerelems * group % fields(i) % compactHaloInfo(j + 1)    ! First dim is at (2), etc.
            end do

            idxsend = 10
            idx = idxsend + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxsend)

                ! Find neighbor in groupSendNeighbors
                do j = 1, group % nGroupSendNeighbors
                    if (group % fields(i) % compactHaloInfo(idx) == group % groupSendNeighbors(j)) then
                        exit
                    end if
                end do
               
                idx = idx + 1                ! skip over endPointID
                sendCounts(j, i) = 0
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    sendCounts(j, i) = sendCounts(j, i) + nlist * ninnerelems

                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

            idxrecv = idx
            idx = idxrecv + 1
            do iendp = 1, group % fields(i) % compactHaloInfo(idxrecv)

                ! Find neighbor in groupRecvNeighbors
                do j = 1, group % nGroupRecvNeighbors
                    if (group % fields(i) % compactHaloInfo(idx) == group % groupRecvNeighbors(j)) then
                        exit
                    end if
                end do
               
                idx = idx + 1                ! skip over endPointID
                recvCounts(j, i) = 0
                do ihalo = 1, group % fields(i) % compactHaloInfo(9)
                    nlist = group % fields(i) % compactHaloInfo(idx)
                    recvCounts(j, i) = recvCounts(j, i) + nlist * ninnerelems

                    idx = idx + 1            ! skip over nList
                    idx = idx + 2 * nlist    ! skip over srcList and destList
                end do
            end do

        end do

        do j = 1, group % nGroupSendNeighbors
            group % groupPackOffsets(j, 1) = 0
            do i = 2, group % nFields
                group % groupPackOffsets(j, i) = group % groupPackOffsets(j, i-1) + sendCounts(j, i-1)
            end do
        end do

        do j = 1, group % nGroupRecvNeighbors
            group % groupUnpackOffsets(j, 1) = 0
            do i = 2, group % nFields
                group % groupUnpackOffsets(j, i) = group % groupUnpackOffsets(j, i-1) + recvCounts(j, i-1)
            end do
        end do

        do j = 1, group % nGroupSendNeighbors
            group % groupSendCounts(j) = 0
            do i = 1, group % nFields
                group % groupSendCounts(j) = group % groupSendCounts(j) + sendCounts(j, i)
            end do

            if (j == 1) then
                group % groupSendOffsets(j) = 1
            else
                group % groupSendOffsets(j) = group % groupSendOffsets(j-1) + group % groupSendCounts(j-1)

                do i = 1, group % nFields
                    group % groupPackOffsets(j, i) = group % groupPackOffsets(j, i) + group % groupSendOffsets(j) - 1
                end do
            end if
        end do

        do j = 1, group % nGroupRecvNeighbors
            group % groupRecvCounts(j) = 0
            do i = 1, group % nFields
                group % groupRecvCounts(j) = group % groupRecvCounts(j) + recvCounts(j, i)
            end do

            if (j == 1) then
                group % groupRecvOffsets(j) = 1
            else
                group % groupRecvOffsets(j) = group % groupRecvOffsets(j-1) + group % groupRecvCounts(j-1)

                do i = 1, group % nFields
                    group % groupUnpackOffsets(j, i) = group % groupUnpackOffsets(j, i) + group % groupRecvOffsets(j) - 1
                end do
            end if
        end do

        deallocate(sendCounts)
        deallocate(recvCounts)

    end subroutine mpas_halo_aggregate_group_info

end module mpas_halo
