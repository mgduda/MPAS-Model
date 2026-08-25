   function num_boundary_cells(block, manager) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager

      integer :: dim_size

      dim_size = num_boundary_elems(block, manager, 'Cell', 'Cells')

   end function num_boundary_cells


   function num_boundary_edges(block, manager) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager

      integer :: dim_size

      dim_size = num_boundary_elems(block, manager, 'Edge', 'Edges')

   end function num_boundary_edges


   function num_boundary_elems(block, manager, elemType, elemTypePlural) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural

      integer :: dim_size

      integer, pointer :: nElemsSolve

      integer, dimension(:), pointer :: indexToElemID
      integer, dimension(:), pointer :: bdyMaskElem

      integer :: iElem
      integer :: nBdyElems_local, nBdyElems


      call read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

!      call mpas_log_write('Min/max bdyMaskElem=$i $i', intArgs=[minval(bdyMaskElem), maxval(bdyMaskElem)])

      nBdyElems_local = 0
      do iElem = 1, nElemsSolve
          if (bdyMaskElem(iElem) > 0) nBdyElems_local = nBdyElems_local + 1
      end do

      call mpas_dmpar_sum_int(block % domain % dminfo, nBdyElems_local, nBdyElems)

!      call mpas_log_write('nBdyElems = $i', intArgs=[nBdyElems])

      deallocate(bdyMaskElem)

      dim_size = nBdyElems

   end function num_boundary_elems


   function boundary_cells(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_timer

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr
      
      call mpas_timer_start('boundary_cells')
      iErr = boundary_elems(block, manager, ownedIndices, 'Cell', 'Cells')
      call mpas_timer_stop('boundary_cells')

   end function boundary_cells


   function boundary_edges(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_timer

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr
      
      call mpas_timer_start('boundary_edges')
      iErr = boundary_elems(block, manager, ownedIndices, 'Edge', 'Edges')
      call mpas_timer_stop('boundary_edges')

   end function boundary_edges


   function boundary_elems(block, manager, ownedIndices, elemType, elemTypePlural) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams
#ifdef MPAS_USE_MPI_F08
      use mpi_f08
#else
      use mpi
#endif
      use mpas_sort

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, dimension(:), pointer :: ownedIndices
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural

      integer :: iErr
      
      integer, pointer :: nElemsSolve

      integer, dimension(:), pointer :: indexToElemID
      integer, dimension(:), pointer :: bdyMaskElem

      integer :: iElem
      integer :: nBdyElems

#ifdef MPAS_USE_MPI_F08
      type (MPI_Comm) :: comm
      type (MPI_Request) :: req
#else
      integer :: comm
      integer :: req
#endif
      integer :: nRanks, myRank

      integer, dimension(:,:), allocatable :: temp

      integer :: i, j, k
      integer :: n, n_next
      integer, dimension(:), pointer :: bdyIDs, bdyIDs_next
      integer :: srcRank, dstRank


      iErr = 0
      
      call read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

      nBdyElems = 0
      do iElem = 1, nElemsSolve
          if (bdyMaskElem(iElem) > 0) nBdyElems = nBdyElems + 1
      end do

!      call mpas_dmpar_sum_int(block % domain % dminfo, nBdyElems_local, nBdyElems)
!      if (nBdyElems > 0) then
!          call mpas_log_write('Again, nBdyElems = $i', intArgs=[nBdyElems])
!      else
!          call mpas_log_write('No boundary Elems')
!      end if

      ! comm = block % domain % dminfo % comm
      if (nBdyElems > 0) then
          call MPI_Comm_split(block % domain % dminfo % comm, 1, block % domain % dminfo % my_proc_id, comm, ierr)
      else
          call MPI_Comm_split(block % domain % dminfo % comm, 0, block % domain % dminfo % my_proc_id, comm, ierr)
      end if
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to split comm', messageType=MPAS_LOG_ERR)
      end if

      if (nBdyElems > 0) then
          call MPI_Comm_size(comm, nRanks, ierr)
          if (ierr /= MPI_SUCCESS) then
              call mpas_log_write('Failed to get comm size', messageType=MPAS_LOG_ERR)
          end if

          call MPI_Comm_rank(comm, myRank, ierr)
          if (ierr /= MPI_SUCCESS) then
              call mpas_log_write('Failed to get comm rank', messageType=MPAS_LOG_ERR)
          end if

          allocate(temp(2,nBdyElems))
          i = 1
          do iElem = 1, nElemsSolve
              if (bdyMaskElem(iElem) > 0) then
                  temp(1,i) = indexToElemID(iElem)   ! Parent element
                  temp(2,i) = i                      ! Compute order on task
                  i = i + 1
              end if
          end do

          call mpas_quicksort(nBdyElems, temp)

          allocate(ownedIndices(nBdyElems))
          do iElem = 1, nBdyElems
              ownedIndices(iElem) = iElem            ! Tentative boundary-only offset
          end do

          srcRank = mod(myRank + nRanks - 1, nRanks)
          dstRank = mod(myRank + 1, nRanks)
          n = nBdyElems
          allocate(bdyIDs(n))

          bdyIDs(:) = temp(1,:)     ! Full domain boundary offsets

          do i = 1, nRanks - 1
              call MPI_Irecv(n_next, 1, MPI_INTEGER, srcRank, srcRank, comm, req, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Recv from $i failed', intArgs=[srcRank], messageType=MPAS_LOG_ERR)

              call MPI_Send(n, 1, MPI_INTEGER, dstRank, myRank, comm, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Send to $i failed', intArgs=[dstRank], messageType=MPAS_LOG_ERR)

              call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Request failed', messageType=MPAS_LOG_ERR)

!              call mpas_log_write('Now handling $i boundary elements...', intArgs=[n_next])

              allocate(bdyIDs_next(n_next))

              call MPI_Irecv(bdyIDs_next, n_next, MPI_INTEGER, srcRank, srcRank, comm, req, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Recv from $i failed', intArgs=[srcRank], messageType=MPAS_LOG_ERR)

              call MPI_Send(bdyIDs, n, MPI_INTEGER, dstRank, myRank, comm, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Send to $i failed', intArgs=[dstRank], messageType=MPAS_LOG_ERR)

              call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
              if (ierr /= MPI_SUCCESS) call mpas_log_write('Request failed', messageType=MPAS_LOG_ERR)

              deallocate(bdyIDs)
              bdyIDs => bdyIDs_next
              n = n_next

              j = 0
              do k = 1, nBdyElems
                  do while (j < n)
                      if (bdyIDs(j+1) > temp(1,k)) exit
                      j = j + 1
                  end do
                  ownedIndices(k) = ownedIndices(k) + j
              end do

          end do
          deallocate(bdyIDs)

          ! Place elements back in original order
          temp(1,:) = ownedIndices(:)                ! Tentative boundary-only offset
          do iElem = 1, nBdyElems
              ownedIndices(temp(2,iElem)) = temp(1,iElem)
          end do

          deallocate(temp)

      else

          allocate(ownedIndices(nBdyElems))

      end if
      deallocate(bdyMaskElem)

      call MPI_Comm_free(comm, ierr)
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to free comm', messageType=MPAS_LOG_ERR)
      end if

   end function boundary_elems


   subroutine read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural
      integer, dimension(:), pointer :: indexToElemID  ! Output
      integer, dimension(:), pointer :: bdyMaskElem  ! Output
      integer, pointer :: nElemsSolve  ! Output

      integer :: ierr
      character(len=StrKIND) :: mesh_stream, mesh_filename

      integer :: fieldtype

      type (MPAS_IO_Handle_type) :: handle


      ierr = block % domain % core % get_mesh_stream(block % configs, block % domain % streamInfo, mesh_stream)
!      call mpas_log_write('Mesh stream is '//trim(mesh_stream))

      call MPAS_get_stream_filename(manager, trim(mesh_stream), filename=mesh_filename)
!      call mpas_log_write('Mesh filename is '//trim(mesh_filename))

      nullify(indexToElemID)
      call mpas_pool_get_array(block % allFields, 'indexTo'//elemType//'ID_blk', indexToElemID)

      nullify(nElemsSolve)
      call mpas_pool_get_dimension(block % dimensions, 'n'//elemTypePlural//'Solve', nElemsSolve)

      if (.not. associated(nElemsSolve)) then
          call mpas_log_write('Could not find nElemsSolve', messageType=MPAS_LOG_ERR)
      end if

      allocate(bdyMaskElem(nElemsSolve))

      if (.not. associated(indexToElemID)) then
          call mpas_log_write('Could not find indexToElemID', messageType=MPAS_LOG_ERR)
      end if

      handle = MPAS_io_open(trim(mesh_filename), MPAS_IO_READ, MPAS_IO_NETCDF, manager % ioContext, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('Error opening '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      end if

      call MPAS_io_inq_var(handle, 'bdyMask'//elemType, fieldtype=fieldtype, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('Error inquiring about bdyMaskElem', messageType=MPAS_LOG_ERR)
      end if

      call MPAS_io_set_var_indices(handle, 'bdyMask'//elemType, indexToElemID(1:nElemsSolve), ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('Error setting indices', messageType=MPAS_LOG_ERR)
      end if

      call MPAS_io_get_var(handle, 'bdyMask'//elemType, bdyMaskElem, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('Error reading bdyMaskElem', messageType=MPAS_LOG_ERR)
      end if

      call MPAS_io_close(handle)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('Error closing '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      end if

   end subroutine read_bdy_mask
