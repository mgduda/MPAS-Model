   integer, parameter :: MPAS_COMPUTED_DIM_NOERR = 0, &
                         MPAS_COMPUTED_DIM_ERROR = -1

   abstract interface
      function mpas_computed_dim_size_function(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)
         import block_type
         import mpas_streamManager_type

         type (block_type), intent(in) :: block
         type (mpas_streamManager_type), intent(inout) :: manager
         integer, intent(in) :: globalDimSize
         integer, intent(in) :: numBlocks
         integer, dimension(:), pointer :: ownedIndices
         integer :: iErr
      end function mpas_computed_dim_size_function
   end interface

   type mpas_dim_size_fn_list
      integer :: nameLen
      character (len=StrKIND) :: decompName
      procedure (mpas_computed_dim_size_function), pointer, nopass :: dimSizeFunc => null()
      type (mpas_dim_size_fn_list), pointer :: next => null()
   end type mpas_dim_size_fn_list
