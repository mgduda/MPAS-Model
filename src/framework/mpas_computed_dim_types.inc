   integer, parameter :: MPAS_COMPUTED_DIM_NOERR = 0, &
                         MPAS_COMPUTED_DIM_ERROR = -1

   abstract interface
      function mpas_computed_dim_size_fn(block) result(dim_size)
         import block_type
         import mpas_streamManager_type

         type (block_type), intent(in) :: block
         integer :: dim_size
      end function mpas_computed_dim_size_fn
   end interface

   type mpas_dim_size_fn_list
      integer :: nameLen = 0
      character (len=StrKIND) :: funcName = ''
      procedure (mpas_computed_dim_size_fn), pointer, nopass :: dimSizeFunc => null()
      type (mpas_dim_size_fn_list), pointer :: next => null()
   end type mpas_dim_size_fn_list
