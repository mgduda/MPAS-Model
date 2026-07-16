! Copyright (c) 2026, the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_computed_dims
!
!> \brief MPAS computed dimensions module
!> \author Michael Duda
!> \date   16 July 2026
!> \details
!>  Work in progress...
!
!-----------------------------------------------------------------------

#define COMMA ,
#define DEBUG_WRITE(M) ! call mpas_log_write(M)
#define WARN_WRITE(M) call mpas_log_write(M , messageType=MPAS_LOG_WARN)
#define ERROR_WRITE(M) call mpas_log_write(M , messageType=MPAS_LOG_CRIT)

module mpas_computed_dims

   use mpas_kind_types
   use mpas_derived_types
   use mpas_stream_manager
   use mpas_log

   implicit none


contains


   function three_hundred(block) result(dim_size)

      type (block_type), intent(in) :: block

      integer :: dim_size

      dim_size = 300

   end function three_hundred


!-----------------------------------------------------------------------
!  routine mpas_computed_dims_create_list
!
!> \brief MPAS method list create routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine creates a decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_computed_dims_create_list(methodList)
      type (mpas_dim_size_fn_list), pointer :: methodList

      DEBUG_WRITE('Creating computed dimensions list')

      if ( .not. associated(methodList) ) then
         allocate(methodList)
      end if

   end subroutine mpas_computed_dims_create_list


!-----------------------------------------------------------------------
!  routine mpas_computed_dims_destroy_list
!
!> \brief MPAS decomp list destruction routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine destroys a decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_computed_dims_destroy_list(methodList)
      type (mpas_dim_size_fn_list), pointer :: methodList

      type (mpas_dim_size_fn_list), pointer :: methodCursor

      DEBUG_WRITE('Destroying decomposition list')

      do while (associated(methodList % next))
         methodCursor => methodList % next

         if ( associated(methodCursor % next) ) then
            methodList % next => methodCursor % next
         else
            nullify(methodList % next)
         end if

         deallocate(methodCursor)
      end do

      deallocate(methodList)
   end subroutine mpas_computed_dims_destroy_list


!-----------------------------------------------------------------------
!  routine mpas_decomp_register_method
!
!> \brief MPAS decomp list register routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine registers a decomposition function in a decomposition list.
!>  If present, the optional argument iErr is set to MPAS_COMPUTED_DIM_ERROR if 
!>  the function is already registered and to MPAS_COMPUTED_DIM_NOERR otherwise.
!
!-----------------------------------------------------------------------
   subroutine mpas_computed_dims_register_method(methodList, methodName, methodFunc, iErr)
      type (mpas_dim_size_fn_list), pointer :: methodList
      character (len=*), intent(in) :: methodName
      procedure (mpas_computed_dim_size_fn), pointer :: methodFunc
      integer, intent(out) :: iErr

      type (mpas_dim_size_fn_list), pointer :: methodCursor
      type (mpas_dim_size_fn_list), pointer :: methodCursorPrev
      integer :: nameLen

      iErr = MPAS_COMPUTED_DIM_NOERR

      nameLen = len_trim(methodName)

      methodCursor => methodList
      nullify(methodCursorPrev)
      do while (associated(methodCursor))
         if ( nameLen == methodCursor % nameLen ) then
            if ( trim(methodName) == trim(methodCursor % funcName) ) then
               DEBUG_WRITE('Decomposition ' // trim(methodName) // ' is already registered.')
               iErr = MPAS_COMPUTED_DIM_ERROR
               return
            end if
         end if
         methodCursorPrev => methodCursor
         methodCursor => methodCursor % next
      end do

      methodCursor => methodCursorPrev

      DEBUG_WRITE('* Registering decomposition '//trim(methodName))

      allocate(methodCursor % next)
      methodCursor => methodCursor % next
      methodCursor % funcName = methodName
      methodCursor % dimSizeFunc => methodFunc
      methodCursor % nameLen = nameLen
      nullify(methodCursor % next)
   end subroutine mpas_computed_dims_register_method


!-----------------------------------------------------------------------
!  routine mpas_decomp_get_method
!
!> \brief MPAS decomp list query routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine querys a decomposition list for a specific decomposition function.
!>  If present, the optional argument iErr is set to MPAS_COMPUTED_DIM_ERROR if 
!>  the requested function is not found and to MPAS_COMPUTED_DIM_NOERR otherwise.
!
!-----------------------------------------------------------------------
   function mpas_computed_dims_get_method(methodFuncList, methodName) result(method)

      ! Arguments
      type (mpas_dim_size_fn_list), pointer :: methodFuncList
      character(len=*), intent(in) :: methodName

      ! Return value
      procedure(mpas_computed_dim_size_fn), pointer :: method

      ! Local variables
      type (mpas_dim_size_fn_list), pointer :: methodCursor
      integer :: nameLen

!      method => three_hundred

!      iErr = MPAS_COMPUTED_DIM_NOERR

      nameLen = len_trim(methodName)

      methodCursor => methodFuncList
      do while (associated(methodCursor))
         if ( nameLen == methodCursor % nameLen ) then
            if ( trim(methodCursor % funcName) == trim(methodName) ) then
               DEBUG_WRITE('Returning method '//trim(methodName))
               method => methodCursor % dimSizeFunc
               return
            end if
         end if
         methodCursor => methodCursor % next
      end do

!      iErr = MPAS_COMPUTED_DIM_ERROR
!      DEBUG_WRITE('Decomposition '//trim(methodName)//' not found in decomposition list.')

   end function mpas_computed_dims_get_method


!-----------------------------------------------------------------------
!  routine mpas_decomp_remove_method
!
!> \brief MPAS decomp list removal routine
!> \author Doug Jacobsen
!> \date   03/11/2015
!> \details
!>  This routine removes a specific decomposition function from a  decomposition list.
!
!-----------------------------------------------------------------------
   subroutine mpas_computed_dims_remove_method(methodList, methodName, iErr)
      type (mpas_dim_size_fn_list), pointer :: methodList
      character (len=*), intent(in) :: methodName
      integer, intent(out) :: iErr

      integer :: nameLen
      type (mpas_dim_size_fn_list), pointer :: methodCursor, methodCursorPrev

      iErr = MPAS_COMPUTED_DIM_NOERR

      nameLen = len_trim(methodName)

      methodCursor => methodList
      nullify(methodCursorPrev)
      do while (associated(methodCursor))
         if ( nameLen == methodCursor % nameLen ) then
            if ( trim(methodCursor % funcName) == trim(methodName) ) then
               DEBUG_WRITE('* Removing methodosition '//trim(methodName))

               if ( associated(methodCursor % next) ) then
                  methodCursorPrev % next => methodCursor % next
               else
                  nullify(methodCursorPrev % next)
               end if

               deallocate(methodCursor)
               return
            end if
         end if
         methodCursorPrev => methodCursor
         methodCursor => methodCursor % next
      end do
   end subroutine mpas_computed_dims_remove_method

end module mpas_computed_dims
