   type mpas_bookkeeper
      type (c_ptr) :: c_bookkeeper = c_null_ptr
      procedure (bookkeeping_init_function), pass, pointer :: init => null()
      procedure (bookkeeping_finalize_function), pass, pointer :: finalize => null()
      procedure (bookkeeping_requisite_function), pass, pointer :: declare_requisite => null()
   end type mpas_bookkeeper

   abstract interface
      function bookkeeping_init_function(this) result(ierr)
         import mpas_bookkeeper
         class (mpas_bookkeeper) :: this
         integer :: ierr
      end function bookkeeping_init_function

      function bookkeeping_finalize_function(this) result(ierr)
         import mpas_bookkeeper
         class (mpas_bookkeeper) :: this
         integer :: ierr
      end function bookkeeping_finalize_function

      function bookkeeping_requisite_function(this, fieldname) result(ierr)
         import mpas_bookkeeper
         class (mpas_bookkeeper) :: this
         character(len=*), intent(in) :: fieldname
         integer :: ierr
      end function bookkeeping_requisite_function
   end interface
