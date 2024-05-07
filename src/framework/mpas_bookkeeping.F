! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_bookkeeping
!
!> \brief Manage copies of Registry-defined fields across host and device
!> \author Michael Duda
!> \date   3 May 2024
!> \details
!>  This module provides ...
!
!-----------------------------------------------------------------------
module mpas_bookkeeping

    public :: mpas_bookkeeping_new_bookkeeper


contains


    !-----------------------------------------------------------------------
    !  routine mpas_bookkeeping_new_bookkeeper
    !
    !> \brief Returns a pointer to a new mpas_bookkeeper instance
    !> \author Michael Duda
    !> \date   3 May 2024
    !> \details
    !>  This routine returns a pointer to a newly allocated instance of an
    !>  mpas_bookkeeper. The new instance has valid methods X(), Y(),
    !>  and Z() that may be called.
    !>
    !>  After all queries via the MPAS_streamInfo_type instance have been
    !>  completed, the instance finalize() method should be called before the
    !>  instance is deallocated.
    !
    !-----------------------------------------------------------------------
    function mpas_bookkeeping_new_bookkeeper() result(new_bookkeeper)

        use mpas_derived_types, only : mpas_bookkeeper

        implicit none

        ! Return value
        type (mpas_bookkeeper), pointer :: new_bookkeeper

        allocate(new_bookkeeper)
        new_bookkeeper % init => bookkeeper_init
        new_bookkeeper % finalize => bookkeeper_finalize
        new_bookkeeper % declare_requisite => bookkeeper_requisite

    end function mpas_bookkeeping_new_bookkeeper


    !-----------------------------------------------------------------------
    !  routine bookkeeper_init
    !
    !> \brief Initializes an MPAS_streamInfo_type instance from a streams XML file
    !> \author Michael Duda
    !> \date   15 November 2023
    !> \details
    !>  This routine should be called as a method within an MPAS_streamInfo_type
    !>  instance, e.g., streaminfo % init(...). Given the name of an MPAS streams
    !>  XML file, this method initializes the instance so that later queries may
    !>  be made with the query() method.
    !
    !-----------------------------------------------------------------------
    function bookkeeper_init(this) result(ierr)

        use mpas_derived_types, only : mpas_bookkeeper

        implicit none

        ! Arguments
        class (mpas_bookkeeper) :: this

        ! Return value
        integer :: ierr

        ! Local variables
!        character(kind=c_char), dimension(len(stream_filename)+1) :: c_stream_filename

        interface
            function bookkeeping_init(bookkeeper) bind(C, name='bookkeeping_init') result(ierr)
                use iso_c_binding, only : c_int, c_ptr
                type(c_ptr) :: bookkeeper
                integer(kind=c_int) :: ierr
            end function bookkeeping_init
        end interface


        ierr = bookkeeping_init(this % c_bookkeeper)

    end function bookkeeper_init


    !-----------------------------------------------------------------------
    !  routine bookkeeper_finalize
    !
    !> \brief Finalizes an instance of the MPAS_streamInfo_type type
    !> \author Michael Duda
    !> \date   15 November 2023
    !> \details
    !>  This routine finalizes an instance of the MPAS_streamInfo_type type
    !>  after all queries about the contents of the streams XML file associated
    !>  with the instance have been completed. This routine should be called as
    !>  a method within an MPAS_streamInfo_type type, e.g.,
    !>  streaminfo % finalize().
    !
    !-----------------------------------------------------------------------
    function bookkeeper_finalize(this) result(ierr)

        use mpas_derived_types, only : mpas_bookkeeper

        implicit none

        ! Arguments
        class (mpas_bookkeeper) :: this

        ! Return value
        integer :: ierr

        interface
            function bookkeeping_finalize(bookkeeper) bind(C, name='bookkeeping_finalize') result(ierr)
                use iso_c_binding, only : c_int, c_ptr
                type(c_ptr) :: bookkeeper
                integer(kind=c_int) :: ierr
            end function bookkeeping_finalize
        end interface


        ierr = bookkeeping_finalize(this % c_bookkeeper)

    end function bookkeeper_finalize


    !-----------------------------------------------------------------------
    !  routine bookkeeper_requisite
    !
    !> \brief ...
    !> \author Michael Duda
    !> \date   15 November 2023
    !> \details
    !>  This routine ...
    !
    !-----------------------------------------------------------------------
    function bookkeeper_requisite(this, fieldname) result(ierr)

        use mpas_derived_types, only : mpas_bookkeeper
        use mpas_c_interfacing, only : mpas_f_to_c_string
        use iso_c_binding, only : c_char

        implicit none

        ! Arguments
        class (mpas_bookkeeper) :: this
        character(len=*), intent(in) :: fieldname

        ! Local variables
        character(kind=c_char), dimension(len(fieldname)+1) :: c_fieldname

        ! Return value
        integer :: ierr

        interface
            function bookkeeping_declare_requisite(bookkeeper, fieldname) bind(C, name='bookkeeping_declare_requisite') result(ierr)
                use iso_c_binding, only : c_int, c_char, c_ptr
                type(c_ptr), value :: bookkeeper
                character(kind=c_char), dimension(*), intent(in) :: fieldname
                integer(kind=c_int) :: ierr
            end function bookkeeping_declare_requisite
        end interface


        call mpas_f_to_c_string(fieldname, c_fieldname)
        ierr = bookkeeping_declare_requisite(this % c_bookkeeper, c_fieldname)

    end function bookkeeper_requisite

end module mpas_bookkeeping
