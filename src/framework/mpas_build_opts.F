! Copyright (c) 2021, The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_build_opts
!
!> \brief   Returns string describing various aspects of the MPAS build environment
!> \author  Michael G. Duda
!> \date    11 November 2021
!> \details
!>  This module defines ...
!
!-----------------------------------------------------------------------
module mpas_build_opts

    implicit none

    private

    public :: mpas_get_build_opts


contains


    !-----------------------------------------------------------------------
    !  routine mpas_get_build_opts
    !
    !> \brief Returns the multi-line build option description from top-level Makefile
    !> \author Michael G. Duda
    !> \date   11 November 2021
    !> \details
    !>  Does something or other...
    !
    !-----------------------------------------------------------------------
    function mpas_get_build_opts() result(opts)

        use iso_c_binding, only : c_f_pointer, c_ptr, c_char, c_int

        ! Return value
        character(len=:), allocatable :: opts

        ! Local variables
        type (c_ptr) :: opts_ptr
        character(kind=c_char), dimension(:), pointer :: opts_c
        integer(c_int) :: opts_len

        interface
            function get_build_opts(len) bind(c, name='get_build_opts')
                use iso_c_binding, only : c_ptr, c_int
                integer(c_int) :: len
                type (c_ptr) :: get_build_opts
            end function get_build_opts
        end interface


        opts_ptr = get_build_opts(opts_len)
        call c_f_pointer(opts_ptr, opts_c, [opts_len])

        allocate(character(len=opts_len) :: opts)

        opts(:) = transfer(opts_c(1:opts_len), opts) 

    end function mpas_get_build_opts

end module mpas_build_opts
