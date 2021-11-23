    integer, parameter :: MPAS_HALO_INVALID   = -1

    integer, parameter :: MPAS_HALO_REAL      = 5001, &
                          MPAS_HALO_INTEGER   = 5002

    type mpas_halo_field
        character(len=StrKIND) :: fieldName = ''
        integer :: nDims = MPAS_HALO_INVALID
        integer :: fieldType = MPAS_HALO_INVALID
        integer, dimension(:), pointer :: compactHaloInfo => null()
    end type mpas_halo_field

    type mpas_halo_group
        character(len=StrKIND) :: groupName = ''
        integer :: nFields = MPAS_HALO_INVALID
        type (mpas_halo_field), dimension(:), pointer :: fields => null()
        type (mpas_halo_group), pointer :: next => null()
    end type mpas_halo_group
