    integer, parameter :: MPAS_HALO_INVALID   = -1

    integer, parameter :: MPAS_HALO_REAL      = 5001, &
                          MPAS_HALO_INTEGER   = 5002

    type mpas_halo_field
        integer :: nDims = MPAS_HALO_INVALID
        integer :: fieldType = MPAS_HALO_INVALID
        type (field2DReal), pointer :: r2d => null()    ! Valid when nDims == 2 and fieldType == MPAS_HALO_REAL
        type (field3DReal), pointer :: r3d => null()    ! Valid when nDims == 3 and fieldType == MPAS_HALO_REAL
    end type mpas_halo_field

    type mpas_halo_group
        character(len=StrKIND) :: groupName = ''
        integer :: nFields = MPAS_HALO_INVALID
        type (mpas_halo_field), dimension(:), pointer :: fields => null()
        type (mpas_halo_group), pointer :: next => null()
    end type mpas_halo_group
