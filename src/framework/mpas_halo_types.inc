    integer, parameter :: MPAS_HALO_INVALID   = -1

    integer, parameter :: MPAS_HALO_REAL      = 5001, &
                          MPAS_HALO_INTEGER   = 5002

    type mpas_halo_field
        character(len=StrKIND) :: fieldName = ''
        integer :: nDims = MPAS_HALO_INVALID
        integer :: fieldType = MPAS_HALO_INVALID
        integer :: timeLevel = MPAS_HALO_INVALID
        integer, dimension(:), pointer :: compactHaloInfo => null()
        real (kind=RKIND), dimension(:,:), pointer :: r2arr => null()
        real (kind=RKIND), dimension(:,:,:), pointer :: r3arr => null()
    end type mpas_halo_field

    type mpas_halo_group
        character(len=StrKIND) :: groupName = ''
        integer :: nFields = MPAS_HALO_INVALID
        type (mpas_halo_field), dimension(:), pointer :: fields => null()

        integer :: nGroupSendNeighbors = MPAS_HALO_INVALID
        integer :: groupSendBufSize = MPAS_HALO_INVALID     ! Total number of elements to be sent in a group exchange
        real (kind=RKIND), dimension(:), pointer :: sendBuf => null()
        integer, dimension(:), pointer :: sendRequests => null()
        integer, dimension(:,:), pointer :: groupPackOffsets => null()
        integer, dimension(:), pointer :: groupSendNeighbors => null()
        integer, dimension(:), pointer :: groupSendOffsets => null()
        integer, dimension(:), pointer :: groupSendCounts => null()

        integer :: nGroupRecvNeighbors = MPAS_HALO_INVALID
        integer :: groupRecvBufSize = MPAS_HALO_INVALID     ! Total number of elements to be received in a group exchange
        real (kind=RKIND), dimension(:), pointer :: recvBuf => null()
        integer, dimension(:), pointer :: recvRequests => null()
        integer, dimension(:,:), pointer :: groupUnpackOffsets => null()
        integer, dimension(:), pointer :: groupRecvNeighbors => null()
        integer, dimension(:), pointer :: groupRecvOffsets => null()
        integer, dimension(:), pointer :: groupRecvCounts => null()

        type (mpas_halo_group), pointer :: next => null()
    end type mpas_halo_group
