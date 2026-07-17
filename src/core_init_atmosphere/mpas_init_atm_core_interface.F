! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_core_interface

   use mpas_attlist
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units
   use mpas_log, only : mpas_log_write

   contains


   !***********************************************************************
   !
   !  routine init_atm_setup_core
   !
   !> \brief   Init atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use init_atm_core, only : init_atm_core_init, init_atm_core_run, init_atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => init_atm_core_init
      core % core_run => init_atm_core_run
      core % core_finalize => init_atm_core_finalize
      core % define_packages => init_atm_define_packages
      core % setup_packages => init_atm_setup_packages
      core % setup_decompositions => init_atm_setup_decompositions
      core % setup_clock => init_atm_setup_clock
      core % setup_log => init_atm_setup_log
      core % get_mesh_stream => init_atm_get_mesh_stream
      core % setup_immutable_streams => init_atm_setup_immutable_streams
      core % setup_derived_dimensions => init_atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => init_atm_setup_decomposed_dimensions
      core % setup_block => init_atm_setup_block
      core % setup_namelist => init_atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "core_variables.inc"

   end subroutine init_atm_setup_core


   !***********************************************************************
   !
   !  routine init_atm_setup_domain
   !
   !> \brief   Init atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atmosphere core.
   !
   !-----------------------------------------------------------------------
   subroutine init_atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "domain_variables.inc"

   end subroutine init_atm_setup_domain


   !***********************************************************************
   !
   !  function init_atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  init_atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_packages(configs, streamInfo, packages, iocontext) result(ierr)

      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr
      logical :: lexist

      logical, pointer :: initial_conds, sfc_update, lbcs
      logical, pointer :: gwd_stage_in, gwd_stage_out, vertical_stage_in, vertical_stage_out, met_stage_in, met_stage_out
      logical, pointer :: gwd_gsl_stage_out
      logical, pointer :: config_native_gwd_static, config_static_interp, config_vertical_grid, config_met_interp
      logical, pointer :: config_native_gwd_gsl_static
      logical, pointer :: first_guess_field
      logical, pointer :: mp_thompson_aers_in
      integer, pointer :: config_init_case
      logical, pointer :: noahmp, config_noahmp_static


      ierr = init_atm_setup_packages_when(configs, packages)
      if (ierr /= 0) then
         return
      end if

      call mpas_pool_get_config(configs, 'config_init_case', config_init_case)
      call mpas_pool_get_config(configs, 'config_static_interp', config_static_interp)
      call mpas_pool_get_config(configs, 'config_native_gwd_static', config_native_gwd_static)
      call mpas_pool_get_config(configs, 'config_native_gwd_gsl_static', config_native_gwd_gsl_static)
      call mpas_pool_get_config(configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(configs, 'config_met_interp', config_met_interp)

      nullify(initial_conds)
      call mpas_pool_get_package(packages, 'initial_condsActive', initial_conds)

      nullify(sfc_update)
      call mpas_pool_get_package(packages, 'sfc_updateActive', sfc_update)

      nullify(lbcs)
      call mpas_pool_get_package(packages, 'lbcsActive', lbcs)

      nullify(gwd_stage_in)
      call mpas_pool_get_package(packages, 'gwd_stage_inActive', gwd_stage_in)

      nullify(gwd_stage_out)
      call mpas_pool_get_package(packages, 'gwd_stage_outActive', gwd_stage_out)

      nullify(gwd_gsl_stage_out)
      call mpas_pool_get_package(packages, 'gwd_gsl_stage_outActive', gwd_gsl_stage_out)

      nullify(vertical_stage_in)
      call mpas_pool_get_package(packages, 'vertical_stage_inActive', vertical_stage_in)

      nullify(vertical_stage_out)
      call mpas_pool_get_package(packages, 'vertical_stage_outActive', vertical_stage_out)

      nullify(met_stage_in)
      call mpas_pool_get_package(packages, 'met_stage_inActive', met_stage_in)

      nullify(met_stage_out)
      call mpas_pool_get_package(packages, 'met_stage_outActive', met_stage_out)

      nullify(mp_thompson_aers_in)
      call mpas_pool_get_package(packages, 'mp_thompson_aers_inActive', mp_thompson_aers_in)

      if (.not. associated(initial_conds) .or. &
          .not. associated(sfc_update) .or. &
          .not. associated(gwd_stage_in) .or. &
          .not. associated(gwd_stage_out) .or. &
          .not. associated(gwd_gsl_stage_out) .or. &
          .not. associated(vertical_stage_in) .or. &
          .not. associated(vertical_stage_out) .or. &
          .not. associated(met_stage_in) .or. &
          .not. associated(met_stage_out) .or. &
          .not. associated(mp_thompson_aers_in)) then
         call mpas_log_write('********************************************************************************', messageType=MPAS_LOG_ERR)
         call mpas_log_write('* Error while setting up packages for init_atmosphere core.',                      messageType=MPAS_LOG_ERR)
         call mpas_log_write('********************************************************************************', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (config_init_case == 8) then
         initial_conds = .false.
         sfc_update = .true.
      else
         initial_conds = .true.
         sfc_update = .false.
      end if

      if (config_init_case == 9) then
         lbcs = .true.
         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if(lexist) mp_thompson_aers_in = .true.
      else
         lbcs = .false.
         mp_thompson_aers_in = .false.
      end if

      if (config_init_case == 7) then

         !
         ! The logic here is a little convoluted
         ! For input, we want to read in fields from all earlier stages, except if those earlier stages are being run now
         ! For output, we want to output the fields that were computed in a stage and all of those from earlier stages
         !
         gwd_stage_in = config_native_gwd_static .and. &
                        (.not. config_static_interp)
         gwd_stage_out = config_native_gwd_static
         gwd_gsl_stage_out = config_native_gwd_gsl_static
         vertical_stage_in = config_vertical_grid .and. &
                             (.not. config_native_gwd_static) .and. &
                             (.not. config_static_interp)
         vertical_stage_out = config_vertical_grid
         met_stage_in = config_met_interp .and. &
                        (.not. config_native_gwd_static) .and. &
                        (.not. config_static_interp) .and. &
                        (.not. config_vertical_grid)
         met_stage_out = config_met_interp

         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if((lexist .and. met_stage_out) .or. (lexist .and. met_stage_in)) mp_thompson_aers_in = .true.

      else if (config_init_case == 8) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .true.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .false.

      !
      ! When interpolating LBC fields, we need all inputs that would be needed for the interpolation
      ! of ICs, so met_stage_in = .true.
      !
      else if (config_init_case == 9) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .true.
         met_stage_out = .true.

         mp_thompson_aers_in = .false.
         inquire(file="QNWFA_QNIFA_SIGMA_MONTHLY.dat",exist=lexist)
         if((lexist .and. met_stage_out) .or. (lexist .and. met_stage_in)) mp_thompson_aers_in = .true.

         initial_conds = .false.   ! Also, turn off the initial_conds package to avoid writing the IC "output" stream

      else if (config_init_case == 13) then
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .true.
         met_stage_in = .false.
         met_stage_out = .false.

      else
         gwd_stage_in = .false.
         gwd_stage_out = .false.
         gwd_gsl_stage_out = .false.
         vertical_stage_in = .false.
         vertical_stage_out = .false.
         met_stage_in = .false.
         met_stage_out = .true.
      end if

      !
      ! Package for 3-d first-guess atmospheric and land-surface fields is only active if
      ! we are interpolating real-data ICs or LBCs
      !
      nullify(first_guess_field)
      call mpas_pool_get_package(packages, 'first_guess_fieldActive', first_guess_field)

      first_guess_field = .false.
      if ((config_init_case == 7 .and. config_met_interp) .or. config_init_case == 9) then
         first_guess_field = .true.
      end if

      !
      ! Noah-MP
      !
      nullify(config_noahmp_static)
      call mpas_pool_get_config(configs, 'config_noahmp_static', config_noahmp_static)

      nullify(noahmp)
      call mpas_pool_get_package(packages, 'noahmpActive', noahmp)

      if (associated(config_noahmp_static) .and. associated(noahmp)) then
         noahmp = config_noahmp_static
      else
         call mpas_log_write('********************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('* Error while setting up packages for init_atmosphere core:',                      &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*    Either the ''noahmp'' package or the ''config_noahmp_static'' namelist',      &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('*    option is not defined.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write('********************************************************************************', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

   end function init_atm_setup_packages


   !***********************************************************************
   !
   !  function init_atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_timekeeping, only : mpas_set_time, mpas_set_timeInterval, mpas_create_clock
      use mpas_derived_types, only : MPAS_Clock_type, MPAS_Time_type, MPAS_Timeinterval_type, mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config
      use mpas_kind_types, only : StrKIND

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_start_time, config_stop_time
      integer, pointer :: config_fg_interval

      type (MPAS_Time_type)         :: start_time, stop_time
      type (MPAS_TimeInterval_type) :: dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_fg_interval', config_fg_interval)

      call mpas_set_time(start_time, dateTimeString=trim(config_start_time))
      call mpas_set_time(stop_time, dateTimeString=trim(config_stop_time))
      call mpas_set_timeInterval(dt, S=config_fg_interval)
      call mpas_create_clock(core_clock, start_time, dt, stopTime=stop_time)

   end function init_atm_setup_clock


   !***********************************************************************
   !
   !  function init_atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types, only : mpas_log_type, domain_type
      use mpas_log, only : mpas_log_init, mpas_log_open
      use mpas_framework, only : mpas_framework_report_settings

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo  !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain  !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call mpas_log_write('')
      call mpas_log_write('MPAS Init-Atmosphere Version '//trim(domain % core % modelVersion))
      call mpas_log_write('')

      call mpas_framework_report_settings(domain)

   end function init_atm_setup_log!}}}


   !***********************************************************************
   !
   !  function init_atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function init_atm_get_mesh_stream(configs, streamInfo, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      ierr = 0

      write(stream,'(a)') 'input'

   end function init_atm_get_mesh_stream


   !***********************************************************************
   !
   !  function init_atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_decompositions(decompList, methodList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR, &
                                     mpas_dim_size_fn_list, mpas_computed_dim_size_fn, MPAS_COMPUTED_DIM_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp 
      use mpas_computed_dims, only : mpas_computed_dims_create_list, mpas_computed_dims_register_method

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      type (mpas_dim_size_fn_list), pointer :: methodList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc
      procedure (mpas_computed_dim_size_fn), pointer :: dimSizeFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      else
         ierr = 1
         return
      end if

      call mpas_computed_dims_create_list(methodList)

      dimSizeFunc => cien

      call mpas_computed_dims_register_method(methodList, 'cien', dimSizeFunc, ierr)

      if ( ierr == MPAS_COMPUTED_DIM_NOERR ) then
         ierr = 0
      else
         ierr = 1
         return
      end if

   end function init_atm_setup_decompositions


   !***********************************************************************
   !
   !  function init_atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be a
   !>  init atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function init_atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call init_atm_generate_structs(block, block % structs, block % dimensions, block % packages)

   end function init_atm_setup_block


   function cien(block, manager) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager

      integer :: dim_size

      integer :: ierr
      character(len=StrKIND) :: mesh_stream, mesh_filename
      type (MPAS_streamInfo_type) :: streamInfo

      integer, pointer :: nCellsSolve

      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:), pointer :: bdyMaskCell
      integer :: fieldtype

      type (MPAS_IO_Handle_type) :: handle

      integer :: iCell
      integer :: nBdyCells_local, nBdyCells


      ierr = init_atm_get_mesh_stream(block % configs, streamInfo, mesh_stream)
      call mpas_log_write('Mesh stream is '//trim(mesh_stream))

      call MPAS_get_stream_filename(manager, trim(mesh_stream), filename=mesh_filename)
      call mpas_log_write('Mesh filename is '//trim(mesh_filename))


      nullify(indexToCellID)
      call mpas_pool_get_array(block % allFields, 'indexToCellID_blk', indexToCellID)

      nullify(nCellsSolve)
      call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

      if (associated(nCellsSolve)) then
          call mpas_log_write('Found nCellsSolve $i', intArgs=[nCellsSolve])
      else
          call mpas_log_write('Could not find nCellsSolve')
      end if

      allocate(bdyMaskCell(nCellsSolve))

      if (associated(indexToCellID)) then
          call mpas_log_write('Found indexToCellID')
      else
          call mpas_log_write('Could not find indexToCellID')
      end if

      handle = MPAS_io_open(trim(mesh_filename), MPAS_IO_READ, MPAS_IO_NETCDF, manager % ioContext, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error opening '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK opening '//trim(mesh_filename))
      end if

      call MPAS_io_inq_var(handle, 'bdyMaskCell', fieldtype=fieldtype, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error inquiring about bdyMaskCell', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK inquiring about bdyMaskCell')
      end if

      call MPAS_io_set_var_indices(handle, 'bdyMaskCell', indexToCellID(1:nCellsSolve), ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error setting indices', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK setting indices')
      end if

      call MPAS_io_get_var(handle, 'bdyMaskCell', bdyMaskCell, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error reading bdyMaskCell', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK reading bdyMaskCell')
      end if

      call MPAS_io_close(handle)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error closing '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK closing '//trim(mesh_filename))
      end if

      call mpas_log_write('Min/max bdyMaskCell=$i $i', intArgs=[minval(bdyMaskCell), maxval(bdyMaskCell)])

      nBdyCells_local = 0
      do iCell = 1, nCellsSolve
          if (bdyMaskCell(iCell) > 0) nBdyCells_local = nBdyCells_local + 1
      end do

      call mpas_dmpar_sum_int(block % domain % dminfo, nBdyCells_local, nBdyCells)

      call mpas_log_write('nBdyCells = $i', intArgs=[nBdyCells])

      deallocate(bdyMaskCell)

      dim_size = nBdyCells

   end function cien


   function uniform_decomp(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr
      
      integer :: blockDimSize, blockStart, i


      ! -- for testing stream reads
      character(len=StrKIND) :: mesh_stream, mesh_filename
      type (MPAS_streamInfo_type) :: streamInfo

      type (field1DInteger), pointer :: bdyMaskCellField
      integer, dimension(:), pointer :: bdyMaskCell

      integer, pointer :: nCellsSolve

      type (MPAS_Stream_type) :: stream

      
      iErr = 0
      
      blockDimSize = globalDimSize / numBlocks
      if ( block % blockID < mod(globalDimSize, numBlocks) ) then
         blockDimSize = blockDimSize + 1
      end if

      blockStart = block % blockID * blockDimSize
      if ( block % blockID >= mod(globalDimSize, numBlocks) ) then
         blockStart = blockStart + mod(globalDimSize, numBlocks)
      end if

      allocate(ownedIndices(blockDimSize))
      do i = 1, blockDimSize
         ownedIndices(i) = blockStart + i
      end do

      ierr = init_atm_get_mesh_stream(block % configs, streamInfo, mesh_stream)
      call mpas_log_write('Mesh stream is '//trim(mesh_stream))

      call MPAS_get_stream_filename(manager, trim(mesh_stream), filename=mesh_filename)
      call mpas_log_write('Mesh filename is '//trim(mesh_filename))

      nullify(bdyMaskCellField)
      call mpas_pool_get_field(block % allFields, 'bdyMaskCell', bdyMaskCellField)

      nullify(bdyMaskCell)
      call mpas_pool_get_array(block % allFields, 'indexToCellID_blk', bdyMaskCell)

      nullify(nCellsSolve)
      call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

      if (associated(nCellsSolve)) then
          call mpas_log_write('Found nCellsSolve $i', intArgs=[nCellsSolve])
      else
          call mpas_log_write('Could not find nCellsSolve')
      end if

      if (associated(bdyMaskCellField)) then
          call mpas_log_write('Found bdyMaskCellField')
      else
          call mpas_log_write('Could not find bdyMaskCellField')
      end if

      if (associated(bdyMaskCell)) then
          call mpas_log_write('Found bdyMaskCell')
      else
          call mpas_log_write('Could not find bdyMaskCell')
      end if

!      call MPAS_createStream(stream, block % domain % ioContext, trim(mesh_filename), MPAS_IO_NETCDF, MPAS_IO_READ)
!      call MPAS_streamAddField(stream, bdyMaskCellField)
!      call MPAS_readStream(stream, 1)
!      call MPAS_closeStream(stream)

      call mpas_log_write('nCellsSolve = $i', intArgs=[nCellsSolve])
      call mpas_log_write('Min/max bdyMaskCell=$i $i', &
                          intArgs=[minval(bdyMaskCell(1:nCellsSolve)), maxval(bdyMaskCell(1:nCellsSolve))])
      call mpas_log_write('size(bdyMaskCell)=$i', intArgs=[size(bdyMaskCell)])
      
   end function uniform_decomp


#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "setup_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module init_atm_core_interface
