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
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method
      use mpas_computed_dims, only : mpas_computed_dims_create_list, mpas_computed_dims_register_method

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      type (mpas_dim_size_fn_list), pointer :: methodList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc
      procedure (mpas_computed_dim_size_fn), pointer :: dimSizeFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => boundary_cells
      call mpas_decomp_register_method(decompList, 'boundary_cells', decompFunc, ierr)

      decompFunc => boundary_edges
      call mpas_decomp_register_method(decompList, 'boundary_edges', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      else
         ierr = 1
         return
      end if

      call mpas_computed_dims_create_list(methodList)

      dimSizeFunc => num_boundary_cells
      call mpas_computed_dims_register_method(methodList, 'num_boundary_cells', dimSizeFunc, ierr)

      dimSizeFunc => num_boundary_edges
      call mpas_computed_dims_register_method(methodList, 'num_boundary_edges', dimSizeFunc, ierr)

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


   function num_boundary_cells(block, manager) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager

      integer :: dim_size

      dim_size = num_boundary_elems(block, manager, 'Cell', 'Cells')

   end function num_boundary_cells


   function num_boundary_edges(block, manager) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager

      integer :: dim_size

      dim_size = num_boundary_elems(block, manager, 'Edge', 'Edges')

   end function num_boundary_edges


   function num_boundary_elems(block, manager, elemType, elemTypePlural) result(dim_size)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural

      integer :: dim_size

      integer :: ierr
      character(len=StrKIND) :: mesh_stream, mesh_filename
      type (MPAS_streamInfo_type) :: streamInfo

      integer, pointer :: nElemsSolve

      integer, dimension(:), pointer :: indexToElemID
      integer, dimension(:), pointer :: bdyMaskElem
      integer :: fieldtype

      type (MPAS_IO_Handle_type) :: handle

      integer :: iElem
      integer :: nBdyElems_local, nBdyElems


      call read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

!      call mpas_log_write('Min/max bdyMaskElem=$i $i', intArgs=[minval(bdyMaskElem), maxval(bdyMaskElem)])

      nBdyElems_local = 0
      do iElem = 1, nElemsSolve
          if (bdyMaskElem(iElem) > 0) nBdyElems_local = nBdyElems_local + 1
      end do

      call mpas_dmpar_sum_int(block % domain % dminfo, nBdyElems_local, nBdyElems)

!      call mpas_log_write('nBdyElems = $i', intArgs=[nBdyElems])

      deallocate(bdyMaskElem)

      dim_size = nBdyElems

   end function num_boundary_elems


   function boundary_cells(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_timer

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr
      
      call mpas_timer_start('boundary_cells')
      iErr = boundary_elems(block, manager, globalDimSize, numBlocks, ownedIndices, 'Cell', 'Cells')
      call mpas_timer_stop('boundary_cells')

   end function boundary_cells


   function boundary_edges(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_timer

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks

      integer, dimension(:), pointer :: ownedIndices

      integer :: iErr
      
      call mpas_timer_start('boundary_edges')
      iErr = boundary_elems(block, manager, globalDimSize, numBlocks, ownedIndices, 'Edge', 'Edges')
      call mpas_timer_stop('boundary_edges')

   end function boundary_edges


   function boundary_elems(block, manager, globalDimSize, numBlocks, ownedIndices, elemType, elemTypePlural) result(iErr)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams
      use mpi_f08
      use mpas_sort

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      integer, intent(in) :: globalDimSize
      integer, intent(in) :: numBlocks
      integer, dimension(:), pointer :: ownedIndices
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural

      integer :: iErr
      
      character(len=StrKIND) :: mesh_stream, mesh_filename
      type (MPAS_streamInfo_type) :: streamInfo

      integer, pointer :: nElemsSolve

      integer, dimension(:), pointer :: indexToElemID
      integer, dimension(:), pointer :: bdyMaskElem
      integer :: fieldtype

      type (MPAS_IO_Handle_type) :: handle

      integer :: iElem
      integer :: nBdyElems_local, nBdyElems

      type (MPI_Comm) :: comm
      type (MPI_Request) :: req
      integer :: nRanks, myRank

      integer, dimension(:,:), allocatable :: temp

      integer :: i, j, k
      integer :: n, n_next
      integer, dimension(:), pointer :: bdyIDs, bdyIDs_next
      integer :: srcRank, dstRank


      iErr = 0
      
      call read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

      nBdyElems = 0
      do iElem = 1, nElemsSolve
          if (bdyMaskElem(iElem) > 0) nBdyElems = nBdyElems + 1
      end do

!      call mpas_dmpar_sum_int(block % domain % dminfo, nBdyElems_local, nBdyElems)

!      if (nBdyElems > 0) then
!          call mpas_log_write('Again, nBdyElems = $i', intArgs=[nBdyElems])
!      else
!          call mpas_log_write('No boundary Elems')
!      end if

      ! comm = block % domain % dminfo % comm
      if (nBdyElems > 0) then
          call MPI_Comm_split(block % domain % dminfo % comm, 1, block % domain % dminfo % my_proc_id, comm, ierr)
      else
          call MPI_Comm_split(block % domain % dminfo % comm, 0, block % domain % dminfo % my_proc_id, comm, ierr)
      end if
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to split comm', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('Successfully split comm', intArgs=[nRanks])
      end if

if (nBdyElems > 0) then
      call MPI_Comm_size(comm, nRanks, ierr)
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to get comm size', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('nRanks = $i', intArgs=[nRanks])
      end if

      call MPI_Comm_rank(comm, myRank, ierr)
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to get comm rank', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('myRank = $i', intArgs=[myRank])
      end if

      allocate(temp(2,nBdyElems))
      i = 1
      do iElem = 1, nElemsSolve
          if (bdyMaskElem(iElem) > 0) then
              temp(1,i) = indexToElemID(iElem)   ! Parent element
              temp(2,i) = i                      ! Compute order on task
              i = i + 1
          end if
      end do

      call mpas_quicksort(nBdyElems, temp)

      allocate(ownedIndices(nBdyElems))
      do iElem = 1, nBdyElems
          ownedIndices(iElem) = iElem            ! Tentative boundary-only offset
      end do

      srcRank = mod(myRank + nRanks - 1, nRanks)
      dstRank = mod(myRank + 1, nRanks)
      n = nBdyElems
      allocate(bdyIDs(n))

      bdyIDs(:) = temp(1,:)     ! Full domain boundary offsets

      do i = 1, nRanks - 1
          call MPI_Irecv(n_next, 1, MPI_INTEGER, srcRank, srcRank, comm, req, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Recv from $i failed', intArgs=[srcRank], messageType=MPAS_LOG_ERR)

          call MPI_Send(n, 1, MPI_INTEGER, dstRank, myRank, comm, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Send to $i failed', intArgs=[dstRank], messageType=MPAS_LOG_ERR)

          call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Request failed', messageType=MPAS_LOG_ERR)

!          call mpas_log_write('Now handling $i boundary elements...', intArgs=[n_next])

          allocate(bdyIDs_next(n_next))

! Optimization: if n/n_next is 0, no need to send/recv
! Once we have a separate communicator for tasks that have bdy elements, this will not apply

          call MPI_Irecv(bdyIDs_next, n_next, MPI_INTEGER, srcRank, srcRank, comm, req, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Recv from $i failed', intArgs=[srcRank], messageType=MPAS_LOG_ERR)

          call MPI_Send(bdyIDs, n, MPI_INTEGER, dstRank, myRank, comm, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Send to $i failed', intArgs=[dstRank], messageType=MPAS_LOG_ERR)

          call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
          if (ierr /= MPI_SUCCESS) call mpas_log_write('Request failed', messageType=MPAS_LOG_ERR)

          deallocate(bdyIDs)
          bdyIDs => bdyIDs_next
          n = n_next

#if 1
          j = 0
          do k = 1, nBdyElems
              do while (j < n)
                  if (bdyIDs(j+1) > temp(1,k)) exit
                  j = j + 1
              end do
              ownedIndices(k) = ownedIndices(k) + j
          end do
#else
          do j = 1, n
              do k = nBdyElems, 1, -1
                  if (bdyIDs(j) < temp(1,k)) then
                      ownedIndices(k) = ownedIndices(k) + 1
                  else
                      exit
                  end if
              end do
          end do
#endif

      end do
      deallocate(bdyIDs)

      ! Place elements back in original order
      temp(1,:) = ownedIndices(:)                ! Tentative boundary-only offset
      do iElem = 1, nBdyElems
          ownedIndices(temp(2,iElem)) = temp(1,iElem)
      end do

      deallocate(temp)

else

      allocate(ownedIndices(nBdyElems))

end if
      deallocate(bdyMaskElem)

      call MPI_Comm_free(comm, ierr)
      if (ierr /= MPI_SUCCESS) then
          call mpas_log_write('Failed to free comm', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('Successfully freed comm', intArgs=[nRanks])
      end if

   end function boundary_elems


   subroutine read_bdy_mask(block, manager, elemType, elemTypePlural, indexToElemID, bdyMaskElem, nElemsSolve)

      use mpas_stream_manager, only : MPAS_get_stream_filename
      use mpas_io_streams

      type (block_type), intent(in) :: block
      type (mpas_streamManager_type), intent(inout) :: manager
      character(len=*), intent(in) :: elemType
      character(len=*), intent(in) :: elemTypePlural
      integer, dimension(:), pointer :: indexToElemID  ! Output
      integer, dimension(:), pointer :: bdyMaskElem  ! Output
      integer, pointer :: nElemsSolve  ! Output

      integer :: ierr
      character(len=StrKIND) :: mesh_stream, mesh_filename
      type (MPAS_streamInfo_type) :: streamInfo

      integer :: fieldtype

      type (MPAS_IO_Handle_type) :: handle


      ierr = init_atm_get_mesh_stream(block % configs, streamInfo, mesh_stream)
      call mpas_log_write('Mesh stream is '//trim(mesh_stream))

      call MPAS_get_stream_filename(manager, trim(mesh_stream), filename=mesh_filename)
      call mpas_log_write('Mesh filename is '//trim(mesh_filename))

      nullify(indexToElemID)
      call mpas_pool_get_array(block % allFields, 'indexTo'//elemType//'ID_blk', indexToElemID)

      nullify(nElemsSolve)
      call mpas_pool_get_dimension(block % dimensions, 'n'//elemTypePlural//'Solve', nElemsSolve)

      if (associated(nElemsSolve)) then
          call mpas_log_write('Found nElemsSolve $i', intArgs=[nElemsSolve])
      else
          call mpas_log_write('Could not find nElemsSolve')
      end if

      allocate(bdyMaskElem(nElemsSolve))

      if (associated(indexToElemID)) then
          call mpas_log_write('Found indexToElemID')
      else
          call mpas_log_write('Could not find indexToElemID')
      end if

      handle = MPAS_io_open(trim(mesh_filename), MPAS_IO_READ, MPAS_IO_NETCDF, manager % ioContext, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error opening '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK opening '//trim(mesh_filename))
      end if

      call MPAS_io_inq_var(handle, 'bdyMask'//elemType, fieldtype=fieldtype, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error inquiring about bdyMaskElem', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK inquiring about bdyMaskElem')
      end if

      call MPAS_io_set_var_indices(handle, 'bdyMask'//elemType, indexToElemID(1:nElemsSolve), ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error setting indices', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK setting indices')
      end if

      call MPAS_io_get_var(handle, 'bdyMask'//elemType, bdyMaskElem, ierr=ierr)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error reading bdyMaskElem', messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK reading bdyMaskElem')
      end if

      call MPAS_io_close(handle)
      if (ierr /= MPAS_IO_NOERR) then
          call mpas_log_write('MGD  - Error closing '//trim(mesh_filename), messageType=MPAS_LOG_ERR)
      else
          call mpas_log_write('MGD - OK closing '//trim(mesh_filename))
      end if

   end subroutine read_bdy_mask

#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "setup_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module init_atm_core_interface
