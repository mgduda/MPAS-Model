! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_core_interface

   use mpas_attlist
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units


   contains


   !***********************************************************************
   !
   !  routine atm_setup_core
   !
   !> \brief   Atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use atm_core, only : atm_core_init, atm_core_run, atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => atm_core_init
      core % core_run => atm_core_run
      core % core_finalize => atm_core_finalize
      core % define_packages => atm_define_packages
      core % setup_packages => atm_setup_packages
      core % setup_decompositions => atm_setup_decompositions
      core % setup_clock => atm_setup_clock
      core % setup_log => atm_setup_log
      core % get_mesh_stream => atm_get_mesh_stream
      core % setup_immutable_streams => atm_setup_immutable_streams
      core % setup_derived_dimensions => atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => atm_setup_decomposed_dimensions
      core % setup_block => atm_setup_block
      core % setup_namelist => atm_setup_namelists
      core % component_role => atm_component_role
      core % component_role_leader => atm_component_role_leader
      core % component_role_is => atm_component_role_is
      core % component_role_includes => atm_component_role_includes

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "core_variables.inc"

   end subroutine atm_setup_core


   !***********************************************************************
   !
   !  routine atm_component_role
   !
   !> \brief Define roles for each MPI task and OpenMP thread
   !> \author Michael Duda, Supreeth Madapur Suresh
   !> \date 1 July 2020
   !> \details 
   !>  This routine returns an integer value representing the "role" for
   !>  any given MPI task and OpenMP thread.
   !>
   !>  Upon its first invocation, this routine will attempt to read two
   !>  environemnt variables: MPAS_RADIATION_RANKS_PER_NODE and
   !>  MPAS_DYNAMICS_RANKS_PER_NODE, which will be used to determine how
   !>  roles are assigned to MPI ranks.
   !>
   !>  If MPAS-A was compiled with OpenACC support, this routine will also
   !>  set OpenACC devices for each MPI rank on its first invocation.
   !
   !-----------------------------------------------------------------------
   integer function atm_component_role(self)

      use mpas_derived_types, only : mpas_coupler_type
      use mpas_kind_types, only : StrKIND
      use mpas_abort, only : mpas_dmpar_global_abort
      use mpas_intracoupler, only : ROLE_INTEGRATE, ROLE_RADIATION

#ifdef _MPI
#ifndef NOMPIMOD
      use mpi, only : MPI_COMM_TYPE_SHARED, MPI_INFO_NULL
#endif
#endif

#ifdef MPAS_OPENACC
      use openacc, only : acc_get_num_devices, acc_set_device_num, acc_set_device_type, &
                          acc_device_nvidia, acc_device_host
#endif

      implicit none

      class (mpas_coupler_type), intent(in) :: self

#ifdef _MPI
#ifdef NOMPIMOD
      include 'mpif.h'
#endif
#endif

      integer, save :: local_rank = -1
      integer, save :: role

      integer :: radiation_ranks_per_node, dynamics_ranks_per_node

      character(len=StrKIND) :: radiation_ranks_per_node_str
      character(len=StrKIND) :: dynamics_ranks_per_node_str

      integer :: ierr
      integer :: istat_rad, istat_dyn
      integer :: local_size
      integer :: on_node_comm
      integer :: socket_rank

#ifdef MPAS_OPENACC
      integer :: num_devices, my_device
      integer :: role_comm, role_rank, role_size
#endif


      !
      ! If this routine has not yet been called, attempt to read environment variables,
      ! and set local variable 'role' to the role of the calling task
      !
      if (local_rank == -1) then
          call get_environment_variable('MPAS_RADIATION_RANKS_PER_NODE', &
                                        value=radiation_ranks_per_node_str, status=istat_rad)
          if (istat_rad == 0) then
              read(radiation_ranks_per_node_str, *) radiation_ranks_per_node
          end if

          call get_environment_variable('MPAS_DYNAMICS_RANKS_PER_NODE', &
                                        value=dynamics_ranks_per_node_str, status=istat_dyn)
          if (istat_dyn == 0) then
              read(dynamics_ranks_per_node_str, *) dynamics_ranks_per_node
          end if

          if (istat_rad == 0 .and. istat_dyn == 0) then  ! Both env vars were successfully read

              !
              ! Determine rank within a shared-memory node
              !
#ifdef _MPI
              call MPI_Comm_split_type(self % intercomm, MPI_COMM_TYPE_SHARED, self % intercomm_rank, &
                                       MPI_INFO_NULL, on_node_comm, ierr)
              call MPI_Comm_rank(on_node_comm, local_rank, ierr)
              call MPI_Comm_size(on_node_comm, local_size, ierr)
#else
              local_rank = 0
              local_size = 1
#endif
              !
              ! Evenly distribute role assignments over two sockets
              !
              if (local_rank < local_size / 2) then
                  socket_rank = local_rank
              else
                  socket_rank = local_rank - (local_size / 2)
              end if

              if (mod(socket_rank, 2) == 0) then
                  role = ROLE_RADIATION
              else
                  role = ROLE_INTEGRATE
              end if

              !
              ! Roles are assigned within a socket alternating between ROLE_RADIATION and
              ! ROLE_INTEGRATE. If the number of dynamics and radiation ranks per node are not
              ! equal, all "excess" ranks on a socket beyond min(dynamics_ranks_per_node, radiation_ranks_per_node)
              ! will be assigned either ROLE_INTEGRATE or ROLE_RADIATION as needed.
              !
              if (socket_rank >= min(dynamics_ranks_per_node, radiation_ranks_per_node)) then
                  if (dynamics_ranks_per_node > radiation_ranks_per_node) then
                      role = ROLE_INTEGRATE
                  else
                      role = ROLE_RADIATION
                  end if
              end if


              !
              ! Assign devices based on rank in the role within a shared-memory node
              !
#ifdef MPAS_OPENACC
#ifdef _MPI
              call MPI_Comm_split(on_node_comm, role, local_rank, role_comm, ierr)
              call MPI_Comm_rank(role_comm, role_rank, ierr)
              call MPI_Comm_size(role_comm, role_size, ierr)
#else
              role_rank = 0
              role_size = 1
#endif

              if (role == ROLE_INTEGRATE) then
                  num_devices = acc_get_num_devices(acc_device_nvidia)
                  my_device = (role_rank * num_devices) / role_size
                  call acc_set_device_num(my_device, acc_device_nvidia)
              else
                  call acc_set_device_type(acc_device_host)
              end if
#endif

          else ! At least one environment variable could not be successfully read

              local_rank = self % intercomm_rank
              role = ior(ROLE_INTEGRATE, ROLE_RADIATION)

          end if

      end if

      atm_component_role = role

   end function atm_component_role


   !***********************************************************************
   !
   !  routine atm_component_role_leader
   !
   !> \brief Identifies the "leader" task for a given role
   !> \author Michael Duda
   !> \date 3 January 2018
   !> \details 
   !>  Given a role, this function returns the MPI task that acts as the leader
   !>  for that role. This is needed to establish MPI intercommunicators
   !>  between two roles.
   !
   !-----------------------------------------------------------------------
   integer function atm_component_role_leader(self, role)

      use mpas_derived_types, only : mpas_coupler_type
      use mpas_intracoupler, only : ROLE_RADIATION

      implicit none

      class (mpas_coupler_type), intent(in) :: self
      integer, intent(in) :: role


      !
      ! Assumption: the first MPI rank always has the role of ROLE_RADIATION
      ! (possibly among other roles), and the second MPI rank always has the role
      ! of ROLE_INTEGRATE (possibly among other roles)
      !

      if (iand(role, ROLE_RADIATION) == ROLE_RADIATION) then
          atm_component_role_leader = 0
      else
          atm_component_role_leader = 1
      end if

   end function atm_component_role_leader


   !***********************************************************************
   !
   !  routine atm_component_role_is
   !
   !> \brief Indicates whether a task has a given role exclusively
   !> \author Michael Duda
   !> \date 11 June 2020
   !> \details 
   !>  Given a role, this function returns .true. if the task is assigned
   !>  that role exclusively, and .false. otherwise. To determine if a task
   !>  is assigned a role, possibly among others, use the role_includes
   !>  function instead.
   !
   !-----------------------------------------------------------------------
   logical function atm_component_role_is(self, role)

      use mpas_derived_types, only : mpas_coupler_type
      use mpas_intracoupler, only : ROLE_INTEGRATE

      implicit none

      class (mpas_coupler_type), intent(in) :: self
      integer, intent(in) :: role

      atm_component_role_is = (self % role() == role)

   end function atm_component_role_is


   !***********************************************************************
   !
   !  routine atm_component_role_includes
   !
   !> \brief Indicates whether a task has a given role
   !> \author Michael Duda
   !> \date 11 June 2020
   !> \details 
   !>  Given a role, this function returns .true. if the task is assigned
   !>  that role, possibly among others, and .false. otherwise. To determine
   !>  if a task is assigned a role exclusively, use the role_is function
   !>  instead.
   !
   !-----------------------------------------------------------------------
   logical function atm_component_role_includes(self, role)

      use mpas_derived_types, only : mpas_coupler_type
      use mpas_intracoupler, only : ROLE_INTEGRATE

      implicit none

      class (mpas_coupler_type), intent(in) :: self
      integer, intent(in) :: role

      atm_component_role_includes = (iand(self % role(), role) == role)

   end function atm_component_role_includes


   !***********************************************************************
   !
   !  routine atm_setup_domain
   !
   !> \brief   Atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "domain_variables.inc"

   end subroutine atm_setup_domain


   !***********************************************************************
   !
   !  function atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function atm_setup_packages(configs, packages, iocontext) result(ierr)

      use mpas_dmpar
      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

#ifdef DO_PHYSICS
      use mpas_atmphys_control
      use mpas_atmphys_packages
#endif

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      logical, pointer :: iauActive
      character(len=StrKIND), pointer :: config_iau_option
      logical, pointer :: limited_areaActive
      logical, pointer :: config_apply_lbcs
      logical, pointer :: config_jedi_da, jedi_daActive
      integer :: local_ierr

      integer, pointer :: config_num_passive_tracers

      logical, pointer :: tracer1Active
      logical, pointer :: tracer2Active
      logical, pointer :: tracer3Active
      logical, pointer :: tracer4Active
      logical, pointer :: tracer5Active
      logical, pointer :: tracer6Active
      logical, pointer :: tracer7Active
      logical, pointer :: tracer8Active
      logical, pointer :: tracer9Active
      logical, pointer :: tracer10Active
      logical, pointer :: tracer11Active
      logical, pointer :: tracer12Active
      logical, pointer :: tracer13Active
      logical, pointer :: tracer14Active
      logical, pointer :: tracer15Active
      logical, pointer :: tracer16Active
      logical, pointer :: tracer17Active
      logical, pointer :: tracer18Active
      logical, pointer :: tracer19Active
      logical, pointer :: tracer20Active
      logical, pointer :: tracer21Active
      logical, pointer :: tracer22Active
      logical, pointer :: tracer23Active
      logical, pointer :: tracer24Active
      logical, pointer :: tracer25Active
      logical, pointer :: tracer26Active
      logical, pointer :: tracer27Active
      logical, pointer :: tracer28Active
      logical, pointer :: tracer29Active
      logical, pointer :: tracer30Active
      logical, pointer :: tracer31Active
      logical, pointer :: tracer32Active
      logical, pointer :: tracer33Active      

      ierr = 0

      !
      ! Incremental analysis update
      !
      nullify(config_iau_option)
      call mpas_pool_get_config(configs, 'config_IAU_option', config_iau_option)

      nullify(iauActive)
      call mpas_pool_get_package(packages, 'iauActive', iauActive)

      if (trim(config_iau_option) /= 'off') then
         iauActive = .true.
      else
         iauActive = .false.
      end if

            !
      ! Set up packages for passive tracers
      !
      nullify(config_num_passive_tracers)
      call mpas_pool_get_config(configs, 'config_num_passive_tracers', config_num_passive_tracers)

      nullify(tracer1Active)
      nullify(tracer2Active)
      nullify(tracer3Active)
      nullify(tracer4Active)
      nullify(tracer5Active)
      nullify(tracer6Active)
      nullify(tracer7Active)
      nullify(tracer8Active)
      nullify(tracer9Active)
      nullify(tracer10Active)
      nullify(tracer11Active)
      nullify(tracer12Active)
      nullify(tracer13Active)
      nullify(tracer14Active)
      nullify(tracer15Active)
      nullify(tracer16Active)
      nullify(tracer17Active)
      nullify(tracer18Active)
      nullify(tracer19Active)
      nullify(tracer20Active)
      nullify(tracer21Active)
      nullify(tracer22Active)
      nullify(tracer23Active)
      nullify(tracer24Active)
      nullify(tracer25Active)
      nullify(tracer26Active)
      nullify(tracer27Active)
      nullify(tracer28Active)
      nullify(tracer29Active)
      nullify(tracer30Active)
      nullify(tracer31Active)
      nullify(tracer32Active)
      nullify(tracer33Active)          

      call mpas_pool_get_package(packages, 'tracer1Active',  tracer1Active)
      call mpas_pool_get_package(packages, 'tracer2Active',  tracer2Active)
      call mpas_pool_get_package(packages, 'tracer3Active',  tracer3Active)
      call mpas_pool_get_package(packages, 'tracer4Active',  tracer4Active)
      call mpas_pool_get_package(packages, 'tracer5Active',  tracer5Active)
      call mpas_pool_get_package(packages, 'tracer6Active',  tracer6Active)
      call mpas_pool_get_package(packages, 'tracer7Active',  tracer7Active)
      call mpas_pool_get_package(packages, 'tracer8Active',  tracer8Active)
      call mpas_pool_get_package(packages, 'tracer9Active',  tracer9Active)
      call mpas_pool_get_package(packages, 'tracer10Active', tracer10Active)
      call mpas_pool_get_package(packages, 'tracer11Active',  tracer11Active)
      call mpas_pool_get_package(packages, 'tracer12Active',  tracer12Active)
      call mpas_pool_get_package(packages, 'tracer13Active',  tracer13Active)
      call mpas_pool_get_package(packages, 'tracer14Active',  tracer14Active)
      call mpas_pool_get_package(packages, 'tracer15Active',  tracer15Active)
      call mpas_pool_get_package(packages, 'tracer16Active',  tracer16Active)
      call mpas_pool_get_package(packages, 'tracer17Active',  tracer17Active)
      call mpas_pool_get_package(packages, 'tracer18Active',  tracer18Active)
      call mpas_pool_get_package(packages, 'tracer19Active',  tracer19Active)
      call mpas_pool_get_package(packages, 'tracer20Active', tracer20Active)
      call mpas_pool_get_package(packages, 'tracer21Active',  tracer21Active)
      call mpas_pool_get_package(packages, 'tracer22Active',  tracer22Active)
      call mpas_pool_get_package(packages, 'tracer23Active',  tracer23Active)
      call mpas_pool_get_package(packages, 'tracer24Active',  tracer24Active)
      call mpas_pool_get_package(packages, 'tracer25Active',  tracer25Active)
      call mpas_pool_get_package(packages, 'tracer26Active',  tracer26Active)
      call mpas_pool_get_package(packages, 'tracer27Active',  tracer27Active)
      call mpas_pool_get_package(packages, 'tracer28Active',  tracer28Active)
      call mpas_pool_get_package(packages, 'tracer29Active',  tracer29Active)
      call mpas_pool_get_package(packages, 'tracer30Active', tracer30Active)
      call mpas_pool_get_package(packages, 'tracer31Active',  tracer31Active)
      call mpas_pool_get_package(packages, 'tracer32Active',  tracer32Active)
      call mpas_pool_get_package(packages, 'tracer33Active',  tracer33Active)

      tracer1Active  = .false.
      tracer2Active  = .false.
      tracer3Active  = .false.
      tracer4Active  = .false.
      tracer5Active  = .false.
      tracer6Active  = .false.
      tracer7Active  = .false.
      tracer8Active  = .false.
      tracer9Active  = .false.
      tracer10Active = .false.
      tracer11Active  = .false.
      tracer12Active  = .false.
      tracer13Active  = .false.
      tracer14Active  = .false.
      tracer15Active  = .false.
      tracer16Active  = .false.
      tracer17Active  = .false.
      tracer18Active  = .false.
      tracer19Active  = .false.
      tracer20Active = .false.
      tracer21Active  = .false.
      tracer22Active  = .false.
      tracer23Active  = .false.
      tracer24Active  = .false.
      tracer25Active  = .false.
      tracer26Active  = .false.
      tracer27Active  = .false.
      tracer28Active  = .false.
      tracer29Active  = .false.
      tracer30Active = .false.
      tracer31Active  = .false.
      tracer32Active  = .false.
      tracer33Active  = .false.

      if (config_num_passive_tracers >= 1) then
         call mpas_log_write('----- Setting up packages for passive tracers -----')
         call mpas_log_write('Activating tracer 1')
         tracer1Active = .true.
      if (config_num_passive_tracers >= 2) then
         call mpas_log_write('Activating tracer 2')
         tracer2Active = .true.
      if (config_num_passive_tracers >= 3) then
         call mpas_log_write('Activating tracer 3')
         tracer3Active = .true.
      if (config_num_passive_tracers >= 4) then
         call mpas_log_write('Activating tracer 4')
         tracer4Active = .true.
      if (config_num_passive_tracers >= 5) then
         call mpas_log_write('Activating tracer 5')
         tracer5Active = .true.
      if (config_num_passive_tracers >= 6) then
         call mpas_log_write('Activating tracer 6')
         tracer6Active = .true.
      if (config_num_passive_tracers >= 7) then
         call mpas_log_write('Activating tracer 7')
         tracer7Active = .true.
      if (config_num_passive_tracers >= 8) then
         call mpas_log_write('Activating tracer 8')
         tracer8Active = .true.
      if (config_num_passive_tracers >= 9) then
         call mpas_log_write('Activating tracer 9')
         tracer9Active = .true.
      if (config_num_passive_tracers >= 10) then
         call mpas_log_write('Activating tracer 10')
         tracer10Active = .true.
      if (config_num_passive_tracers >= 11) then
         call mpas_log_write('Activating tracer 11')
         tracer11Active = .true.
      if (config_num_passive_tracers >= 12) then
         call mpas_log_write('Activating tracer 12')
         tracer12Active = .true.
      if (config_num_passive_tracers >= 13) then
         call mpas_log_write('Activating tracer 13')
         tracer13Active = .true.
      if (config_num_passive_tracers >= 14) then
         call mpas_log_write('Activating tracer 14')
         tracer14Active = .true.
      if (config_num_passive_tracers >= 15) then
         call mpas_log_write('Activating tracer 15')
         tracer15Active = .true.
      if (config_num_passive_tracers >= 16) then
         call mpas_log_write('Activating tracer 16')
         tracer16Active = .true.
      if (config_num_passive_tracers >= 17) then
         call mpas_log_write('Activating tracer 17')
         tracer17Active = .true.
      if (config_num_passive_tracers >= 18) then
         call mpas_log_write('Activating tracer 18')
         tracer18Active = .true.
      if (config_num_passive_tracers >= 19) then
         call mpas_log_write('Activating tracer 19')
         tracer19Active = .true.
      if (config_num_passive_tracers >= 20) then
         call mpas_log_write('Activating tracer 20')
         tracer20Active = .true.
      if (config_num_passive_tracers >= 21) then
         call mpas_log_write('Activating tracer 21')
         tracer21Active = .true.
      if (config_num_passive_tracers >= 22) then
         call mpas_log_write('Activating tracer 22')
         tracer22Active = .true.
      if (config_num_passive_tracers >= 23) then
         call mpas_log_write('Activating tracer 23')
         tracer23Active = .true.
      if (config_num_passive_tracers >= 24) then
         call mpas_log_write('Activating tracer 24')
         tracer24Active = .true.
      if (config_num_passive_tracers >= 25) then
         call mpas_log_write('Activating tracer 25')
         tracer25Active = .true.
      if (config_num_passive_tracers >= 26) then
         call mpas_log_write('Activating tracer 26')
         tracer26Active = .true.
      if (config_num_passive_tracers >= 27) then
         call mpas_log_write('Activating tracer 27')
         tracer27Active = .true.
      if (config_num_passive_tracers >= 28) then
         call mpas_log_write('Activating tracer 28')
         tracer28Active = .true.
      if (config_num_passive_tracers >= 29) then
         call mpas_log_write('Activating tracer 29')
         tracer29Active = .true.
      if (config_num_passive_tracers >= 30) then
         call mpas_log_write('Activating tracer 30')
         tracer30Active = .true.
      if (config_num_passive_tracers >= 31) then
         call mpas_log_write('Activating tracer 31')
         tracer31Active = .true.
      if (config_num_passive_tracers >= 32) then
         call mpas_log_write('Activating tracer 32')
         tracer32Active = .true.
      if (config_num_passive_tracers >= 33) then
         call mpas_log_write('Activating tracer 33')
         tracer33Active = .true.
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if 
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if      

      !
      ! Limited-area
      !
      nullify(config_apply_lbcs)
      call mpas_pool_get_config(configs, 'config_apply_lbcs', config_apply_lbcs)

      nullify(limited_areaActive)
      call mpas_pool_get_package(packages, 'limited_areaActive', limited_areaActive)

      if (config_apply_lbcs) then
         limited_areaActive = .true.
      else
         limited_areaActive = .false.
      end if

      !
      ! JEDI data assimilation
      !
      nullify(config_jedi_da)
      call mpas_pool_get_config(configs, 'config_jedi_da', config_jedi_da)

      nullify(jedi_daActive)
      call mpas_pool_get_package(packages, 'jedi_daActive', jedi_daActive)

      if (associated(config_jedi_da) .and. associated(jedi_daActive)) then
         jedi_daActive = config_jedi_da
      else
         ierr = ierr + 1
         call mpas_log_write('Package setup failed for ''jedi_da''. '// &
              'Either ''jedi_da'' is not a package, or ''config_jedi_da'' is not a namelist option.', &
              messageType=MPAS_LOG_ERR)
      end if

#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at
      !least one physics parameterization is called (using the logical moist_physics):
      call physics_namelist_check(configs)

      local_ierr = atmphys_setup_packages(configs, packages, iocontext)
      if (local_ierr /= 0) then
         ierr = ierr + 1
         call mpas_log_write('Package setup failed for atmphys in core_atmosphere', messageType=MPAS_LOG_ERR)
      end if
#endif

   end function atm_setup_packages


   !***********************************************************************
   !
   !  function atm_setup_clock
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
   function atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_derived_types, only : MPAS_Clock_type, mpas_pool_type
      use atm_core, only : atm_simulation_clock_init

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      ierr = 0

      call atm_simulation_clock_init(core_clock, configs, ierr)

   end function atm_setup_clock


   !***********************************************************************
   !
   !  function atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function atm_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types, only : mpas_log_type, domain_type, mpas_coupler_type
      use mpas_log, only : mpas_log_init, mpas_log_open

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo    !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain          !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err
      type (mpas_coupler_type), pointer :: mpas_cpl

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      !
      ! Make log filenames specific to the role of each MPI task
      !
      mpas_cpl => domain % mpas_cpl
      write(logInfo % outputLog % filename, '(a,i2.2,a,i4.4,a)') 'log.atmosphere.role', mpas_cpl % role(), '.', domain % dminfo % my_proc_id, '.out'
      write(logInfo % errorLog % filename, '(a,i2.2,a,i4.4,a)') 'log.atmosphere.role', mpas_cpl % role(), '.', domain % dminfo % my_proc_id, '.err'

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call mpas_log_write('')
#ifdef SINGLE_PRECISION
      call mpas_log_write('Using default single-precision reals')
#else
      call mpas_log_write('Using default double-precision reals')
#endif
      call mpas_log_write('')

   end function atm_setup_log!}}}


   !***********************************************************************
   !
   !  function atm_get_mesh_stream
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
   function atm_get_mesh_stream(configs, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         call mpas_log_write('config_do_restart was not found when defining mesh stream.', messageType=MPAS_LOG_ERR)
         ierr = 1
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end function atm_get_mesh_stream


   !***********************************************************************
   !
   !  function atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function atm_setup_decompositions


   !***********************************************************************
   !
   !  function atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be
   !>  an atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type
      use mpas_pool_routines, only : mpas_pool_get_config
      use mpas_log, only : mpas_log_write

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      integer, pointer :: cam_pcnst
      integer :: err_level

      ierr = 0

      call atm_generate_structs(block, block % structs, block % dimensions, block % packages)

      !
      ! When MPAS-A is operating as a dycore in CAM, the scalars/scalars_tend var_arrays are
      ! allocated by the call to atm_allocate_scalars, below. The CAM-MPAS interface layer
      ! should have added a config, cam_pcnst, to the configs pool to indicate how many scalars
      ! are to be allocated.
      !
      nullify(cam_pcnst)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, 'cam_pcnst', cam_pcnst)
      call mpas_pool_set_error_level(err_level)
      if (associated(cam_pcnst)) then
         call mpas_log_write('')
         call mpas_log_write('** Config ''cam_pcnst'' is defined with a value of $i', intArgs=[cam_pcnst])
         call mpas_log_write('   Scalars will be allocated separately from Registry-defined variables')
         call mpas_log_write('')
         ierr = atm_allocate_scalars(block, cam_pcnst)
      end if

   end function atm_setup_block


   !***********************************************************************
   !
   !  function atm_allocate_scalars
   !
   !> \brief   Allocate scalars and scalars_tend var_arrays
   !> \author  Michael G. Duda
   !> \date    20 May 2020
   !> \details 
   !>  When MPAS-A is operating as a dycore for CAM, the scalars and
   !>  scalars_tend var_arrays are allocated separately from other Registry-
   !>  defined variables, since the set of scalars to be handled by the dycore
   !>  is not known until runtime. This routine allocates these var_arrays,
   !>  but it does not define which constituent is at which position in
   !>  var_arrays; this is defined later in the CAM-MPAS interface layer.
   !
   !-----------------------------------------------------------------------
   function atm_allocate_scalars(block, num_scalars) result(ierr)

      use mpas_derived_types, only : block_type

      use mpas_derived_types, only : mpas_pool_type, field3dReal, MPAS_LOG_ERR
      use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_add_dimension, mpas_pool_add_field
      use mpas_log, only : mpas_log_write

      implicit none

      ! Arguments
      type (block_type), pointer :: block
      integer, intent(in) :: num_scalars

      ! Return value
      integer :: ierr

      ! Local variables
      integer :: i, j, timeLevs
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tendPool
      type (field3dReal), dimension(:), pointer :: scalarsField


      ierr = 0

      !
      ! Allocate scalars var_array
      !
      nullify(statePool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)

      if (.not. associated(statePool)) then
         call mpas_log_write('No pool named ''state'' was found in atm_allocate_scalars', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      timeLevs = 2

      call mpas_pool_add_dimension(statePool, 'num_scalars', num_scalars)

      allocate(scalarsField(timeLevs))

      do i = 1, timeLevs
         scalarsField(i) % block => block
         scalarsField(i) % fieldName = 'scalars'
         scalarsField(i) % dimNames(1) = 'num_scalars'
         scalarsField(i) % dimNames(2) = 'nVertLevels'
         scalarsField(i) % dimNames(3) = 'nCells'
         scalarsField(i) % defaultValue = 0.0
         scalarsField(i) % missingValue = -1.0
         scalarsField(i) % isDecomposed = .true.
         scalarsField(i) % hasTimeDimension = .true.
         scalarsField(i) % isActive = .true.
         scalarsField(i) % isVarArray = .true.
         scalarsField(i) % isPersistent = .true.

         allocate(scalarsField(i) % constituentNames(num_scalars))

         allocate(scalarsField(i) % attLists(num_scalars))
         do j = 1, num_scalars
            allocate(scalarsField(i) % attLists(j) % attList)
         end do

      end do

      call mpas_pool_add_field(statePool, 'scalars', scalarsField)
      call mpas_pool_add_field(block % allFields, 'scalars', scalarsField)


      !
      ! Allocate scalars_tend var_array
      !
      nullify(tendPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      if (.not. associated(tendPool)) then
         call mpas_log_write('No pool named ''tend'' was found in atm_allocate_scalars', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      timeLevs = 1

      call mpas_pool_add_dimension(tendPool, 'num_scalars_tend', num_scalars)

      allocate(scalarsField(timeLevs))

      do i = 1, timeLevs
         scalarsField(i) % block => block
         scalarsField(i) % fieldName = 'scalars_tend'
         scalarsField(i) % dimNames(1) = 'num_scalars_tend'
         scalarsField(i) % dimNames(2) = 'nVertLevels'
         scalarsField(i) % dimNames(3) = 'nCells'
         scalarsField(i) % defaultValue = 0.0
         scalarsField(i) % missingValue = -1.0
         scalarsField(i) % isDecomposed = .true.
         scalarsField(i) % hasTimeDimension = .true.
         scalarsField(i) % isActive = .true.
         scalarsField(i) % isVarArray = .true.
         scalarsField(i) % isPersistent = .true.

         allocate(scalarsField(i) % constituentNames(num_scalars))

         allocate(scalarsField(i) % attLists(num_scalars))
         do j = 1, num_scalars
            allocate(scalarsField(i) % attLists(j) % attList)
         end do

      end do

      call mpas_pool_add_field(tendPool, 'scalars_tend', scalarsField)
      call mpas_pool_add_field(block % allFields, 'scalars_tend', scalarsField)

   end function atm_allocate_scalars

#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module atm_core_interface
