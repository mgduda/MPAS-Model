! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core

   use mpas_framework
   use mpas_timekeeping

   use test_core_halo_exch
   use test_core_field_tests
   use test_core_timekeeping_tests
   use mpas_stream_manager

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  function test_core_init
   !
   !> \brief   MPAS Core Initialization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs the necessary initialization of a core. This can
   !>   involve things like setting up coeffiecients for remapping and advection,
   !>   but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   function test_core_init(domain, startTimeStamp) result(iErr)!{{{
   
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      type (mpas_pool_type), pointer :: modelPool
      character (len=StrKIND), pointer :: xtime
      integer :: iErr

      !
      ! For testing alarms
      !
      integer :: istep
      type (MPAS_Clock_type) :: test_clock
      type (MPAS_Time_type) :: test_startTime
      type (MPAS_Time_type) :: test_stopTime
      type (MPAS_Time_type) :: test_currTime
      type (MPAS_Time_type) :: test_alarmTime
      type (MPAS_Time_type) :: test_refTime
      type (MPAS_TimeInterval_type) :: test_timeStep
      type (MPAS_TimeInterval_type) :: test_alarmTimeInterval


      iErr = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, iErr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

      call mpas_pool_get_subpool(domain % blocklist % structs, 'model', modelPool)
      call mpas_pool_get_array(modelPool, 'xtime', xtime)
      xtime = startTimeStamp

      call mpas_stream_mgr_read(domain % streamManager, ierr=iErr)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=iErr)

      !
      ! Alarm tests
      !
      call mpas_log_write('==== Begin alarm tests ====')

      call mpas_set_time(test_startTime, YYYY=2000, MM=01, DD=01, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_set_time(test_stopTime, YYYY=2100, MM=01, DD=01, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_set_timeInterval(test_timeStep, dt=3600.0_RKIND, ierr=ierr)

      call mpas_set_time(test_alarmTime, YYYY=2000, MM=01, DD=01, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_set_timeInterval(test_alarmTimeInterval, dt=86400.0_RKIND, ierr=ierr)

      call mpas_create_clock(test_clock, test_startTime, test_timeStep, test_stopTime, ierr=ierr)

      call mpas_add_clock_alarm(test_clock, 'foobar', test_alarmTime, test_alarmTimeInterval, ierr)

      do istep = 1, 24*365
          if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
              call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)
              test_currTime = mpas_get_clock_time(test_clock, MPAS_NOW, iErr)
              call mpas_get_time(test_currTime, dateTimeString=startTimeStamp) 
              call mpas_log_write('**ALARM** '//trim(startTimeStamp))
          end if
          call mpas_advance_clock(test_clock, ierr=ierr)
      end do

      test_currTime = mpas_get_clock_time(test_clock, MPAS_NOW, iErr)
      call mpas_get_time(test_currTime, dateTimeString=startTimeStamp) 

      call mpas_log_write('')
      call mpas_log_write('Now it is '//trim(startTimeStamp))
      call mpas_log_write('')

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2000-01-01_08')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2000, MM=01, DD=01, H=8, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2020-01-01_18')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2020, MM=01, DD=01, H=18, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2000-06-15_00')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2000, MM=06, DD=15, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2020-01-01_00')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2020, MM=01, DD=01, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_set_clock_direction(test_clock, MPAS_BACKWARD, ierr)

      call mpas_log_write('')
      call mpas_log_write('>>>> Running the clock backwards <<<<<')
      call mpas_log_write('')

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2000-01-01_08')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2000, MM=01, DD=01, H=8, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2020-01-01_18')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2020, MM=01, DD=01, H=18, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2000-06-15_00')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2000, MM=06, DD=15, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_log_write('=================================')
      call mpas_log_write('Setting ref time to 2020-01-01_00')
      call mpas_log_write('---------------------------------')
      call mpas_set_time(test_refTime, YYYY=2020, MM=01, DD=01, H=0, M=0, S=0, S_n=0, S_d=0, ierr=ierr)
      call mpas_adjust_alarm_to_reference_time(test_clock, 'foobar', test_refTime, ierr)
      call mpas_print_alarm(test_clock, 'foobar', ierr)
      if (mpas_is_alarm_ringing(test_clock, 'foobar', ierr=ierr)) then
          call mpas_log_write('-> Alarm is RINGING')
      else
          call mpas_log_write('-> Alarm is NOT ringing')
      end if
      call mpas_reset_clock_alarm(test_clock, 'foobar', ierr=ierr)

      call mpas_set_clock_direction(test_clock, MPAS_FORWARD, ierr)

      call mpas_log_write('==== End alarm tests ====')

   end function test_core_init!}}}


   !***********************************************************************
   !
   !  function test_core_run
   !
   !> \brief   MPAS Core Run Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   function test_core_run(domain) result(iErr)!{{{
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_vector_operations
      use mpas_geometry_utils
      use test_core_streams, only : test_core_streams_test
      use test_core_sorting, only : test_core_test_sorting
      use mpas_halo_testing, only : mpas_halo_tests
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: iErr

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: pool
      type (mpas_pool_iterator_type) :: itr
      integer :: numThreads, threadLimit, maxThreads
      integer, dimension(:), allocatable :: threadErrs

      maxThreads = mpas_threading_get_max_threads()
      allocate(threadErrs(maxThreads))
      threadErrs(:) = 0

#if 0
      !
      ! Test performance of framework sorting routines
      !
      call test_core_test_sorting(domain, iErr)
      if (iErr == 0) then
         call mpas_log_write(' * Sorting tests: SUCCESS')
      else
         call mpas_log_write(' * Sorting tests: FAILURE', MPAS_LOG_ERR)
      end if

      !
      ! Test functionality of mpas_halo module
      !
      call mpas_log_write('')
      call mpas_log_write('Testing mpas_halo module:')
      call mpas_halo_tests(domain, iErr)
      if (iErr == 0) then
         call mpas_log_write('* mpas_halo tests: SUCCESS')
      else
         call mpas_log_write('* mpas_halo tests: FAILURE', MPAS_LOG_ERR)
      end if
      call mpas_log_write('')

      iErr = 0

      call mpas_unit_test_fix_periodicity(iErr)
      call mpas_unit_test_triangle_signed_area_sphere(iErr)

      call mpas_unit_test_velocity_conversion(iErr)
      call mpas_unit_test_wachspress_hexagon(iErr)
      call mpas_unit_test_wachspress_triangle(iErr)

      !$omp parallel default(firstprivate) shared(domain, threadErrs)
      call test_core_halo_exch_test(domain, threadErrs, iErr)
      !$omp end parallel
      if ( iErr == 0 ) then
         call mpas_log_write(' * Halo Exchange Test: SUCCESS')
      else
         call mpas_log_write(' * Halo Exchange Test: FAILURE', MPAS_LOG_ERR)
      end if

      !$omp parallel default(firstprivate) shared(domain, threadErrs)
      call test_core_test_fields(domain, threadErrs, ierr)
       if ( iErr == 0 ) then
         call mpas_log_write(' * Field Tests: SUCCESS')
      else
         call mpas_log_write(' * Field Tests: FAILURE', MPAS_LOG_ERR)
      end if
      !$omp end parallel

      call test_core_streams_test(domain, threadErrs, iErr)
      if ( iErr == 0 ) then
         call mpas_log_write('Stream I/O tests: SUCCESS')
      else
         call mpas_log_write('Stream I/O tests: FAILURE', MPAS_LOG_ERR)
      end if


      call test_core_test_intervals(domain, threadErrs, iErr)

      ! Test writing of block write streams, which have the prefix 'block_'
!      block => domain % blocklist
!      do while (associated(block))
!         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID="block_.*")
!         call mpas_stream_mgr_block_write(domain % streamManager, block, streamID="block_.*", forceWriteNow=.true.)
!         block => block % next
!      end do
!
!      call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true.)
#endif

      !
      ! Test functionality of adjustments to alarm reference time
      !
      call mpas_log_write('')
      call mpas_log_write('Testing mpas_adjust_alarm_to_reference_time:')
      call mpas_adjust_alarm_tests(domain, iErr)
      if (iErr == 0) then
         call mpas_log_write('* mpas_adjust_alarm_tests tests - all tests passed: SUCCESS')
      else
         call mpas_log_write('* mpas_adjust_alarm_tests tests - $i failed tests: FAILURE', intArgs=[iErr])
      end if
      call mpas_log_write('')

      deallocate(threadErrs)

   end function test_core_run!}}}

   !***********************************************************************
   !
   !  function test_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   function test_core_finalize(domain) result(iErr)!{{{
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function test_core_finalize!}}}

end module test_core
