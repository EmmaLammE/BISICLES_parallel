module pfasst_main
  use pfasst            !< contains all neccessary mods for pfasst rountines. pfasst/src/pfasst.f90
  use pf_mod_mpi
  use pf_mod_comm_mpi
  use pf_mod_comm
 
   use, intrinsic :: iso_c_binding
   implicit none
 
 contains
 
   !!! still working on it
   subroutine Pf_Main(AmrIceHolderPtr,pf_comm_fromBisicles,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,evolve_velocity_bisicles,numGridPointsBisicles,pf_num_repeats,PF_VERBOSE) bind(c, name="Pf_Main")
     use pfasst             !< contains all neccessary mods for pfasst rountines. pfasst/src/pfasst.f90
     use pf_mod_dtype
     use pf_mod_mpi
     use pf_mod_comm_mpi
     use pf_mod_comm
     use probin            !< Local module reading/parsing problem parameters
     use hooks             !< Local module for diagnostics and i/o
     use pf_my_sweeper     !< Local module for sweeper
     use pf_my_level       !< Local module for level
     use encap
     use pf_space_comm
     use pfasst_bisicles
     use mpi
 
     type(c_ptr),value :: AmrIceHolderPtr
     ! type(c_ptr),value :: pf_commPtr
     integer(c_int), value :: pf_comm_fromBisicles
     integer(c_int), value :: crse_nsteps
     real(c_double), value :: dt_bisicles
     real(c_double), value :: Tfin_bisicles
     integer(c_int), value :: maxStep_bisicles
     integer(c_int), value :: numGridPointsBisicles
     integer(c_int), value :: pf_num_repeats
    !  character(kind=c_char), intent(IN) :: pf_plot_prefix(:)
     logical(c_bool), value :: PF_VERBOSE
     logical(c_bool), value :: evolve_velocity_bisicles
 
     !> local vars
     type(pf_pfasst_t), target      :: pf           !<  the main pfasst structure
     type(pf_comm_t)                :: comm         !<  the communicator (here it is mpi)
     type(bisicles_vector_encap), allocatable:: y_0       !<  the initial condition
     type(bisicles_vector_encap), allocatable:: y_end     !<  the solution at the final time
     type(bisicles_vector_encap), allocatable:: y_test     !<  test solution 
     class(pf_encap_t), allocatable :: y_0_base
     class(pf_encap_t), allocatable :: y_end_base
     class(pf_encap_t), allocatable :: y_test_base
     type(bisicles_vector_factory)           :: bvf !< bisicle_vector
     character(256)                 :: pf_fname     !<  file name for input of PFASST parameters
     class(pf_level_t),  pointer    :: lev          !<  Level to set up
 
     integer, allocatable           :: lev_shape(:,:)
     integer                        ::  l   !  loop variable over levels
     integer                        :: level_index, nnodes
     integer ::  ierror
     integer :: nproc, rank, error
     integer :: space_comm, time_comm, space_color, time_color
     type(pf_comm_t) :: pf_time_comm      !<  the communicator (here it is mpi)
     integer :: f_comm
 
     real(c_double) :: norm
     real(pfdp), pointer :: v(:)
     type(my_sweeper_t) :: s
     real(pfdp) :: t, y_0_IC, y_end_IC
     type(c_ptr) :: z_c_ptr
     real :: pf_start_time, pf_finish_time, elapse_time, elapse_time_max,elapse_time_min,elapse_time_sum, elapse_time_mean
     real :: elapse_time_per_repeat, elapse_time_per_repeat_max, elapse_time_per_repeat_min, elapse_time_per_repeat_sum
     integer :: i, num_repeat_pf_run

    ! setup IO stripe width
    !  MPI_Info :: 
     integer :: ierr, info

    ! Create a new info object
    call MPI_Info_create(info, ierr)
    ! Set the stripe width hint
    call MPI_Info_set(info, "striping_factor", "24", ierr)
 
 
 
     !> Read problem parameters
     call probin_init(pf_fname)
 
     call create_simple_communicators(nspace, ntime, pf_comm_fromBisicles, space_comm, time_comm, space_color, time_color, space_dim)
     ! f_comm = MPI_Comm_c2f(pf_commPtr)
     ! print *,'MPI_COMM_WORLD ',MPI_COMM_WORLD
     ! print *,'time communicator in pfasst ',pf_comm_fromBisicles,', pf mpi time comm ',time_comm,', pf mpi space comm ',nspace
     call mpi_comm_size(pf_comm_fromBisicles, nproc, error)
     call mpi_comm_rank(pf_comm_fromBisicles, rank,  error)
     print *, 'time communicator in pfasst, num of procs for one time step ',nproc, ', rank ', rank
 
     !>  Set up communicaton
     call pf_mpi_create(comm, pf_comm_fromBisicles)
     ! print *, '-------------- done assigning mpi communicator from bisicles to pfasst ---------------'
     ! print *, "after assigning mpi comm, comm%nproc ",comm%nproc
     !>  Create the pfasst structure
     call pf_pfasst_create(pf, comm, fname=pf_fname) !< o
     ! print *, "after create pf obj, pf%comm%nproc ",pf%comm%nproc
     ! print *, '--------------------------- done creating pfasst obj ---------------------------------'
 
     !> ----- initialize the vectors & solvers for later, not initial condition of setting values -----
      call PfasstBisiclesInit(pf, lev_shape,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,evolve_velocity_bisicles,numGridPointsBisicles, AmrIceHolderPtr)
     !  print *, 'after pfasst init'
     !  call pf%levels(pf%state%finest_level)%Q(1)%eprint()
      !> PfasstBisiclesInit=allocate lev_shape+ulevel+factory+sweeper+level_set_size+pf_setup
     !> ----- end of initialize the vectors & solvers for later, including params change -----
     !print *, 'pfasst_main.f90 0000 grid size from bisicles to pfasstpf%levels(2)%lev_shape ',pf%levels(2)%lev_shape 
     ! print *, 'check if temporal params are passed in correctly:'
     ! print *, 'pfasst_main.f90 0000 T final ',Tfin
     ! print *, 'pfasst_main.f90 0000 dt ',dt
     ! print *, 'pfasst_main.f90 0000 nsteps ',nsteps
     !  print *,'5 pf%comm%nproc ', pf%comm%nproc
     !>  Output run parameters
    !  if (PF_VERBOSE) then
    !    call print_loc_options(pf, pf_comm_fromBisicles)
    !    ! check if num of time steps assigned correctly
 
    !    ! print *, '--------------------------- done print out ------------------------------------------'
    !  end if
 
     
     !>  Add some hooks for output
     call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
     ! print *, '--------------------------- done pf adding hooks -------------------------------------'
     ! print *, 'after adding hooks'
     ! call pf%levels(pf%state%finest_level)%Q(1)%eprint()
     ! print *,'6 pf%comm%nproc ', pf%comm%nproc
     !> setup initial condition
     level_index = 1
     !> Set the initial condition
     call bvf%create_single(y_0_base, level_index, lev_shape(level_index,:))
     !print *,'y_end create single '
     call bvf%create_single(y_end_base, level_index, lev_shape(level_index,:))
     call bvf%create_single(y_test_base, level_index, lev_shape(level_index,:))
 
     y_0 = cast_as_bisicles_vector(y_0_base)
     y_end = cast_as_bisicles_vector(y_end_base)
     y_test = cast_as_bisicles_vector(y_test_base)
 
     !> intialize y_0 and y_end to some value
     y_0_IC = 1.0_pfdp
     y_end_IC = 10.7_pfdp
     !call BisiclesAssignIC(pf,y_0,y_end,y_0_IC,y_end_IC)
     !print *,'y_end assign single '
     call y_end%setval(1000.0_pfdp)
     call y_0%setval(1000.0_pfdp)
    !  print *, "Bisicles Disjoint Boxes in current PFASST temporal processor:"
    !  call y_0%eprintLevelDataBox(y_0)
     ! packing test
     ! call y_0%eprint()
     ! allocate(v(num_grid_points))
     ! call y_0%pack(v)
     ! call y_test%unpack(v)
     ! call y_test%eprint()
     ! call exit(0)
 
     !call y_0%copy(y_end)
     !print *,'y_end axpy single '
     !call y_end%axpy(-1.0_pfdp, y_0)
     !print *,'y_end done axpy single '
     !call y_0%pack(v)
 
     !call y_end%eprint()
     !norm = y_end%norm()
    !  t=0.1
     z_c_ptr=pf%cptr_AmrIceHolder
 
     !call y_end%savesnap()
     !call ABORT
 
     s = cast_as_my_sweeper_t(pf%levels(level_index)%ulevel%sweeper)
     !call s%f_eval(y_0, t, level_index, y_end, z_c_ptr)
     !call Initialize_H(y_0) ! NEED TO BE FIXED
     
     print *, '--------------------------- Start pfasst run --------------------------------------'
     ! call PfasstBisiclesPrintAmr(y_0,pf%cptr_AmrIceHolder)
     !> run pfasst
     print *,"Start pfasst timing... for ",pf_num_repeats," repeats"
     num_repeat_pf_run = pf_num_repeats
     elapse_time_per_repeat_sum = 0
     elapse_time_per_repeat_max = 0
     elapse_time_per_repeat_min = 1e8
     do i = 1,num_repeat_pf_run
      print *,"  pf repeats # ",i
       call cpu_time(pf_start_time)
       !> pfasst run
       call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps,y_end)
       !> end pfasst run
       call cpu_time(pf_finish_time)
       elapse_time_per_repeat = pf_finish_time-pf_start_time
       elapse_time_per_repeat_sum = elapse_time_per_repeat_sum + elapse_time_per_repeat
       if(elapse_time_per_repeat>elapse_time_per_repeat_max) then
         elapse_time_per_repeat_max = elapse_time_per_repeat
       end if
       if (elapse_time_per_repeat<elapse_time_per_repeat_min) then
         elapse_time_per_repeat_min = elapse_time_per_repeat
       end if
     enddo
     elapse_time = elapse_time_per_repeat_sum/num_repeat_pf_run
    !  print '("Averaged time for pfasst run on rank "I0" = ",f8.3," seconds for "I0" repeats.")',rank,elapse_time,num_repeat_pf_run
 
     !> calculate average time of pf run
     call MPI_REDUCE(elapse_time_per_repeat_max,elapse_time_max,1,MPI_REAL,mpi_max,0,pf_comm_fromBisicles,ierror) ! compute mpi global sum and put to root processor 0
     call MPI_REDUCE(elapse_time_per_repeat_min,elapse_time_min,1,MPI_REAL,mpi_min,0,pf_comm_fromBisicles,ierror) ! compute mpi global sum and put to root processor 0
     call MPI_REDUCE(elapse_time,elapse_time_sum,1,MPI_REAL,mpi_sum,0,pf_comm_fromBisicles,ierror) ! compute mpi global sum and put to root processor 0
     if(pf%rank==0) then
       print '("Global average run time = ",f8.3,", max = ",f8.3,", min = ",f8.3," seconds.")',elapse_time_sum/nproc,elapse_time_max,elapse_time_min
     end if
 
     !> save the results
     call pfasst_bisicles_save_results(pf, AmrIceHolderPtr)
 
     
     !> Close mpi
    !  call mpi_finalize(ierror)
 
   end subroutine Pf_Main
 
 
   !subroutine Pf_Run() bind(c, name="Pf_Run")
 
 
   !end subroutine Pf_Run
 
 
 
 
 
 end module pfasst_main