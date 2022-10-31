module pfasst_main
 use pfasst            !< contains all neccessary mods for pfasst rountines. pfasst/src/pfasst.f90
 use pf_mod_mpi
 use pf_mod_comm_mpi
 use pf_mod_comm

  use, intrinsic :: iso_c_binding
  implicit none

contains

  !!! still working on it
  subroutine Pf_Main(AmrIceHolderPtr,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,numGridPointsBisicles,pf_plot_prefix,PF_VERBOSE) bind(c, name="Pf_Main")
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

    type(c_ptr),value :: AmrIceHolderPtr
    integer(c_int), value :: crse_nsteps
    real(c_double), value :: dt_bisicles
    real(c_double), value :: Tfin_bisicles
    integer(c_int), value :: maxStep_bisicles
    integer(c_int), value :: numGridPointsBisicles
    character(kind=c_char), intent(IN) :: pf_plot_prefix(:)
    logical(c_bool), value :: PF_VERBOSE

    !> local vars
    type(pf_pfasst_t), target      :: pf           !<  the main pfasst structure
    type(pf_comm_t)                :: comm         !<  the communicator (here it is mpi)
    type(bisicles_vector_encap), allocatable:: y_0       !<  the initial condition
    type(bisicles_vector_encap), allocatable:: y_end     !<  the solution at the final time
    class(pf_encap_t), allocatable :: y_0_base
    class(pf_encap_t), allocatable :: y_end_base
    type(bisicles_vector_factory)           :: bvf !< bisicle_vector
    character(256)                 :: pf_fname     !<  file name for input of PFASST parameters
    class(pf_level_t),  pointer    :: lev          !<  Level to set up

    integer, allocatable           :: lev_shape(:,:)
    integer                        ::  l   !  loop variable over levels
    integer                        :: level_index, nnodes
    integer ::  ierror
    integer :: nproc, rank, error
    integer :: space_comm, time_comm, space_color, time_color

    real(c_double) :: norm
    real(pfdp), pointer :: v(:)
    type(my_sweeper_t) :: s
    real(pfdp) :: t, y_0_IC, y_end_IC
    type(c_ptr) :: z_c_ptr

    !> first initialize mpi
    !> Initialize MPI
    !call mpi_init(ierror)
    !if (ierror /= 0) &
    !   stop "ERROR: Can't initialize MPI."


    ! check size
    nproc = 2
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)


    !> Read problem parameters
    call probin_init(pf_fname)

    call create_simple_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)

    !>  Set up communicaton
    call pf_mpi_create(comm, time_comm)
    print *, '-------------- done assigning mpi communicator from bisicles to pfasst ---------------'

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname) !< o
    print *, '--------------------------- done creating pfasst obj ---------------------------------'


    !> ----- initialize the vectors & solvers for later, not initial condition of setting values -----
     call PfasstBisiclesInit(pf, lev_shape,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,numGridPointsBisicles, AmrIceHolderPtr)
    !> PfasstBisiclesInit=allocate lev_shape+ulevel+factory+sweeper+level_set_size+pf_setup
    !> ----- end of initialize the vectors & solvers for later, including params change -----
    !print *, 'pfasst_main.f90 0000 grid size from bisicles to pfasstpf%levels(2)%lev_shape ',pf%levels(2)%lev_shape 
    print *, 'check if temporal params are passed in correctly:'
    print *, 'pfasst_main.f90 0000 T final ',Tfin
    print *, 'pfasst_main.f90 0000 dt ',dt
    print *, 'pfasst_main.f90 0000 nsteps ',nsteps

    !>  Output run parameters
    if (PF_VERBOSE) then
      call print_loc_options(pf)
      ! check if num of time steps assigned correctly

      print *, '--------------------------- done print out ------------------------------------------'
    end if

    
    !>  Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    print *, '--------------------------- done pf adding hooks -------------------------------------'
   

    !> setup initial condition
    level_index = 1
    !> Set the initial condition
    call bvf%create_single(y_0_base, level_index, lev_shape(level_index,:))
    !print *,'y_end create single '
    call bvf%create_single(y_end_base, level_index, lev_shape(level_index,:))
    y_0 = cast_as_bisicles_vector(y_0_base)
    y_end = cast_as_bisicles_vector(y_end_base)

    !> intialize y_0 and y_end to some value
    y_0_IC = 1.0_pfdp
    y_end_IC = 10.7_pfdp
    !call BisiclesAssignIC(pf,y_0,y_end,y_0_IC,y_end_IC)
    !print *,'y_end assign single '
    call y_end%setval(1000.0_pfdp)
    call y_0%setval(1000.0_pfdp)
    !call y_0%copy(y_end)
    !print *,'y_end axpy single '
    !call y_end%axpy(-1.0_pfdp, y_0)
    !print *,'y_end done axpy single '
    !call y_0%pack(v)

    !call y_end%eprint()
    !norm = y_end%norm()
    t=0.1
    z_c_ptr=pf%cptr_AmrIceHolder

    !call y_end%savesnap()
    !call ABORT

    s = cast_as_my_sweeper_t(pf%levels(level_index)%ulevel%sweeper)
    !call s%f_eval(y_0, t, level_index, y_end, z_c_ptr)
    !call Initialize_H(y_0) ! NEED TO BE FIXED
    
    print *, '--------------------------- done setting up IC --------------------------------------'

    !> run pfasst
    call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps,y_end)

    !> save the results
    call pfasst_bisicles_save_results(pf, AmrIceHolderPtr)

    
    !> Close mpi
    !call mpi_finalize(ierror)

  end subroutine Pf_Main


  !subroutine Pf_Run() bind(c, name="Pf_Run")


  !end subroutine Pf_Run





end module pfasst_main