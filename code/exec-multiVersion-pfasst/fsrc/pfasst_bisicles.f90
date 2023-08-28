
!>  Module for overwrite bisicles parameters f
module pfasst_bisicles

  use pfasst
  use pf_my_level
  use encap
  use probin
  use pf_my_sweeper

  use pf_mod_mpi
  use pf_mod_comm_mpi
  use pf_mod_comm
  implicit none


  !real(pfdp), save :: dt     ! time step
  !real(pfdp), save :: Tfin   ! Final time
  !integer, save :: nsteps    ! number of time steps
  !type(c_ptr), save :: cptr_AmrIceHolder

  !namelist /params/  dt, Tfin, nsteps
  !namelist /params/ cptr_AmrIceHolder


!type, extends(pf_pfasst_t) :: bisicles_holder_ptr
!      type(c_ptr) :: c_test_ptr = c_null_ptr ! c pointer
!end type bisicles_holder_ptr


contains

  subroutine PfasstBisiclesInit(pf, lev_shape,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,numGridPointsBisicles, AmrIceHolderPtr) 
    type(pf_pfasst_t), intent(inout) :: pf
    integer, allocatable, intent(inout) :: lev_shape(:,:)
    integer, intent(in)                        :: crse_nsteps
    real*8, intent(in)                        :: dt_bisicles
    real*8, intent(in)                        :: Tfin_bisicles
    integer, intent(in)                        :: maxStep_bisicles
    integer, intent(in)                        :: numGridPointsBisicles
    type(c_ptr), intent(in)                 :: AmrIceHolderPtr
    
    type(pf_comm_t) :: comm
    type(my_sweeper_t) :: sw_finest, sw_lev
    type(my_level_t) :: my_lev
    class(pf_level_t),  pointer :: lev

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: level_index

    integer :: nproc, rank, error
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: n_init, refine_factor, FComp_setup_flag
    logical :: setup_start_coarse_flag


    !> call mpi ??
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    !> first update the params from bisicles
    call update_bisicles_params(pf,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,numGridPointsBisicles, AmrIceHolderPtr)

    !> Loop over levels and set some level specific parameters
    allocate(lev_shape(pf%nlevels,1))
    ! print *, 'in pfasst bisicles init'
    ! call pf%levels(pf%state%finest_level)%Q(1)%eprint()
    do l = 1, pf%nlevels
       lev_shape(l,1) = num_grid_points
       !lev_shape(l,2) = AmrIceHolderPtr
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(bisicles_vector_factory::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !> setup sweeper solver??
       sw_lev = cast_as_my_sweeper_t(pf%levels(l)%ulevel%sweeper)
       call BisiclesSolverInit(sw_lev%c_bisicles_solver_ptr, &
                         l, &
                         lev_shape(l,1), &
                         cptr_AmrIceHolder)
      !  print *, 'after bisicles solver init'
      !  call pf%levels(pf%state%finest_level)%Q(1)%eprint()
       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,lev_shape(l,:)) ! shape_in=1, meaning grid size=1
       ! check if grid size is assigned correctly
       ! print *, 'pfasst_main.f90 0000 grid size from bisicles to pfasst ',pf%levels(l)%lev_shape
    end do
    !>  Set up some pfasst stuff
    ! level = 1 -> tauQ=8; level = 2 -> tauQ=16, no third level???
    ! print *, 'before pfasst set up'
    ! call pf%levels(pf%nlevels)%Q(1)%eprint()
    call pf_pfasst_setup(pf)
    ! print *, '--------------------------- done pf additional setup ---------------------------------'
    ! lev => pf%levels(2)
    ! print *, 'after pfasst set up, tauQ allocated ',allocated(pf%levels(2)%tauQ),'  print Q(1)'
    ! call pf%levels(pf%nlevels)%Q(1)%eprint()

  end subroutine PfasstBisiclesInit




  subroutine update_bisicles_params(pf,crse_nsteps,dt_bisicles,Tfin_bisicles,maxStep_bisicles,numGridPointsBisicles, AmrIceHolderPtr)
    type(pf_pfasst_t), intent(inout)           :: pf  
    integer, intent(in)                        :: crse_nsteps
    real*8, intent(in)                        :: dt_bisicles
    real*8, intent(in)                        :: Tfin_bisicles
    integer, intent(in)                        :: maxStep_bisicles
    integer, intent(in)                        :: numGridPointsBisicles
    type(c_ptr), intent(in)                 :: AmrIceHolderPtr

    type(bisicles_holder_ptr), pointer           :: pf_bisicles


    !  Some local variables for reading
    !print *,'pfasst_bisicles.f90 00000, crse_nsteps ',crse_nsteps
    !print *,'pfasst_bisicles.f90 00000, pf%state%nsteps ',pf%state%nsteps
    pf%state%nsteps = crse_nsteps
    nsteps = crse_nsteps
    dt = dt_bisicles
    Tfin = Tfin_bisicles
    maxStep = maxStep_bisicles

    num_grid_points=numGridPointsBisicles


      !  print *,'pfasst_bisicles.f90 pf '
       pf_bisicles => cast_as_pf_bisicles_t(pf)
!       pf_bisicles%c_test_ptr=AmrIceHolderPtr

    
    pf%cptr_AmrIceHolder=AmrIceHolderPtr
    cptr_AmrIceHolder=AmrIceHolderPtr
    !print *,'nsteps after changed ',nsteps
    !print *,'pfasst_bisicles.f90 pf%cptr_AmrIceHolder ', pf%cptr_AmrIceHolder
    !print *,'pfasst_bisicles.f90 pf%c_test_ptr ', pf_bisicles%c_test_ptr
    ! up till here, amr pointer is successfully assigned to pf_bisicles
    !call ABORT
  
  end subroutine update_bisicles_params


  subroutine pfasst_bisicles_save_results(pf, AmrIceHolderPtr)
    type(pf_pfasst_t), intent(in),target           :: pf 
    type(c_ptr), intent(in)                 :: AmrIceHolderPtr
    integer :: temp_level_to_save
    integer :: snapshot_to_save

    type(pf_level_t), pointer :: lev    !!  points to current level
    class(bisicles_vector_encap), pointer :: y_encap
    type(c_ptr)      :: c_AmrIceHolderPtr

    temp_level_to_save=pf%state%finest_level    ! 3 is the finest level
    snapshot_to_save=3 ! lev%nnodes is some huge num, need to be fixed

    lev => pf%levels(temp_level_to_save)   !  Assign level pointer
    c_AmrIceHolderPtr=pf%cptr_AmrIceHolder

    y_encap => cast_as_bisicles_vector(lev%Q(snapshot_to_save))
    !print *, ' saving results... y_encap ID y_encap%c_encap_ptr,',y_encap%c_encap_ptr
    print *, '                    temp_level_to_save,snapshot_to_save',temp_level_to_save,snapshot_to_save
    call y_encap%eprint()
    

    call y_encap%savesnap()
    !call PfasstBisiclesSaveResults(y_encap%c_encap_ptr,c_AmrIceHolderPtr)

  end subroutine pfasst_bisicles_save_results
  


  subroutine BisiclesAssignIC(pf,y_0,y_end,y_0_IC,y_end_IC)
    type(pf_pfasst_t), intent(inout), target   :: pf 
    type(bisicles_vector_encap), allocatable:: y_0, y_end
    real(pfdp) :: y_0_IC, y_end_IC

    !  local var
    class(pf_level_t), pointer :: lev
    integer                    :: level_index

    !>  first assign the value to y_0 and y_end bisicle vector
    call y_0%setval(y_0_IC)
    call y_end%setval(y_end_IC)

    !>  loop over levels
    !  pointer to finest  level to start

    if (pf%state%finest_level > 1) then ! fine level=3
       do level_index = pf%state%finest_level, 1, -1
          lev => pf%levels(level_index)
          call lev%q0%copy(y_0, flags=0) ! assign IC to level data q0 on each level
          call lev%ulevel%sweeper%spreadq0(pf,level_index, pf%state%t0) ! spread q0 to all Q in time
          !call lev%sweeper%spreadq0(pf,level_index, pf%state%t0) ! spread q0 to all Q in time


       end do  !  level_index = pf%state%finest_level, 2, -1
    end if
  end subroutine BisiclesAssignIC


  !function cast_as_pf_bisicles_t(level_polymorph) result(pf_bisicles_obj)
  !  class(pf_pfasst_t), intent(in), target :: level_polymorph
  !  type(bisicles_holder_ptr), pointer :: pf_bisicles_obj

  !  select type(level_polymorph)
  !  type is (bisicles_holder_ptr)
  !     pf_bisicles_obj => level_polymorph
  !  end select
  !end function cast_as_pf_bisicles_t
  
  

  subroutine PfasstBisiclesPrintAmr(y, AmrIceHolderPtr) 
    type(bisicles_vector_encap), allocatable:: y
    type(c_ptr), intent(in)          :: AmrIceHolderPtr

    call PfasstPrintAmr(y%c_encap_ptr,AmrIceHolderPtr)


  end subroutine PfasstBisiclesPrintAmr

end module pfasst_bisicles
