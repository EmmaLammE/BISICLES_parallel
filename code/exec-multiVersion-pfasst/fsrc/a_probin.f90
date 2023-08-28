!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem  y'=lam1*y+lam2*y
module probin
  use pfasst
  use pf_mod_mpi
  use omp_lib

  !  The namlist for local variables
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time
  integer, save :: nsteps    ! number of time steps
  integer, save :: num_grid_points, nspace, ntime, space_dim
  integer, save :: nproc_per_time ! num of processors for each time step, this is for bisicles spatial mpi

  character(len=128), save :: pfasst_nml  ! file for reading pfasst parameters
  type(c_ptr), save :: cptr_AmrIceHolder
  integer, save :: maxStep     ! num of max time step to run

  !namelist /params/  dt, Tfin, nsteps, num_grid_points,pfasst_nml,cptr_AmrIceHolder
  namelist /params/  pfasst_nml, nspace, ntime, nproc_per_time, space_dim

contains
  
  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname

    !  Some local variables for reading
    character(len=128) :: arg  !  command line argument
    character(len=128)    :: probin_fname   !<  file name for input parameters
    character(len=256) :: istring           ! stores command line argument
    character(len=1024) :: message          ! use for I/O error messages
    integer :: ios,iostat
    integer :: i   !  loop variable
    integer :: un  !  file read unit

    !> Set the name of the input file
    probin_fname = "pfasst_params.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)

    !> set defaults
    !nsteps  = -1
    nspace = 1
    ntime = 1 !nproc
    nproc_per_time = 1
    space_dim = 2

    !dt      = 0.01_pfdp
    !Tfin    = 1.0_pfdp
    pfasst_nml=probin_fname
    
    !>  Read in stuff from input file
    un = 9
    write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
    open(unit=un, file = probin_fname, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    !>  Read the command line
    i = 0
    do
       call get_command_argument(i, arg)
       if (LEN_TRIM(arg) == 0) EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          READ(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    !if (Tfin .gt. 0.0) dt = Tfin/dble(nsteps)

    !  Return the name of the file from which to read PFASST parameters
    pf_fname=pfasst_nml
  end subroutine probin_init

  !>  Subroutine to output run parameters 
  subroutine print_loc_options(pf,pf_comm, un_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer, intent(in) :: pf_comm
    integer,           intent(in   ), optional :: un_opt
    integer :: un = 6

    integer :: num_time_procs, rank_time_procs,num_global_procs, error

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt

    !>  Output the PFASST options with the LibPFASST routine
    call pf_print_options(pf,un_opt=un)

    call mpi_comm_rank(pf_comm, rank_time_procs,  error)
    call mpi_comm_size(pf_comm, num_time_procs,  error)
    call mpi_comm_size(MPI_COMM_WORLD, num_global_procs,  error)
    nspace = num_global_procs/num_time_procs

    !  Print out the local parameters
    write(un,*) '=================================================='
    write(un,*) ' '
    write(un,*) 'Local Variables'
    write(un,*) '----------------'
    !write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    !write(un,*) 'Dt:     ', Dt, '! Time step size'
    !write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'  
    write(un,*) 'num spatial procs per temporal proc:   ',   nspace
    write(un,*) 'num temporal procs:   ',   num_time_procs, ', rank ',rank_time_procs
    write(un,*) 'Number of spacial dimensions:   ', space_dim
    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  
  
  
end module probin
