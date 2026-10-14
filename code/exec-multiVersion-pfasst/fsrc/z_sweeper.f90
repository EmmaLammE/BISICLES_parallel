!
! This file is part of LIBPFASST.
!
module pf_my_sweeper
  use encap
  use pf_mod_imex_sweeper_bisicles
  use pfasst
  !use pf_my_level
  use probin
  implicit none
 
 
  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_bisicles_t) :: my_sweeper_t
   type(c_ptr) :: c_bisicles_solver_ptr= c_null_ptr
   !type(c_ptr) :: c_AmrIceHolderPtr= c_null_ptr
   contains
 
 
     procedure :: f_eval    !  Computes the explicit rhs terms
     procedure :: f_comp    !  Does implicit solves
 
 
     procedure :: initialize !  Overwrites imex sweeper initialize
     procedure :: destroy    !  Overwrites imex sweeper destroy
 
 
 end type my_sweeper_t
 
 
 
 
  interface
    subroutine BisiclesSolverInit(bisicles_solver_ptr, level_index,nx,c_AmrIceHolderPtr) bind(c, name="BisiclesSolverInit")
        use iso_c_binding
        type(c_ptr) :: bisicles_solver_ptr
        integer, value :: nx, level_index
        type(c_ptr), value :: c_AmrIceHolderPtr
     end subroutine BisiclesSolverInit
 
 
     subroutine BisiclesVectorSetHIC(y_0,c_AmrIceHolderPtr) bind(c, name="BisiclesVectorSetHIC")
        use iso_c_binding
        type(c_ptr), value :: y_0
        type(c_ptr), value :: c_AmrIceHolderPtr
     end subroutine BisiclesVectorSetHIC
 
 
 
 
    subroutine BisiclesSolverFEval(bisicles_dHdt,y, t, level_index, f, dt,maxStep,evolve_velocity, c_AmrIceHolderPtr) bind(c, name="BisiclesSolverFEval")
        use iso_c_binding
        type(c_ptr), value :: bisicles_dHdt
        type(c_ptr), value :: y
        real(c_double), value :: t
        integer, value :: level_index
        type(c_ptr), value :: f
        real(c_double), value :: dt
        integer(c_int), value :: maxStep
        logical(c_bool), value :: evolve_velocity
        type(c_ptr), value :: c_AmrIceHolderPtr
    end subroutine BisiclesSolverFEval
 
 
    subroutine PfasstPrintAmr(y,c_AmrIceHolderPtr) bind(c, name="PfasstPrintAmr")
      use iso_c_binding
      type(c_ptr), value :: y
      type(c_ptr), value :: c_AmrIceHolderPtr
   end subroutine PfasstPrintAmr
 
 
 
 
  end interface
 
 
 contains
 
 
  !>  Routine to set up sweeper variables and operators
  subroutine initialize(this, pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf
    integer, intent(in) :: level_index
    integer :: nx
    type(c_ptr) :: c_AmrIceHolderPtr
 
 
    !>  Call the imex sweeper initialization
    call this%imex_bisicles_initialize(pf,level_index)
 
 
    !>  Set variables for explicit and implicit parts (just to show you can)
    this%implicit=.FALSE.
    this%explicit=.TRUE.
 
 
    ! Space variables
    nx = pf%levels(level_index)%lev_shape(1)
    c_AmrIceHolderPtr=pf%cptr_AmrIceHolder
 
 
    !> Call the Bisicles solver initialization
    call BisiclesSolverInit(this%c_bisicles_solver_ptr, &
                         level_index, &
                         nx,c_AmrIceHolderPtr)
 
 
  end subroutine initialize
 
 
 
 
  subroutine Initialize_H(y_0)
    type(bisicles_vector_encap), intent(inout) :: y_0
    type(c_ptr) :: c_AmrIceHolderPtr
   
    call BisiclesVectorSetHIC(y_0%c_encap_ptr,c_AmrIceHolderPtr)
  end subroutine Initialize_H
 
 
  !>  destroy the sweeper type
  subroutine destroy(this, pf,level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout), target :: pf
    integer, intent(in) :: level_index
 
 
    !>  Call the imex sweeper destroy
    call this%imex_bisicles_destroy(pf,level_index)
 
 
    !  Nothing to do
 
 
  end subroutine destroy
 
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, c_AmrIceHolderPtr)
    !!!!!!!!!! commented out temporarily
    !use probin, only:  lam1, lam2
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    !integer,             intent(in)    :: piece
    type(c_ptr),         intent(in   ) :: c_AmrIceHolderPtr
    real(pfdp) :: val
   
    class(bisicles_vector_encap), pointer :: y_encap, f_encap
   
    !!!!!!!!!! commented out temporarily
    y_encap => cast_as_bisicles_vector(y)
    f_encap => cast_as_bisicles_vector(f)
 
 
    call BisiclesSolverFEval(this%c_bisicles_solver_ptr,y_encap%c_encap_ptr, t, level_index, f_encap%c_encap_ptr, dt, maxStep,evolve_velocity,c_AmrIceHolderPtr)
 
 
  end subroutine f_eval
 
 
  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    !!!!!!!!!! commented out temporarily
    !use probin, only:  lam1, lam2
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
 
 
    !class(bisicles_vector_encap), pointer :: y_encap, f_encap, rhs_encap
    real(pfdp) :: val
    ! EL - not implemented
  end subroutine f_comp
 
 
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !>  Here are some extra routines to help out
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial(y_0)
    type(bisicles_vector_encap), intent(inout) :: y_0
    real(pfdp) :: val
    !call exact(0.0_pfdp, val)
    !call y_0%setval(val)
  end subroutine initial
 
 
  !> Routine to return the exact solution
  subroutine exact(t, yex)
    !use probin, only: lam1,lam2
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex
 
 
    !yex=exp((lam1+lam2)*t)
 
 
  end subroutine exact
 
 
 
 
  function cast_as_my_sweeper_t(pf_sweeper_t_polymorph) result(my_sweeper_t_obj)
      class(pf_sweeper_t), intent(in), target :: pf_sweeper_t_polymorph
      type(my_sweeper_t), pointer :: my_sweeper_t_obj
 
 
      select type(pf_sweeper_t_polymorph)
      type is (my_sweeper_t)
         my_sweeper_t_obj => pf_sweeper_t_polymorph
      end select
    end function cast_as_my_sweeper_t
 
 
 end module pf_my_sweeper
 
 
 