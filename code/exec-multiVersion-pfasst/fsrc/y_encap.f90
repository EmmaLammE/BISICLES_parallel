!
! This file is part of LIBPFASST.
!

!> Simple bisicles_vector encapsulation example
!!

module encap
   use pfasst
   use probin
   implicit none
 
   !>  Type to create and destroy the local data encapsulation
   type, extends(pf_factory_t) :: bisicles_vector_factory
      type(c_ptr) :: c_factory_ptr = c_null_ptr ! c pointer
    contains
      procedure :: create_single  => bisicles_vector_create_single
      procedure :: create_array  => bisicles_vector_create_array
      procedure :: destroy_single => bisicles_vector_destroy_single
      procedure :: destroy_array => bisicles_vector_destroy_array
   end type bisicles_vector_factory
 
   !>  Type to extend the abstract encap and set procedure pointers
   type, extends(pf_encap_t) :: bisicles_vector_encap
      integer :: vector_size
      type(c_ptr) :: c_encap_ptr = c_null_ptr ! c pointer
    contains
      procedure :: setval => bisicles_vector_setval
      procedure :: copy => bisicles_vector_copy
      procedure :: norm => bisicles_vector_norm
      procedure :: pack => bisicles_vector_pack
      procedure :: unpack => bisicles_vector_unpack
      procedure :: axpy => bisicles_vector_axpy
      procedure :: eprint => bisicles_vector_eprint
      procedure :: eprintLevelDataBox => bisicles_vector_eprintLevelDataBox
      procedure :: getval
      procedure :: savesnap
   end type bisicles_vector_encap

   type, extends(pf_pfasst_t) :: bisicles_holder_ptr
      type(c_ptr) :: c_test_ptr = c_null_ptr ! c pointer
   end type bisicles_holder_ptr

   interface

      subroutine BisiclesVectorCreate(x, &
                                   num_grid_points, c_AmrIceHolderPtr) bind(c, name="BisiclesVectorCreate")
         use iso_c_binding
         type(c_ptr) :: x
         integer, value :: num_grid_points
         type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine BisiclesVectorCreate
    
      subroutine BisiclesVectorDestroy(x) bind(c, name="BisiclesVectorDestroy")
         use iso_c_binding
         type(c_ptr), value :: x
      end subroutine BisiclesVectorDestroy
    
      subroutine BisiclesVectorSetVal(x, val, c_AmrIceHolderPtr) bind(c, name="BisiclesVectorSetVal")
         use iso_c_binding
         type(c_ptr), value:: x
         real(c_double), value :: val
         type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine BisiclesVectorSetVal
    
      subroutine BisiclesVectorCopy(dest, src, c_AmrIceHolderPtr) bind(c, name="BisiclesVectorCopy")
         use iso_c_binding
         type(c_ptr), value :: dest, src
         type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine BisiclesVectorCopy
   
      function BisiclesCurrentVectorSize(x) result(vector_size) bind(c, name="BisiclesCurrentVectorSize")
         use iso_c_binding
         type(c_ptr), value :: x ! vector to be packed
         integer(c_int) :: vector_size
      end function
      
      function BisiclesVectorPack(x,c_AmrIceHolderPtr,level_id) result(z) bind(c, name="BisiclesVectorPack")
         use iso_c_binding
         type(c_ptr), value :: x ! vector to be packed
         type(c_ptr), value :: c_AmrIceHolderPtr
         type(c_ptr) :: z ! packed x is stored in z
         integer, value :: level_id
      end function
 
      subroutine BisiclesVectorUnpack(x, z, c_AmrIceHolderPtr) bind(c, name="BisiclesVectorUnpack")
         use iso_c_binding
         type(c_ptr), value :: x
         type(c_ptr), value :: z
         type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine BisiclesVectorUnpack
    
      function BisiclesVectorNorm(x, c_AmrIceHolderPtr) result(norm) bind(c, name="BisiclesVectorNorm")
        use iso_c_binding
        type(c_ptr), value :: x
        type(c_ptr), value :: c_AmrIceHolderPtr
        real(c_double) :: norm
      end function
    
      subroutine BisiclesVectorAxpy(y, a, x, c_AmrIceHolderPtr) bind(c, name="BisiclesVectorAxpy")
        use iso_c_binding
        type(c_ptr), value :: x, y
        real(c_double), value  :: a
        type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine BisiclesVectorAxpy
    
      subroutine BisiclesVectorL2Print(x) bind(c, name="BisiclesVectorL2Print")
        use iso_c_binding
        type(c_ptr), value :: x
      end subroutine BisiclesVectorL2Print

      subroutine BisiclesVectorLevelDataBox(x,y) bind(c, name="BisiclesVectorLevelDataBox")
         use iso_c_binding
         type(c_ptr), value :: x,y
       end subroutine BisiclesVectorLevelDataBox

      function BisiclesVectorGetVal(x, c_AmrIceHolderPtr) result(val) bind(c, name="BisiclesVectorGetVal")
        use iso_c_binding
        type(c_ptr), value :: x
        type(c_ptr), value :: c_AmrIceHolderPtr
        real(c_double) :: val
      end function

      subroutine PfasstBisiclesSaveResults(y, c_AmrIceHolderPtr) bind(c, name="PfasstBisiclesSaveResults")
          use iso_c_binding
          type(c_ptr), value :: y
          type(c_ptr), value :: c_AmrIceHolderPtr
      end subroutine PfasstBisiclesSaveResults

      integer function MPI_Comm_c2f(c_handle) bind(C, name="PfasstBisicles_MPI_Comm_c2f")
         use iso_c_binding
         type(c_ptr), value :: c_handle
      end function

   end interface

contains

  !>  The following are the base subroutines that encapsulation factories need to provide
  
  !>  Subroutine to allocate one encap
  subroutine bisicles_vector_create_single(this, x, level_index, lev_shape)
    class(bisicles_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: ierr

    integer :: num_grid_points
    real(pfdp) :: intialization_factor
    intialization_factor = -1.0
    
    !print *, 'cptr_AmrIceHolder', cptr_AmrIceHolder

    num_grid_points = lev_shape(1)
    !cptr_AmrIceHolder = lev_shape(2)

    allocate(bisicles_vector_encap::x, stat=ierr)
    !if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)

    select type(x)
    type is (bisicles_vector_encap)
      !  print *,'encap.f90 0000 calling BisiclesVectorCreate......................, num grid point ', num_grid_points
       call BisiclesVectorCreate(x%c_encap_ptr,num_grid_points, cptr_AmrIceHolder)
       !call BisiclesVectorAxpy(x%c_encap_ptr, intialization_factor, x%c_encap_ptr, cptr_AmrIceHolder)
       !print *, 'encap.f90 1111 finish BisiclesVectorCreate......................', x%c_encap_ptr
       x%vector_size = num_grid_points
    end select
  end subroutine bisicles_vector_create_single

  !> Subroutine to create an array of encaps
  subroutine bisicles_vector_create_array(this, x, n, level_index,lev_shape)
    class(bisicles_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:) ! array to be created
    integer,               intent(in   )              :: n  ! size of array to build
    integer,               intent(in   ) ::  level_index ! passed by default,  not needed here
    integer,               intent(in   ) ::  lev_shape(:) ! passed by default, not needed here
    integer :: i, ierr

    integer :: num_grid_points
    real(pfdp) :: intialization_factor
    intialization_factor = -1.0

    num_grid_points = lev_shape(1)
    !cptr_AmrIceHolder = lev_shape(2)

    allocate(bisicles_vector_encap::x(n),stat=ierr)
    !if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',n)
    select type(x)
    type is (bisicles_vector_encap)
       do i = 1, n
           call BisiclesVectorCreate(x(i)%c_encap_ptr,num_grid_points, cptr_AmrIceHolder)
           !call BisiclesVectorAxpy(x(i)%c_encap_ptr, intialization_factor, x(i)%c_encap_ptr, cptr_AmrIceHolder)
           x%vector_size = num_grid_points
       end do
    end select
  end subroutine bisicles_vector_create_array

  !> Subroutine to destroy a single array encap
  subroutine bisicles_vector_destroy_single(this, x)
    class(bisicles_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x
    integer                                           :: ierr
   
    select type(x)
    type is (bisicles_vector_encap)
       call BisiclesVectorDestroy(x%c_encap_ptr)
    end select
    deallocate(x,stat=ierr)
  end subroutine bisicles_vector_destroy_single

  !> Subroutine to destroy an array of arrays
  subroutine bisicles_vector_destroy_array(this, x)
    class(bisicles_vector_factory), intent(inout)              :: this
    class(pf_encap_t),     intent(inout), allocatable :: x(:)
    integer                                           :: i, ierr
    class(bisicles_vector_encap), pointer :: x_ptr

    select type(x)
    type is (bisicles_vector_encap)
       do i = 1, size(x)
          call BisiclesVectorDestroy(x(i)%c_encap_ptr)
       end do
    end select
    deallocate(x,stat=ierr)
  end subroutine bisicles_vector_destroy_array

  !>  The following are the base subroutines that all encapsulations must provide

  !> Subroutine to set array to a bisicles_vector  value.
  subroutine bisicles_vector_setval(this, val, flags)
    class(bisicles_vector_encap), intent(inout)      :: this
    real(c_double),     intent(in   )       :: val
    integer,        intent(in   ), optional :: flags
    !print *, '---------------- encap.f90 bisicles_vector_setval --------------------'
    !print *,'val ',val
    call BisiclesVectorSetVal(this%c_encap_ptr, val, cptr_AmrIceHolder)
    !print *, '---------------- bisicles_vector_setval done --------------------'
  end subroutine bisicles_vector_setval

  !> Subroutine to copy an array
  subroutine bisicles_vector_copy(this, src, flags)
    class(bisicles_vector_encap),    intent(inout)      :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    class(bisicles_vector_encap), pointer               :: src_bisicles_vector_encap

    select type(src)
    type is (bisicles_vector_encap)
       call BisiclesVectorCopy(this%c_encap_ptr, src%c_encap_ptr, cptr_AmrIceHolder)
    end select
  end subroutine bisicles_vector_copy

  !> Subroutine to pack into a flat array for sending
  subroutine bisicles_vector_pack(this, z, flags)
     class(bisicles_vector_encap), intent(in) :: this
    real(pfdp), intent(out) :: z(:)
    integer, intent(in), optional :: flags
    real(pfdp), pointer :: z_ptr(:)
    type(c_ptr) :: z_ptr_test
    type(c_ptr) :: z_c_ptr
    integer(c_int) :: num_grid_points
    integer :: level_id

   !  do level_id = 1, 1
    ! EL - packing should be correct now, without amr
      z_c_ptr = BisiclesVectorPack(this%c_encap_ptr,cptr_AmrIceHolder,flags)
      num_grid_points = BisiclesCurrentVectorSize(this%c_encap_ptr)
      call c_f_pointer(z_c_ptr, z_ptr, [num_grid_points]) ! convert z_c_ptr to z_ptr
   !  end do
      ! print *,'in packing num_grid_points ',num_grid_points,", size of z ",size(z_ptr)
      
      ! do i = 1,1024
      ! print *,'size of z_ptr ',i,z_ptr(i)
      ! end do
    z = z_ptr
   !  print *, size(z)
  end subroutine bisicles_vector_pack

  !> Subroutine to unpack  after receiving
  subroutine bisicles_vector_unpack(this, z, flags)
     class(bisicles_vector_encap), intent(inout) :: this
     real(pfdp),     intent(in   ) :: z(:)
     integer,     intent(in   ), optional :: flags
     real(pfdp), target :: z2(size(z))
     type(c_ptr) :: z_c_ptr
      ! print *, "before unpacking z in ",size(z)
     z2 = z
     z_c_ptr = c_loc(z2(1))
     
     call BisiclesVectorUnpack(this%c_encap_ptr, z_c_ptr, cptr_AmrIceHolder);
  end subroutine bisicles_vector_unpack

  !> Subroutine to define the norm of the array (here the abs value)
  function bisicles_vector_norm(this, flags) result (norm)
    class(bisicles_vector_encap), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(c_double) :: norm
    ! call BisiclesVectorNorm(this%c_encap_ptr, norm)
   !  print *,'     ... in bisicles vector norm'
    norm = BisiclesVectorNorm(this%c_encap_ptr, cptr_AmrIceHolder)
  end function bisicles_vector_norm

  !> Subroutine to compute y = a x + y where a is a bisicles_vector and x and y are arrays
  subroutine bisicles_vector_axpy(this, a, x, flags)
    class(bisicles_vector_encap), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(c_double), intent(in) :: a
    integer, intent(in), optional :: flags

    class(pf_encap_t), allocatable:: y_0, y_0_base

    !y_0 = cast_as_bisicles_vector(y_0_base)
    !y_0=x
    !call y_0%eprint()

    select type(x)
    type is (bisicles_vector_encap) 
       !print *, 'before calling BisiclesVectorAxpy this ptr ',this%c_encap_ptr, 'x ptr ',x%c_encap_ptr
       call BisiclesVectorAxpy(this%c_encap_ptr, a, x%c_encap_ptr, cptr_AmrIceHolder)
    end select
  end subroutine bisicles_vector_axpy

  !> Jordi stopped here
  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine bisicles_vector_eprint(this,flags)
    class(bisicles_vector_encap), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Print the  value
    call BisiclesVectorL2Print(this%c_encap_ptr)
  end subroutine bisicles_vector_eprint

  subroutine bisicles_vector_eprintLevelDataBox(this,x,flags)
   class(bisicles_vector_encap), intent(inout) :: this
   class(pf_encap_t), intent(in) :: x ! array to be printed
   integer,           intent(in   ), optional :: flags
   !  Print the  value
   select type(x)
    type is (bisicles_vector_encap)
      call BisiclesVectorLevelDataBox(this%c_encap_ptr,x%c_encap_ptr)
   end select
  end subroutine bisicles_vector_eprintLevelDataBox

  function getval(this) result(val)
     class(bisicles_vector_encap), intent(inout) :: this
     real(pfdp) :: val
     val = BisiclesVectorGetVal(this%c_encap_ptr, cptr_AmrIceHolder)
  end function

  !  Helper function to cast an abstract encap to the bisicles_vector_encap
  function cast_as_bisicles_vector(encap_polymorph) result(bisicles_vector_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(bisicles_vector_encap), pointer :: bisicles_vector_obj
    
    select type(encap_polymorph)
    type is (bisicles_vector_encap)
       bisicles_vector_obj => encap_polymorph
    end select
  end function cast_as_bisicles_vector


  !> Subroutine to save the snapshot 
  subroutine savesnap(this) 
    class(bisicles_vector_encap), intent(inout) :: this
    call PfasstBisiclesSaveResults(this%c_encap_ptr, cptr_AmrIceHolder)
  end subroutine savesnap

  function cast_as_pf_bisicles_t(level_polymorph) result(pf_bisicles_obj)
    class(pf_pfasst_t), intent(in), target :: level_polymorph
    type(bisicles_holder_ptr), pointer :: pf_bisicles_obj

    select type(level_polymorph)
    type is (bisicles_holder_ptr)
       pf_bisicles_obj => level_polymorph
    end select
  end function cast_as_pf_bisicles_t

end module encap
