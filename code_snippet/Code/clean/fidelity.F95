!*********************************************
!Finds "fidelity" of two states
!
!
!Can be compiled with or without BLAS
!
!
!*********************************************
function fidelity(state_one, state_two)
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Precision of numerical values
  integer,          parameter :: dp=selected_real_kind(15, 300)	!IEEE 754 Double Precision
    
  !Function variables
  integer :: n                    									!integer size of state vectors
  real(kind=dp) :: fidelity       									!Fidelity variable for function output
  complex(kind=dp) :: ZDOTC										!ZDOT variable for blas library function

  !Input state vectors
  complex(kind=dp), dimension(:), Allocatable :: state_one		!First input state vector
  complex(kind=dp), dimension(:), Allocatable :: state_two		!Second input state vector

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! ARRAY SIZE CHECK !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *, SIZE(state_one), SIZE(state_two)

  !Checks to see if array sizes of two state vectors match
  !If they do assigns value to n for inner product calculation
  !Otherwise stops program and prints error message
  if(SIZE(state_one).eq.SIZE(state_two)) then
    n = SIZE(state_one)
  else
    stop 'Array size mismatch in fidelity function'
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !If BLAS is enabled when compiling (-lblas), this code segment
  !will be compiled into final program using BLAS functions
#ifdef lblas
    
  fidelity = (abs(ZDOTC(n, state_one,1, state_two, 1)))**2
    
  !If BLAS is not enabled when compiling, this section of intrinsic functions
  !will be used instead.
#else
    
  fidelity = (abs(DOT_PRODUCT(state_one, state_two)))**2
    
!ends the preprocessor if statement
#endif
    
  !returns value of fidelity from subroutine to main program
  return
 
end function fidelity
  
