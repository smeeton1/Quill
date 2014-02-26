module class
use until
use solver
use Densutil
use func
implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
 
 subroutine pagerank(D, p, t, alpha, filename, v, r)
   complex(kdp), dimension(:,:), intent(in)    :: D
  complex(kdp), dimension(:), intent(inout)    :: p
  real(kdp), intent(in)                        :: alpha
  integer, intent(in)                          :: t
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm, moment1, moment2
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
 end subroutine
 
 
 end module