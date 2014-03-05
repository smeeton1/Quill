module class
use until
use solver
use Densutil
use func
use eigenSolve
implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
 
 subroutine pagerank_it(D, p, t, alpha, filename, v, r)
   complex(kdp), dimension(:,:), intent(in)    :: D
  complex(kdp), dimension(:), intent(inout)    :: p
  real(kdp), intent(in)                        :: alpha
  integer, intent(in)                          :: t
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm, moment1, moment2
  complex(kdp), dimension(:), allocatable      :: rho_out ,alpha2 
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(rho_out(n),alphs2(n)) 
  rho_out = p
  
  alpha2(:)=cmplx((1.0-alpha),0.0)
  
  rho_out = alpha*matmul(D,rho_out) + alphs2
  
  
 end subroutine
 
 
 
  subroutine pagerank_ei(D, p, alpha, filename, work)
   complex(kdp), dimension(:,:), intent(in)    :: D
  complex(kdp), dimension(:), intent(inout)    :: p
  real(kdp), intent(in)                        :: alpha
  logical, intent(out)                         :: work
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n,j
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: eig  
  
  n=size(D,1)
  allocate(SO(n,n),psi(n,n),eig(n))
  SO=D
  do i=1,n
    do j=1,n
      D(i,j)=alpha*D(i,j)+(1.0-alpha)/real(n)
    enddo
  enddo
  
  call eigenVecR(SO,eig,psi)
  j=0
  do i=1,n
   if (1-eig(i)).lt.(0.000001)then
    j=i
    goto 99
   endif
  enddo
  
  99 continue
  if (j).ne.0 then
   call write_Vec(filename,psi(:,j)) 
   p(:)=psi(:,j)
   norm=Sum(psi(:,j))
   open(4,file=filename,ACCESS='append',ACTION='write')
   write(4,*)'Norm= ',norm
   close(4)   
   work=.true.
  else
   open(4,file=filename,ACCESS='append',ACTION='write')
   write(4,*)'NO eigenvalue equal to 1'
   close(4)
   work=.false.
  endif
  
  
 end subroutine
 
 end module