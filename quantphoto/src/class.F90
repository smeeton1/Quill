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
 
 subroutine pagerank_it(D, p, err, alpha, filename)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: p
  real(kdp), intent(in)                        :: alpha, err
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  real(kdp)                                    :: error
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: rho_out, rho_out2 ,alpha2 
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(rho_out(n),rho_out2(n),alpha2(n)) 
  rho_out = p
  error=1
  alpha2(:)=cmplx((1.0-alpha),0.0)
  do while (error.gt.err)

   rho_out2 = alpha*matmul(D,rho_out) + alpha2
  
   error = maxval(abs(rho_out-rho_out2))

   rho_out = rho_out2
  enddo
  call write_Vec(filename,rho_out2) 
  norm=Sum(rho_out2)
  call write_Vec(filename,rho_out2/norm)
  open(4,file=filename,ACCESS='append',ACTION='write')
  write(4,*)' Norm= ',real(norm)
  close(4)   
!  p=rho_out2
  deallocate(rho_out,rho_out2,alpha2)
 end subroutine
 
 
 
 subroutine pagerank_ei(D, p, alpha, filename, work)
   complex(kdp), dimension(:,:), intent(in)     :: D
   complex(kdp), dimension(:), intent(inout)    :: p
   real(kdp), intent(in)                        :: alpha
   logical, intent(out)                         :: work
   character(len=80),intent(in)                 :: filename
   integer                                      :: i,n,j
   complex(kdp), dimension(:,:), allocatable    :: SO, psi
   complex(kdp)                                 :: norm
   complex(kdp), dimension(:), allocatable      :: eigs  
  
   n=size(D,1)
   allocate(SO(n,n),psi(n,n),eigs(n))
   SO=cmplx(0,0)
   do i=1,n
     do j=1,n
       SO(i,j)=alpha*D(i,j)+(1.0-alpha)/real(n)
     enddo
   enddo
   
   call eigenVecR(SO,eigs,psi)
   call write_Vec(filename,eigs)
   j=0
   do i=1,n
    if ( (1.0-real(eigs(i))) .lt. (0.000001)) then
     j=i
     goto 99
    endif
   enddo
   
   99 continue
   if (j .ne. 0 )then
    call write_Vec(filename,psi(:,j)) 
!    p(:)=psi(:,j)
    norm=Sum(psi(:,j)) 
    call write_Vec(filename,psi(:,j)/norm)
    open(4,file=filename,ACCESS='append',ACTION='write')
    write(4,*)' Norm= ',real(norm)
    close(4)   
    work=.true.
   else
    open(4,file=filename,ACCESS='append',ACTION='write')
    write(4,*)'NO eigenvalue equal to 1'
    close(4)
    work=.false.
   endif
  
   deallocate(psi,SO,eigs)
 end subroutine
 
 end module