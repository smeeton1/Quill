module diswalk
use eigenSolve
use until
use func
use Densutil
use supperoperator
implicit none


 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp

 contains



 subroutine L_make_DCDG(SO,alpha)
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: H,D
  real(kdp)                                    :: sca

  n=16 !size(D,1)
  allocate(H(n,n),D(n,n))
  
  sca=1.0/sqrt(2.0)
  H(:,:)=cmplx(0,0)
  D(:,:)=cmplx(0,0)

  do i=1,n
    if(mod(i,2).eq.0)then
       H(i,i)=-1*sca
       H(i,i-1)=1*sca
    else
       H(i,i)=1*sca
       H(i,i+1)=1*sca
    endif
  enddo
  
  do i=1,n
    if(mod(i,2).eq.0)then
       D(i,i)=-1*sca
    else
       D(i,i+1)=1*sca
    endif
  enddo


 
  call B_SO_for_H((1-alpha)*H,SO)
   
  call Stand_to_SO(alpha*D,SO)
 

  deallocate(H,D)

 end subroutine
 
 !makes the transition operator for a 8 node line
 subroutine L_make_DTDG(SO,alpha)
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: H, D
  real(kdp)                                    :: sca

  n=16!size(D,1)
  allocate(H(n,n),D(n,n))
  
  sca=1.0
  H(:,:)=cmplx(0,0)
  D(:,:)=cmplx(0,0)

  
  H(1,15)=1
  do i=2,n
    if(mod(i,2).eq.0)then
       H(i,i+2)=1*sca
    else
       H(i,i-2)=1*sca
    endif
  enddo
  
  do i=1,n
    if(mod(i,2).eq.0)then
       D(i,i+2)=1*sca
    endif
  enddo
  call write_Mat_real('H',H)
 
  call B_SO_for_H((1-alpha)*H,SO)    
  call Stand_to_SO(alpha*D,SO)

  deallocate(H,D)

 end subroutine 
  
! function for TC to geather only T has direction in it 
 subroutine L_make_DDsF(SO,alpha)
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: Hc,Ht,D
  real(kdp)                                    :: sca

  n=16 !size(D,1)
  allocate(Hc(n,n),D(n,n),Ht(n,n))
  
  sca=1.0/sqrt(2.0)
  Hc(:,:)=cmplx(0,0)
  Ht(:,:)=cmplx(0,0)
  D(:,:)=cmplx(0,0)

  do i=1,n
    if(mod(i,2).eq.0)then
       H(i,i)=-1*sca
       H(i,i-1)=1*sca
    else
       H(i,i)=1*sca
       H(i,i+1)=1*sca
    endif
  enddo

  Ht(1,15)=1
  do i=2,n
    if(mod(i,2).eq.0)then
       Ht(i,i+2)=1*sca
    else
       Ht(i,i-2)=1*sca
    endif
  enddo
  
  do i=1,n
    if(mod(i,2).eq.0)then
       D(i,i+2)=1*sca
    endif
  enddo

  Ht=matmul(Ht,Hc)
  D=matmul(D,Hc)
 
  call B_SO_for_H((1-alpha)*Ht,SO)
   
  call Stand_to_SO(alpha*D,SO)
 

  deallocate(H,D)

 end subroutine
  
   
 !Performs a directional walk for a discreet walker on an 8 node line
 subroutine Dis_Dir_8line(rho, t, alpha, beta, filename, r)
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha,beta,err
  logical, intent(in)                          :: r
  character(len=80),intent(in)                 :: filename
  real(kdp)                                    :: error
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SOC, SOT, psi
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
  
  n=16
  allocate(SOC(n*n,n*n),SOT(n*n,n*n),rho_out(n*n),psi(t,n))
  
  SOC(:,:)=cmplx(0,0)
  SOT(:,:)=cmplx(0,0)
  call L_make_DDsF(SOC,alpha)
  !call L_make_DTDG(SOT,beta)
  
 ! call write_Mat('SO',SO)
  call extract_pointerS(rho, psi(1,:))
 ! write(*,*)size(rho),size(SOC,1),size(SOT,1)
  i=1
  error=1
  rho_out=rho
  do i=0,t
    rho_out= matmul(SOC,rho_out)
    call extract_pointerS(rho_out, psi(i+1,1:n))
  enddo
  if(r)then
    call write_Mat_real(filename,psi)
  else
    call write_Mat(filename,psi)
  endif
  
  !call write_moments(filename,psi)
  

  call write_Vec(filename,psi(2,1:n))
  norm= get_Norm(rho_out)
  open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
  write(4,*)''
  write(4,"(a,2F7.3)")' norm= ',norm
  close(4)
  rho=rho_out
  
  deallocate(SOC,SOT,rho_out,psi)

 end subroutine
 
 
 end module