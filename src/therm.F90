module therm
use until
use solver
use Densutil
use func
use supperoperator
implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains

 
 !genteates the hamiltioan for a lattic of size m*p=n
 subroutine make_lattic(H,E,V,alpha,m)
  complex(kdp), dimension(:,:), intent(inout)  :: H
  complex(kdp), intent(in)                     :: E,V,alpha
  complex(kdp)                                 :: eval
  integer, intent(in)                          :: m
  integer                                      :: n, i,j,p
  
  n=size(H,1)
  H=cmplx(0.0,0.0)
  p=n/m
  eval=alpha*exp(E)
  
  H(1,1)=eval-V/(cmplx(p,0.0))
  H(1,2)=eval-V/(cmplx(p,0.0))
  H(1,m+1)=eval-V/(cmplx(p,0.0))
  
  H(n,n)=eval-V
  H(n,n-1)=eval-V
  H(n,n-m)=eval-V
  
  do i=2,n-1
     if(i.lt.m)then
       H(i,i)=eval-V/(cmplx(p,0.0))
       H(i,i+1)=eval-V/(cmplx(p,0.0))
       H(i,i-1)=eval-V/(cmplx(p,0.0))
       H(i,m+i)= eval-V/(cmplx(p,0.0))    
     elseif((i.ge.m).and.(i.le.(n-m)))then
       H(i,i)=eval-V/(cmplx(p,0.0)-cmplx(int(i/m),0.0))
       H(i,i+1)=eval-V/(cmplx(p,0.0)-cmplx(int(i/m),0.0))
       H(i,i-1)=eval-V/(cmplx(p,0.0)-cmplx(int(i/m),0.0))
       H(i,m+i)=eval-V/(cmplx(p,0.0)-cmplx(int(i/m),0.0))
       H(i,i-m)=eval-V/(cmplx(p,0.0)-cmplx(int(i/m),0.0))
     elseif(i.gt.(n-m))then
       H(i,i)=eval-V
       H(i,i+1)=eval-V
       H(i,i-1)=eval-V
       H(i,i-m)=eval-V
     endif
  enddo
  call write_Mat_real('H',H)
 
 end subroutine
 
 
 ! Makes the So for when there is tempuarture and colligions
 subroutine ele_L_Tcol(H,strength,T,coli,sink,conect,SO)
  complex(kdp), dimension(:,:), intent(in)     :: H
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  integer, intent(in)                          :: sink
  integer, dimension(:), intent(in)            :: conect
  complex(kdp), intent(in)                     :: strength, T, coli
  complex(kdp),dimension(:,:),allocatable      :: D
  integer                                      :: n
  
  n=size(H,1)
  allocate(D(n,n))
  
  call B_SO_for_H(H,SO)
  
  D(:,:)=T
  call Stand_to_SO(D,SO)
  
  D(:,:)=coli
  call De_en_cojn(D,SO)
  
  call De_sink(sink,conect,real(strength),SO)
  
  deallocate(D)
 end subroutine

 ! Makes the So for when there is tempuarture only
 subroutine ele_L_T(H,strength,T,sink,conect,SO)
  complex(kdp), dimension(:,:), intent(in)     :: H
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  integer, intent(in)                          :: sink
  integer, dimension(:), intent(in)            :: conect
  complex(kdp), intent(in)                     :: strength, T
  complex(kdp),dimension(:,:),allocatable      :: D
  integer                                      :: n
  
  n=size(H,1)
  allocate(D(n,n))
  
  call B_SO_for_H(H,SO)
  
  D(:,:)=T
  call Stand_to_SO(D,SO)
 
  call De_sink(sink,conect,real(strength),SO)
  
  deallocate(D) 
 
 end subroutine
 
 subroutine ele_tran(H,rho,Tend,dt,v,r,filename,strength,T,coli,sink,conect)
  complex(kdp), dimension(:,:), intent(in)     :: H
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: dt
  integer, intent(in)                          :: Tend, sink
  integer, dimension(:), intent(in)            :: conect
  complex(kdp), intent(in)                     :: strength, T, coli
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: psi, SO
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
  
  n=size(H,1)
  allocate(rho_out(n*n),psi(Tend+1,n),SO(n*n,n*n))
  
  if(coli.eq.0.0)then
    call ele_L_T(H,strength,T,sink,conect,SO)
  else
    call ele_L_Tcol(H,strength,T,coli,sink,conect,SO)
  endif
  !call write_Mat_real('SO',SO)
  open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
  if(v)then
    write(4,"(a)")'In the before'
    close(4)
    call write_rho(rho,filename)
    call extract_pointerS(rho, psi(1,1:n))
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a)")'Pointer States'
    close(4)
    call write_Vec(filename,psi(1,1:n))
    norm= get_Norm(rho)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a,2F7.3)")'norm= ',norm
    close(4)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a)")'In the after'
  endif
  close(4) 
 ! call write_Mat('SO',SO)
  call extract_pointerS(rho, psi(1,1:n))
 
  do i=1,Tend
    call expm(SO,dt*real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(i+1,1:n))
  enddo
  if(r)then
    call write_Mat_real(filename,psi)
  else
    call write_Mat(filename,psi)
  endif
  
  !call write_moments(filename,psi)
  
  if(v)then
    call write_rho(rho_out,filename)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a)")'Pointer States'
    close(4)
    call write_Vec(filename,psi(Tend,1:n))
    norm= get_Norm(rho_out)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a,2F7.3)")'norm= ',norm
    close(4)
  endif
    norm= get_Norm(rho_out)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a,2F7.3)")'norm= ',norm
    close(4)
 
  deallocate(rho_out,psi,SO)
 
 end subroutine


end module therm
