module minum
use eigenSolve
use until
use func
use Densutil
use supperoperator
implicit none


 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp

 contains
 
 ! routine to take a transition matrix and make the full super operator
! this will generate the H matrix so it dose not need to be put in
 subroutine L_make_D_min(SO,D,alpha)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: H
  real(kdp)                                    :: sca

  n=size(D,1)
  allocate(H(n,n))
  
  sca=1.0
  H(:,:)=cmplx(0,0)

  do i=1,n
    do j=1,n
     if(D(i,j).ne.0)then
       H(i,j)=1*sca
       H(j,i)=1*sca
     end if
    enddo
  enddo

 
  call B_SO_for_H((1-alpha)*H,SO)    
  call Stand_to_SO(alpha*D,SO)


  deallocate(H)

 end subroutine
 
  subroutine L_make_D_min_C(SO,D,alpha)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: H
  real(kdp)                                    :: sca

  n=size(D,1)
  allocate(H(n,n))
  
  sca=1.0
  H(:,:)=cmplx(0,0)

  do i=1,n
    do j=1,n
     if(D(i,j).ne.0)then
       H(i,j)=1*sca
       H(j,i)=1*sca
     end if
    enddo
  enddo

 
  call B_SO_for_H((1-alpha)*H,SO)    
  call De_en_cojn(alpha*D,SO)


  deallocate(H)

 end subroutine
 
 
 ! performs a walk on a directed graph using open system methods and dose the walk from time 1 to t
 ! extract the pointer states for each step.
 ! if v is true more details are wrten to the file 
 ! if r is true only writes out the real part of the pointer states.
 subroutine Min_Run(D, rho, t, alpha, filename, v, r)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha
  integer, intent(in)                          :: t
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm, moment1, moment2
  complex(kdp), dimension(:), allocatable      :: rho_out,ent  
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(SO(n*n,n*n),rho_out(n*n),psi(t+1,n),ent(t+1))
  write(filename2,'(a,a)')trim(filename),'-ent'
  
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
  
  
  SO(:,:)=cmplx(0,0)
  call L_make_D_min(SO,D,alpha)
  call extract_pointerS(rho, psi(1,1:n))
  ent(1)= VonNueE_vec(rho)
 
  do i=1,t
    call expm(SO,real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(i+1,1:n))
    ent(i+1)= abs(VonNueE_vec(rho_out))
  enddo
  
  if(r)then
    call write_Mat_real(filename,psi)
    call write_Vec_real(filename2,ent)
  else
    call write_Mat(filename,psi)
    call write_Vec_real(filename2,ent)
  endif
  
  
  if(v)then
    call write_rho(rho_out,filename)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a)")'Pointer States'
    close(4)
    call write_Vec(filename,psi(t,1:n))
    norm= get_Norm(rho_out)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a,2F7.3)")'norm= ',norm
    close(4)
  endif
  rho=rho_out
  
  deallocate(SO,rho_out,psi)

 end subroutine


 subroutine Min_Con(D, rho, err, alpha, filename, r)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha,err
  logical, intent(in)                          :: r
  character(len=80),intent(in)                 :: filename
  real(kdp)                                    :: error
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(SO(n*n,n*n),rho_out(n*n),psi(2,n))
  
  SO(:,:)=cmplx(0,0)
  call L_make_D_min(SO,D,alpha)
 ! call write_Mat('SO',SO)
  call extract_pointerS(rho, psi(1,:))
  i=1
  error=1
  rho_out(:)=cmplx(0,0)
  do while ((error.gt.err).or.(i.gt.10000))
    i=i+1
    call expm(SO,real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(2,1:n))
    error = maxval(abs(psi(1,1:n)-psi(2,1:n)))
    psi(1,:)=psi(2,:)
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
  
  deallocate(SO,rho_out,psi)

 end subroutine
 
  subroutine Min_Run(D, rho, t, alpha, filename, v, r)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha
  integer, intent(in)                          :: t
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm, moment1, moment2
  complex(kdp), dimension(:), allocatable      :: rho_out,ent  
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(SO(n*n,n*n),rho_out(n*n),psi(t+1,n),ent(t+1))
  write(filename2,'(a,a)')trim(filename),'-ent'
  
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
  
  
  SO(:,:)=cmplx(0,0)
  call L_make_D_min_C(SO,D,alpha)
  call extract_pointerS(rho, psi(1,1:n))
  ent(1)= VonNueE_vec(rho)
 
  do i=1,t
    call expm(SO,real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(i+1,1:n))
    ent(i+1)= abs(VonNueE_vec(rho_out))
  enddo
  
  if(r)then
    call write_Mat_real(filename,psi)
    call write_Vec_real(filename2,ent)
  else
    call write_Mat(filename,psi)
    call write_Vec_real(filename2,ent)
  endif
  
  
  if(v)then
    call write_rho(rho_out,filename)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a)")'Pointer States'
    close(4)
    call write_Vec(filename,psi(t,1:n))
    norm= get_Norm(rho_out)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)''
    write(4,"(a,2F7.3)")'norm= ',norm
    close(4)
  endif
  rho=rho_out
  
  deallocate(SO,rho_out,psi)

 end subroutine


 subroutine Min_Con(D, rho, err, alpha, filename, r)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha,err
  logical, intent(in)                          :: r
  character(len=80),intent(in)                 :: filename
  real(kdp)                                    :: error
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(SO(n*n,n*n),rho_out(n*n),psi(2,n))
  
  SO(:,:)=cmplx(0,0)
  call L_make_D_min_C(SO,D,alpha)
 ! call write_Mat('SO',SO)
  call extract_pointerS(rho, psi(1,:))
  i=1
  error=1
  rho_out(:)=cmplx(0,0)
  do while ((error.gt.err).or.(i.gt.10000))
    i=i+1
    call expm(SO,real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(2,1:n))
    error = maxval(abs(psi(1,1:n)-psi(2,1:n)))
    psi(1,:)=psi(2,:)
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
  
  deallocate(SO,rho_out,psi)

 end subroutine
 
 
 
 end module