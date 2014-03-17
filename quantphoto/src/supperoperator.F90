!code to build the Superoperator SO from the hamiltonian H and using damping terms. 
! Both function require that SO, H, D are squar matrixs and SO is n^2 in size withn being the size of both H and D

module supperoperator
use until
use solver
use Densutil
use func
implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
!make SO for the Hamiltonian only 
subroutine B_SO_for_H(H,SO)
 complex(kdp), dimension(:,:), intent(in)     :: H
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n


 n=size(H,1)
 SO(:,:)=cmplx(0,0)
 

!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1 
   do j=0,n-1
    do k=1,n
     do l=1,n
        if(i.eq.j)then
	    if(k.eq.l)then
             SO(i*n+k,j*n+l)=cmplx(0,1)*(H(i+1,i+1)-H(k,k))
            elseif(k.lt.l)then
             SO(i*n+k,j*n+l)=-cmplx(0,1)*H(l,k)
            else
	     SO(i*n+k,j*n+l)=-cmplx(0,1)*H(l,k)
	    endif 
        elseif((i.ne.j).and.(k.eq.l))then
            if((i.ne.1).and.(j.eq.1))then
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            elseif((i.eq.1).and.(j.ne.1))then
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            else
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            endif
        endif
     enddo
    enddo
   enddo
 enddo



end subroutine


!adding in dephasing for the operator style <O|rho|O>
subroutine A_D_to_SO(D,SO)
 complex(kdp), dimension(:,:), intent(in)     :: D
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n


 n=size(D,1)
!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1
   do j=0,n-1
    do k=1,n
     do l=1,n
      if((i.eq.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,i+1)
      elseif((i.lt.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,j+1)
      elseif((i.gt.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,j+1)
      endif
     enddo
    enddo
   enddo
 enddo

end subroutine


! adds the dephasing for when it has the form: Ot rho O -1/2{Ot O , rho}
subroutine Stand_to_SO(D,SO)
 complex(kdp), dimension(:,:), intent(in)     :: D
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,q,n


 n=size(D,1)
!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1
   do j=0,n-1
    do k=1,n
     do l=1,n
      if((i.eq.j).and.(k.eq.l).and.((i+1).eq.k))then
        do  q=1,n
          if(q.ne.(i+1))then
	    SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,q)
          endif
        enddo
      elseif((i.eq.j).and.(k.eq.l).and.((i+1).ne.k))then
        do  q=1,n
	    SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-0.5*D(i+1,q)-0.5*D(k,q)
        enddo
      elseif((i.lt.j).and.(k.eq.(i+1)).and.(l.eq.(j+1)))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)+D(l,k)
      elseif((i.gt.j).and.(k.eq.(i+1)).and.(l.eq.(j+1)))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)+D(l,k)
      endif
     enddo
    enddo
   enddo
 enddo

end subroutine


! routine to take a transition matrix and make the full super operator
! this will generate the H matrix so it dose not need to be put in
 subroutine L_make_DG(SO,D,alpha)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:,:), intent(inout)  :: SO
  real(kdp),intent(in)                         :: alpha
  integer                                      :: i,j,k,l,n
  complex(kdp), dimension(:,:), allocatable    :: H

  n=size(D,1)
  allocate(H(n,n))
!       call write_Mat('h',D)
   
  H(:,:)=cmplx(0,0)

  do i=1,n
    do j=1,n
     if(D(i,j).ne.0)then
       H(i,j)=1
       H(j,i)=1
     end if
    enddo
  enddo

 
  call B_SO_for_H((1-alpha)*H,SO)
  
!     call write_Mat('h',H)
!     call write_Mat('h',SO)
    
    
  call Stand_to_SO(alpha*D,SO)
  
!     call write_Mat('h',SO)  

  deallocate(H)

 end subroutine
 
 
 ! performs a walk on a directed graph using open system methods and dose the walk from time 1 to t
 ! extract the pointer states for each step.
 ! if v is true more details are wrten to the file 
 ! if r is true only writes out the real part of the pointer states.
 subroutine Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
  complex(kdp), dimension(:,:), intent(in)     :: D
  complex(kdp), dimension(:), intent(inout)    :: rho
  real(kdp), intent(in)                        :: alpha
  integer, intent(in)                          :: t
  logical, intent(in)                          :: v, r
  character(len=80),intent(in)                 :: filename
  integer                                      :: i,n
  complex(kdp), dimension(:,:), allocatable    :: SO, psi
  complex(kdp)                                 :: norm, moment1, moment2
  complex(kdp), dimension(:), allocatable      :: rho_out  
  character(len=90)                            :: filename2
  
  n=size(D,1)
  allocate(SO(n*n,n*n),rho_out(n*n),psi(t,n))
  
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
  call L_make_DG(SO,D,alpha)
 ! call write_Mat('SO',SO)

  do i=1,t
    call expm(SO,real(i,kdp),rho,rho_out)
    call extract_pointerS(rho_out, psi(i,1:n))
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
    call write_Vec(filename,psi(t,1:n))
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
  rho=rho_out
  
  deallocate(SO,rho_out,psi)

 end subroutine


 subroutine Dir_Gra_Con(D, rho, err, alpha, filename, r)
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
  call L_make_DG(SO,D,alpha)
 ! call write_Mat('SO',SO)
  call extract_pointerS(rho, psi(1,:))
  i=1
  error=1
  rho_out(:)=cmplx(0,0)
  do while (error.gt.err)
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