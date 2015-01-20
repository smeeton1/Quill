module Densutil

!some function for working with density matirices.
!all functions are assuming that the density matrix rho is in the vectorised form of lenght n*n and that the colomes follow each other in order.

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains

 !gets the norm of the density matrix Tr(rho)
 complex(kdp) function get_Norm(rho) result(norm)
  complex(kdp),dimension(:),intent(in):: rho
  integer                             :: n, sqn, i, j
  
  n=size(rho)
  sqn=int(sqrt(real(n)))

  norm=rho(1)
  j=1
  do i=2,n
    if(i.eq.(j+sqn+1))then  
      norm=norm+rho(i)
      j=i
    endif
  enddo


 end function

 !writes out the density matrix
 subroutine write_rho(rho, filename)
  complex(kdp),dimension(:),intent(in):: rho
  character(len=80),intent(in)        :: filename 
  integer                             :: n, i, j

  n=int(sqrt(real(size(rho))))

  open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
  write(4,*)''
  write(4,*)'Density Matrix Rho'
  do i=1,n
   do j=0,n-1
       write(4,"(F7.3,a,a,F7.3,a)",ADVANCE='no')real(rho(i+j*n)),' ','i',aimag(rho(i+j*n)),' '
   enddo 
   write(4,*)''
  enddo
  close(4)

 end subroutine

 !gets the pointer states from the density matrix
 !assumes these are the states along the diagonal
 subroutine extract_pointerS(rho, point_state)
  complex(kdp),dimension(:),intent(in)   :: rho
  complex(kdp),dimension(:),intent(inout):: point_state
  integer                                :: n, sqn, i, j,k
  
  n=size(rho)
  sqn=int(sqrt(real(n)))
  
  point_state(1)=rho(1)
  k=1
  j=1
  do i=2,n
    if(i.eq.(j+sqn+1))then 
      k=k+1 
      point_state(k)=rho(i)
      j=i
    endif
  enddo
  
 end subroutine
  
 complex(kdp) function Norm_vec(A) result(norm)
  complex(kdp),dimension(:),intent(in)       :: A
  integer				     :: i,n
  
  n=size(A)
  norm=cmplx(0.0,0.0)
  do i=1,n
   norm =norm+conjg(A(i))*A(i)
  enddo
 
 end function
 
 real(kdp) function vec_norm_real(A) result(norm)
  complex(kdp),dimension(:),intent(inout)    :: A
  integer				     :: j,n
  n=size(A)

  norm=0.0
  do j=1,n
    norm = norm + real(conjg(A(j))*A(j))
  enddo
  
 end function




end module