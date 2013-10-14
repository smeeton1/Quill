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

  open(4,file=filename,ACCESS='append',ACTION='write')
  
  do i=1,n
   do j=0,n-1
       write(4,"(F7.3,a,a,F7.3,a)")real(rho(i+j*n)),' ','i',aimag(rho(i+j*n)),' '
   enddo 
   write(4,*)'\n'
  enddo
  close(4)

 end subroutine

 !gets the pointer states from the density matrix
 !assumes these are the states along the diagonal
 subroutine extract_pointerS(rho, point_state)
  complex(kdp),dimension(:),intent(in)   :: rho
  complex(kdp),dimension(:),intent(inout):: point_state
  integer                                :: n, sqn, i, j
  
  n=size(rho)
  sqn=int(sqrt(real(n)))

  norm=rho(1)

  do i=2,n
    if(i.eq.(j+sqn+1))then  
      point_state(i)=rho(i)
      j=i
    endif
  enddo


 end subroutine

end module