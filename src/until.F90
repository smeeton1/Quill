module until
 use supperoperator
 implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains



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
  call Stand_to_SO(alpha*D,SO)

  deallocate(H)

 end subroutine


! routine to read in a matrix from a file.
 subroutine Read_Mat(Filename,mat)
  complex(kdp),dimension(:,:),intent(in):: mat
  character(len=80),intent(in)          :: filename 
 
 end subroutine

 subroutine write_Mat(filename,mat)
  complex(kdp),dimension(:,:),intent(in):: mat
  character(len=80),intent(in)          :: filename 
  integer                               :: n, i, j

  n=size(mat)

  open(4,file=filename,ACCESS='append',ACTION='write')
  
  do i=1,n
   do j=1,n
       write(4,"(F7.3,a,a,F7.3,a)")real(mat(i,j)),' ','i',aimag(mat(i,j)),' '
   enddo 
   write(4,*)'\n'
  enddo
  close(4)
 
 end subroutine

 subroutine write_Vec(filename,vec)
  complex(kdp),dimension(:),intent(in)  :: vec
  character(len=80),intent(in)          :: filename 
  integer                               :: n, i

  n=size(vec)

  open(4,file=filename,ACCESS='append',ACTION='write')

  do i=1,n
    write(4,"(F7.3,a,a,F7.3,a)")real(vec(i)),' ','i',aimag(vec(i)),' '
  enddo
  close(4)
 
 end subroutine

end module