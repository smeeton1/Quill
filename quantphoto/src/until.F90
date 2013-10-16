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

  n=size(mat,1)

  open(4,file=filename,ACCESS='append',ACTION='write')
  write(4,*)''
  do i=1,n
   do j=1,n
       write(4,"(F6.3,a,a,F6.3,a)",ADVANCE='no')real(mat(i,j)),' ','i',aimag(mat(i,j)),' '
   enddo 
   write(4,*)''
  enddo
  close(4)
 
 end subroutine

 subroutine write_Vec(filename,vec)
  complex(kdp),dimension(:),intent(in)  :: vec
  character(len=80),intent(in)          :: filename 
  integer                               :: n, i

  n=size(vec)

  open(4,file=filename,ACCESS='append',ACTION='write')
  write(4,*)''
  do i=1,n
    write(4,"(F6.3,a,a,F6.3,a)",ADVANCE='no')real(vec(i)),' ','i',aimag(vec(i)),' '
  enddo
  close(4)
 
 end subroutine

  
 subroutine make_dir(dirname)
  character(len=80),intent(in)          :: dirname 
  character(len=90)			:: dirnames,directmk
  logical                               :: dirtest

  write(dirnames,'(a)')trim(dirname)
  inquire(file=dirnames,exist=dirtest)
 if(.not.dirtest)then
  write(directmk,'(a,a)')'mkdir ',trim(dirname)
  Call system(directmk)
 endif


 end subroutine


end module