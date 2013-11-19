module until
 implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains

! routine to read in a matrix from a file. Assumes that the file has only number. 
! mat should be sized for expected size of matrix.
 subroutine Read_Mat(Filename,mat)
  complex(kdp),dimension(:,:),intent(inout):: mat
  character(len=80),intent(in)             :: filename 
  integer                                  :: n, i, j
  
  n=size(mat)
  open(4,file=filename,ACCESS='append',ACTION='read')
  do i=1,n
   do j=1,n
      read(4,"(F15.10)") mat(i,j)
   enddo
  enddo
  
  
  close(4)
  
 
 end subroutine

 subroutine write_Mat(filename,mat)
  complex(kdp),dimension(:,:),intent(in):: mat
  character(len=80),intent(in)          :: filename 
  integer                               :: n, m, i, j

  n=size(mat,1)
  m=size(mat,2)

  open(4,file=filename,ACCESS='append',ACTION='write')
  write(4,*)''
  do i=1,n
   do j=1,m
       write(4,"(F15.10,a,a,F15.10,a)",ADVANCE='no')real(mat(i,j)),' ','i',aimag(mat(i,j)),' '
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
    write(4,"(F15.10,a,a,F15.10,a)",ADVANCE='no')real(vec(i)),' ','i',aimag(vec(i)),' '
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