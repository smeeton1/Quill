module func
use eigenSolve
use until
!use special
implicit none


 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp

 contains
 
 
 ! computes the first moment for a time series
 subroutine first_moment(in, moment)
  complex(kdp), dimension(:), intent(in)    :: in
  complex(kdp),intent(out)                  :: moment
  complex(kdp), dimension(:), allocatable   :: out 
  integer(kdp)				    :: N, i
  complex(kdp)                              :: top, bottom
  
  N=size(in)
  allocate(out(N))
  call fftw_1D(in, out)
  
  do i=1,N
    top=top+out(i)*CMPLX(i,0)*CMPLX(i,0)
  enddo
  
  do i=1,N
    bottom=bottom+out(i)
  enddo
   
  moment=sqrt(top/bottom)
  
  deallocate(out)
 
 end subroutine
 
 
  ! computes the second moment for a time series
 subroutine second_moment(in, moment)
  complex(kdp), dimension(:), intent(in)    :: in
  complex(kdp),intent(out)                  :: moment
  complex(kdp), dimension(:), allocatable   :: out 
  integer(kdp)				    :: N, i
  complex(kdp)                              :: top, bottom
  
  N=size(in)
  allocate(out(N))
  call fftw_1D(in, out)
  
  do i=1,N
    top=top+out(i)*CMPLX(i,0)
  enddo
  
  do i=1,N
    bottom=bottom+out(i)
  enddo
   
  moment=top/bottom
  
  deallocate(out)
 
 end subroutine
 
 subroutine write_moments(filename,in)
  complex(kdp), dimension(:,:), intent(in)  :: in
  complex(kdp)                              :: moment1, moment2
  integer(kdp)				    :: n,n2, i
  character(80), intent(in)                 :: filename
  character(90)                             :: filename2
 
  n=size(in,1)
  n2=size(in,2)
  write(filename2,'(a,a)')trim(filename),'_moments'
  open(5,file=filename2,STATUS='replace',ACCESS='append',ACTION='write')
  write(5,'(a,a,a,a,a)')'Node',' ','first moment',' ','second moment'
  do i=1,n2
   call first_moment(in(1:n,i), moment1)
   call second_moment(in(1:n,i), moment2)
   write(5,'(I3,2F8.5)')i,real(moment1),real(moment2)
  enddo
  close(5)
 end subroutine
 
 subroutine Rand_Weight(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    do j=1,m
      if(D(i,j).ne.0)then
       D(i,j)=5*rand()
      endif
    enddo
  enddo
  
  
 end subroutine
 
 subroutine Mat_Gen(D,porb)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  real,intent(in)                            :: porb
  integer				     :: i,j,n,m
  n=size(D,1)
  m=size(D,2)
  
  D(:,:)=cmplx(0,0);
  do i=1,n
    do j=1,m
      if((rand().gt.porb).and.(i.ne.j))then
       D(i,j)=5*rand()
      endif
    enddo
  enddo
  
  
 end subroutine
 
 subroutine row_norm(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  complex(kdp)  			     :: tot
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    tot=cmplx(0.0,0.0)
    do j=1,m
       tot = tot + D(i,j)
    enddo
    if(tot.ne.cmplx(0.0,0.0))then
     do j=1,m
        D(i,J) = D(i,j)/tot
     enddo
    endif
  enddo
    
 end subroutine
 
  subroutine col_norm(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  complex(kdp)  			     :: tot
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    tot=cmplx(0.0,0.0)
    do j=1,m
       tot = tot + D(j,i)
    enddo
    if(tot.ne.cmplx(0.0,0.0))then
     do j=1,m
        D(j,i) = D(j,i)/tot
     enddo
    endif
  enddo
    
 end subroutine
 
 subroutine k_product()
 
 end subroutine
 
 subroutine par_trace()
 
 end subroutine
 
end module