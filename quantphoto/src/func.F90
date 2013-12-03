module solver
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
  integer(kdp)				    :: N
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
 
 end subroutine
 
 
  ! computes the second moment for a time series
 subroutine second_moment(in, moment)
  complex(kdp), dimension(:), intent(in)    :: in
  complex(kdp),intent(out)                  :: moment
  complex(kdp), dimension(:), allocatable   :: out 
  integer(kdp)				    :: N
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
 
 end subroutine
 
 subroutine k_product()
 
 end subroutine
 
 subroutine par_trace()
 
 end subroutine
 
end module