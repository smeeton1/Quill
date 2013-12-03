module eigenSolve
  use, intrinsic :: iso_c_binding
  
  implicit none 
#include <fftw3.f>
!     interface 
!        subroutine dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
!           complex(8), dimension(:)                :: in
!           complex(8), dimension(:)                :: out 
!           integer(8)				  :: plan,N
!           integer                                 :: FFTW_FORWARD,FFTW_ESTIMATE
!        end subroutine dfftw_plan_dft_1d
!        
!        subroutine dfftw_execute_dft(plan, in, out)
!           complex(8), dimension(:)                :: in
!           complex(8), dimension(:)                :: out 
!           integer(8)				  :: plan
!        end subroutine dfftw_execute_dft
!        
!        subroutine dfftw_destroy_plan(plan)
!           integer(8)				  :: plan
!        end subroutine dfftw_destroy_plan       
!     end interface
  
  
  
  
 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
 
 
 ! function to calculate the eiginvaules of a matrix using lapack
  subroutine eigen(mat,v)
    complex(8), dimension(:,:), intent(in)  :: mat
    complex(8), dimension(:), intent(out)   :: v 
    interface 
       subroutine zgeev(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, & 
            & LWORK, RWORK, INFO) 
         character :: JOBVL, JOBVR 
         integer :: INFO, LDA, LDVL, LDVR, LWORK, N 
         real(8) :: RWORK(*) 
         complex(8) :: A(*), VL(*), VR(*), W(*), WORK(*) 
       end subroutine zgeev
    end interface
    complex(8), dimension(:, :), allocatable :: array, leftvectors, rightvectors 
    complex(8), dimension(:), allocatable :: eigenvalues, work1 
    real(8), dimension(:), allocatable :: work2 
    integer :: n, info
    !
    n = size(mat(1,:))
    allocate(array(n, n), eigenvalues(n), work1(2*n), work2(2*n), & 
         & leftvectors(n, n), rightvectors(n, n)) 
    array(:,:) = mat(:,:)
    call zgeev('V', 'V', n, array, n, eigenvalues, leftvectors, n, rightvectors,& 
         & n, work1, 2*n, work2, info) 
    v(:) = eigenvalues(:)
    deallocate(array, eigenvalues, work1, work2, leftvectors, rightvectors)
  end subroutine eigen
  
  
  ! dose the forward fft for complex vectors in 1 dimention using fftw it is not normallized.
  subroutine fftw_1D(in, out)
    complex(8), dimension(:), intent(in)    :: in
    complex(8), dimension(:), intent(out)   :: out 
    integer(8)				    :: plan, N
    
    N=size(in)
    call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, in, out)
    call dfftw_destroy_plan(plan)
  
  end subroutine
  
  
  
  
  
  
  
end module eigenSolve