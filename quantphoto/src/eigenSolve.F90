module eigenSolve
  implicit none 
 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
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
end module eigenSolve