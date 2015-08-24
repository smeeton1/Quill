program discreettest
use eigenSolve
use until
use func
use Densutil
use supperoperator
use diswalk
implicit none

  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n,t,l,start
  real(EP)                                     :: alpha,  beta
  complex(EP), dimension(:), allocatable       :: rhod
  character(len=90)                            :: filename
  logical                                      :: r


 t=50
 n=16
 allocate(rhod(n*n))
 do l=1,11
     alpha=real(l-1)*0.1
     do k=1,11
       beta=real(k-1)*0.1
       rhod(:)=cmplx(0,0)
       rhod(1)=1
       do i=1,n
        !rhod(i+(i-1)*n)=cmplx(1./n,0.0)
        !rhod(i+(i-1)*n+1)=cmplx(1./n,0.0)
       enddo
       
       !call write_Vec_real('rhod',rhod)
       
       write(filename,'(a,2F4.2)')'results/distest/test',alpha,beta
       write(*,*)filename
  
       call Dis_Dir_8line(rhod, t, alpha, beta, filename, r)
       
    enddo
 enddo
 
end program
