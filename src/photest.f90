! small program to perform the movent of a electron on the network of a phot reciving molicule


module eigenSolve
  implicit none 

contains
  subroutine eigen(mat,v)
    complex(8), dimension(:,:), intent(in)  :: mat
    real(8), dimension(:), intent(out) :: v 
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
    v(:) = real(eigenvalues(:))
    deallocate(array, eigenvalues, work1, work2, leftvectors, rightvectors)
  end subroutine eigen
end module eigenSolve



program photest
  use eigenSolve
  implicit none
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n, M, tmax,t
  real(8)                                      :: val
  real(EP), dimension(:,:), allocatable        :: theta0,theta1,theta2,muthold
  complex(EP), dimension(:,:), allocatable     ::tot,H
  complex(EP), dimension(:), allocatable       :: psi,psir
  real(8), dimension(:), allocatable           :: eighold
  real(EP)                                     :: emax, emin, alpha
  real(8)                                      :: cpu_start, cpu_end
  character(len=80)                            :: fmt,filename,filename2,directmk
  CHARACTER(len=20)                            :: iter,step,node,pha,direct,typ
  logical                                      :: dirtest
  complex(8), dimension(:, :), allocatable     :: array, leftvectors, rightvectors 
  complex(8), dimension(:), allocatable        ::  work1 
  real(8), dimension(:), allocatable           :: work2 
  integer                                      ::  info

  !open(6,file='test',STATUS='replace',ACTION='write')


  n=7
  tmax=500.0



  allocate(H(n,n),theta0(n,n),theta1(n,n),theta2(n,n),tot(n,n),muthold(n,n),psi(n),psir(n),eighold(n))
    allocate(array(n, n), work1(2*n), work2(2*n), leftvectors(n, n), rightvectors(n, n)) 


!   H(1,1)=200;H(1,2)=-96;H(1,3)=5;H(1,4)=-4.4;H(1,5)=4.7;H(1,6)=-12.6;H(1,7)=-6.2
!   H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=6.8;H(2,5)=4.5;H(2,6)=7.4;H(2,7)=-0.3
!   H(3,1)=5;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0.8;H(3,6)=-8.4;H(3,7)=7.6
!   H(4,1)=-4.4;H(4,2)=6.8;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=-14.2;H(4,7)=-67
!   H(5,1)=4.7;H(5,2)=4.5;H(5,3)=0.8;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=-0.1
!   H(6,1)=12.6;H(6,2)=7.4;H(6,3)=-8.4;H(6,4)=-14.2;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
!   H(7,1)=-6.2;H(7,2)=-0.3;H(7,3)=7.6;H(7,4)=-67;H(7,5)=-0.1;H(7,6)=-38.3;H(7,7)=230

  H(1,1)=200;H(1,2)=-96;H(1,3)=0;H(1,4)=0;H(1,5)=0;H(1,6)=0;H(1,7)=0
  H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=0;H(2,5)=0;H(2,6)=0;H(2,7)=0
  H(3,1)=0;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0;H(3,6)=0;H(3,7)=0
  H(4,1)=0;H(4,2)=0;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=0;H(4,7)=0
  H(5,1)=0;H(5,2)=0;H(5,3)=0;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=0
  H(6,1)=0;H(6,2)=0;H(6,3)=0;H(6,4)=0;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
  H(7,1)=0;H(7,2)=0;H(7,3)=0;H(7,4)=0;H(7,5)=0;H(7,6)=-38.3;H(7,7)=230

  
  !write(*,*)'H=',H

  psi(:)=cmplx(0,0);psi(1)=cmplx(1,0)

  !write(*,*)'psi=',psi

  call eigen(H,eighold)
  !array(:,:)=H(:,:)
  !call zgeev('V', 'V', n, array, n, eighold, leftvectors, n, rightvectors,& 
         !& n, work1, 2*n, work2, info) 

  emax = maxval(eighold)
  emin = minval(eighold)

  !write(*,*)eighold
  write(*,*)'emax=',emax
  write(*,*)'emin=',emin
  open(3,file='results_Hm',STATUS='replace',ACTION='write')

do t=50,tmax,50
  write(*,*)'run for T=',t
  alpha=(emax-emin)*t/2
  M=CEILING(1.2*alpha+30)
  theta0(:,:)=0
  theta1(:,:)=0
  theta2(:,:)=0
 do i=1,n
  theta0(i,i)=1
 enddo

  do i=1,n
   do j=1,n
      theta1(i,j)=-(2*H(i,j)-(emax+emin)*theta0(i,j))/(emax-emin)
   enddo
  enddo


 do i=1,n
   do j=1,n
      tot(i,j)=BESJN(0,alpha)*theta0(i,j)+2*cmplx(0,1)*BESJN(1,alpha)*theta1(i,j)
   enddo
 enddo


 do i=2,M
  muthold=matmul(H,theta1)
  do j=1,n
    do k=1,n
      theta2(j,k)=-2*(2*muthold(j,k)-(emax+emin)*theta1(j,k))&
         &/(emax-emin)-theta0(j,k)
    enddo
  enddo
  do j=1,n
    do k=1,n
      tot(j,k)=tot(j,k)+2*(cmplx(0,1)**i)*BESJN(i,alpha)*theta2(j,k)
    enddo
  enddo
  theta0=theta1
  theta1=theta2
 end do
 psir=exp(cmplx(0,-1)*(emax+emin)*t/2)*matmul(tot,psi)

 write(3,*)'t=',t,'Psi',psir*CONJG(psir)
 
enddo

 close(3)

!  deallocate(array, eighold, work1, work2, leftvectors, rightvectors)
!  deallocate(psi, H, theta0, theta1, theta2, tot)
end program



 