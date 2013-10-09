! small program to perform the movent of a electron on the network of a phot reciving molicule with environmental infulances

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  implicit none
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n
  real(EP)                                     :: tua,Temp,Kb,Hbar,pi,wc,Er,t,tmax,err,alpha
  complex(EP), dimension(:,:), allocatable     :: H,D1,D2,SO,D3
  complex(EP), dimension(:), allocatable       :: psi
  real(8)                                      :: cpu_start, cpu_end,dt
  character(len=80)                            :: fmt,filename,filename2,directmk
  CHARACTER(len=20)                            :: iter,step,node,pha,direct,typ
  logical                                      :: dirtest
  integer                                      :: info

  !open(6,file='test',STATUS='replace',ACTION='write')


  n=4
  err = 0.0001
  Kb = 1.3806488e-23
  pi = 3.14159265359
  Hbar =  6.62606957e-34


  allocate(H(n,n),psi(n))
  allocate(SO(n*n,n*n),D1(n,n),D2(n,n),D3(n,n))


!   H(1,1)=200;H(1,2)=-96;H(1,3)=5;H(1,4)=-4.4;H(1,5)=4.7;H(1,6)=-12.6;H(1,7)=-6.2
!   H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=6.8;H(2,5)=4.5;H(2,6)=7.4;H(2,7)=-0.3
!   H(3,1)=5;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0.8;H(3,6)=-8.4;H(3,7)=7.6
!   H(4,1)=-4.4;H(4,2)=6.8;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=-14.2;H(4,7)=-67
!   H(5,1)=4.7;H(5,2)=4.5;H(5,3)=0.8;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=-0.1
!   H(6,1)=12.6;H(6,2)=7.4;H(6,3)=-8.4;H(6,4)=-14.2;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
!   H(7,1)=-6.2;H(7,2)=-0.3;H(7,3)=7.6;H(7,4)=-67;H(7,5)=-0.1;H(7,6)=-38.3;H(7,7)=230


  Temp = 73; wc = 150; Er= 35
  tua = (2*pi*kb*Temp*Er)/(hbar*hbar*wc)

!   D1(1,1)=0;D1(1,2)=tua;D1(1,3)=tua;D1(1,4)=tua;D1(1,5)=tua;D1(1,6)=tua;D1(1,7)=tua
!   D1(2,1)=tua;D1(2,2)=0;D1(2,3)=tua;D1(2,4)=tua;D1(2,5)=tua;D1(2,6)=tua;D1(2,7)=tua
!   D1(3,1)=tua;D1(3,2)=tua;D1(3,3)=0;D1(3,4)=tua;D1(3,5)=tua;D1(3,6)=tua;D1(3,7)=tua
!   D1(4,1)=tua;D1(4,2)=tua;D1(4,3)=tua;D1(4,4)=0;D1(4,5)=tua;D1(4,6)=tua;D1(4,7)=tua
!   D1(5,1)=tua;D1(5,2)=tua;D1(5,3)=tua;D1(5,4)=tua;D1(5,5)=0;D1(5,6)=tua;D1(5,7)=tua
!   D1(6,1)=tua;D1(6,2)=tua;D1(6,3)=tua;D1(6,4)=tua;D1(6,5)=tua;D1(6,6)=0;D1(6,7)=tua
!   D1(7,1)=tua;D1(7,2)=tua;D1(7,3)=tua;D1(7,4)=tua;D1(7,5)=tua;D1(7,6)=tua;D1(7,7)=0
! 
!   D2(1,1)=0;D2(1,2)=0;D2(1,3)=0.5;D2(1,4)=0;D2(1,5)=0;D2(1,6)=0;D2(1,7)=0
!   D2(2,1)=0;D2(2,2)=0;D2(2,3)=0.5;D2(2,4)=0;D2(2,5)=0;D2(2,6)=0;D2(2,7)=0
!   D2(3,1)=0.5;D2(3,2)=0.5;D2(3,3)=1;D2(3,4)=0.5;D2(3,5)=0.5;D2(3,6)=0.5;D2(3,7)=0.5
!   D2(4,1)=0;D2(4,2)=0;D2(4,3)=0.5;D2(4,4)=0;D2(4,5)=0;D2(4,6)=0;D2(4,7)=0
!   D2(5,1)=0;D2(5,2)=0;D2(5,3)=0.5;D2(5,4)=0;D2(5,5)=0;D2(5,6)=0;D2(5,7)=0
!   D2(6,1)=0;D2(6,2)=0;D2(6,3)=0.5;D2(6,4)=0;D2(6,5)=0;D2(6,6)=0;D2(6,7)=0
!   D2(7,1)=0;D2(7,2)=0;D2(7,3)=0.5;D2(7,4)=0;D2(7,5)=0;D2(7,6)=0;D2(7,7)=0

!   H(1,1)=200;H(1,2)=-96;H(1,3)=0;H(1,4)=0;H(1,5)=0;H(1,6)=0;H(1,7)=0
!   H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=0;H(2,5)=0;H(2,6)=0;H(2,7)=0
!   H(3,1)=0;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0;H(3,6)=0;H(3,7)=0
!   H(4,1)=0;H(4,2)=0;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=0;H(4,7)=0
!   H(5,1)=0;H(5,2)=0;H(5,3)=0;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=0
!   H(6,1)=0;H(6,2)=0;H(6,3)=0;H(6,4)=0;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
!   H(7,1)=0;H(7,2)=0;H(7,3)=0;H(7,4)=0;H(7,5)=0;H(7,6)=-38.3;H(7,7)=230
  D3(1,1)=0;D3(1,2)=1;D3(1,3)=0.5;D3(1,4)=0;
  D3(2,1)=0;D3(2,2)=0;D3(2,3)=1;D3(2,4)=0;
  D3(3,1)=0;D3(3,2)=0;D3(3,3)=0;D3(3,4)=1;
  D3(4,1)=1;D3(4,2)=0;D3(4,3)=0;D3(4,4)=0;  

  !call B_SO_for_H(H,SO)
  !call A_D_to_SO(D1,SO)
  !call A_D_to_SO(D2,SO)

  !write(*,*)'H=',H

  psi(:)=cmplx(0,0);psi(1)=cmplx(1,0)
  t=1
  tmax=0.0
  dt=0.1
  SO(:,:)=cmplx(0,0)
  call L_make_DG(SO,D3,alpha)

!   open(3,file='results',STATUS='replace',ACTION='write')
!   call CBM(H,psi,t)
!   write(3,*)'t=',t,'Psi',psi*CONJG(psi)
!   psi(:)=cmplx(0,0);psi(1)=cmplx(1,0)
!   call rk45(H,psi,tmax,t,dt,err)
!   write(3,*)'t=',t,'Psi',psi*CONJG(psi)

  !
  !call Stand_to_SO(D3,SO)

  !write(*,*)'psi=',psi
!   open(4,file='SO',STATUS='replace',ACTION='write')
!   write(4,"(16F5.0)")real(SO)
!   close(4)




 
!enddo

 close(3)

!  deallocate(array, eighold, work1, work2, leftvectors, rightvectors)
!  deallocate(psi, H, theta0, theta1, theta2, tot)
end program



 