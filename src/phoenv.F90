! small program to perform the movent of a electron on the network of a phot reciving molicule with environmental infulances

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  use Densutil 
  use therm
  implicit none
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j, l, k, n, m, sink, Tend,p
  integer, dimension(:), allocatable           :: conect
  real(EP)                                     :: Temp,Kb,Hbar,pi,wc,Er,t,tmax,err
  complex(EP), dimension(:,:), allocatable     :: H,SO
  complex(EP), dimension(:), allocatable       :: rho
  real(8)                                      :: cpu_start, cpu_end,dt
  complex(EP)				       :: norm, coli, strength, alpha, E, Ve, tua
  character(len=80)                            :: fmt,filename,filename2,dirname
  CHARACTER(len=20)                            :: iter,step,node,pha,direct,typ
  logical                                      :: dirtest,work,v,r
  integer                                      :: info

  !open(6,file='test',STATUS='replace',ACTION='write')


  n=7
  err = 0.00001
  Kb = 1.3806488e-23
  pi = 3.14159265359
  Hbar = 1.0! 6.62606957e-34
  alpha = 0.8
  
  Tend = 50
  v = .false.
  r = .true.
  

  dt = 0.1

  allocate(H(n,n),rho(n*n),conect(7))

  conect(1)=1
  conect(2)=2
  conect(3)=3
  conect(1)=4
  conect(2)=5
  conect(3)=6
  conect(1)=7
  
  sink=3
  strength=cmplx(0.5,0.0)
  
  rho(5+4*n)=cmplx(1.0,0.0)
  
  H(1,1)=200;H(1,2)=-96;H(1,3)=5;H(1,4)=-4.4;H(1,5)=4.7;H(1,6)=-12.6;H(1,7)=-6.2
  H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=6.8;H(2,5)=4.5;H(2,6)=7.4;H(2,7)=-0.3
  H(3,1)=5;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0.8;H(3,6)=-8.4;H(3,7)=7.6
  H(4,1)=-4.4;H(4,2)=6.8;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=-14.2;H(4,7)=-67
  H(5,1)=4.7;H(5,2)=4.5;H(5,3)=0.8;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=-0.1
  H(6,1)=12.6;H(6,2)=7.4;H(6,3)=-8.4;H(6,4)=-14.2;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
  H(7,1)=-6.2;H(7,2)=-0.3;H(7,3)=7.6;H(7,4)=-67;H(7,5)=-0.1;H(7,6)=-38.3;H(7,7)=230
 
  wc = 150; Er= 35
  do i=-2,2 
   Temp = 273.0-real(i)*20; 
   tua = cmplx((2*pi*kb*Temp*Er)/(hbar*hbar*wc),0.0)
   write(dirname,'(a,I3.3)')'results/photoFull',Int(Temp)
   call make_dir(dirname)
   do l=1,11
     coli=cmplx(real(l-1)*0.1,0.0)
     write(filename,'(a,a,a,F4.2)')trim(dirname),'/','start5w',real(coli)
     call ele_tran(H,rho,Tend,dt,v,r,filename,strength,tua,coli,sink,conect)
   enddo
  enddo

! 
!   H(1,1)=200;H(1,2)=-96;H(1,3)=0;H(1,4)=0;H(1,5)=0;H(1,6)=0;H(1,7)=0
!   H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=0;H(2,5)=0;H(2,6)=0;H(2,7)=0
!   H(3,1)=0;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0;H(3,6)=0;H(3,7)=0
!   H(4,1)=0;H(4,2)=0;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=0;H(4,7)=0
!   H(5,1)=0;H(5,2)=0;H(5,3)=0;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=0
!   H(6,1)=0;H(6,2)=0;H(6,3)=0;H(6,4)=0;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
!   H(7,1)=0;H(7,2)=0;H(7,3)=0;H(7,4)=0;H(7,5)=0;H(7,6)=-38.3;H(7,7)=230
! ! 
!   write(dirname,'(a)')'photodag'  
!   call make_dir(dirname)
!   do l=1,11
!     coli=cmplx(real(l-1)*0.1,0.0)
!     write(filename,'(a,a,a,F4.2)')trim(dirname),'/','start5w',real(coli)
!     call ele_tran(H,rho,Tend,v,r,filename,strength,tua,coli,sink,conect)
!   enddo

  deallocate(H,rho,conect)
  
  m=5
  p=5
  n=m*p
  sink=3+m*(p-1)

  
  allocate(H(n,n),rho(n*n),conect(3))
  rho(:)=cmplx(0.0,0.0)
    do i=1,m
      rho(i+(i-1)*n)=cmplx(1./m,0.0)
  enddo
  conect(1)=2+m*(p-1)-1
  conect(2)=2+m*(p-1)+1
  conect(3)=2+m*(p-2)
!   
!   
! 
!  
  do i=-2,2 
   Temp = 273.0-real(i)*20; 
   tua = cmplx((2*pi*kb*Temp*Er)/(hbar*hbar*wc),0.0)
   do k=1,6
    do j=1,6
   
     E=cmplx(100*(k-1)/Temp,0.0) 
     Ve=cmplx(100*(j-1),0.0)
     alpha= cmplx(1.0/Temp,0.0)
     write(dirname,'(a,I2.2,a,I2.2,a,I3.3)')'results/lattex-',k,'-',j,'-',int(Temp)
     call make_dir(dirname)
     call make_lattic(H,E,Ve,alpha,m)
     do l=1,11
       coli=cmplx(real(l-1)*0.1,0.0)
       write(filename,'(a,a,a,F4.2)')trim(dirname),'/','start1r-',real(coli)
       call ele_tran(H,rho,Tend,dt,v,r,filename,strength,tua,coli,sink,conect)
     enddo
    enddo
   enddo
  enddo
  


!  deallocate(array, eighold, work1, work2, leftvectors, rightvectors)
!  deallocate(rho, H, theta0, theta1, theta2, tot)
end program



 