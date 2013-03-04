! small program to perform the movent of a electron on the network of a phot reciving molicule

program photest

  integer                                      :: i, j,k, n, M
  integer                                      :: steps, v1, v2, c1, c2
  integer, dimension(:), allocatable           :: p
  complex(8), dimension(:,:,:,:), allocatable  :: space,spaceCoin
  real(8)                                      :: val
  complex(8), dimension(:,:), allocatable      :: H,theta0,theta1,theta2,tot
  complex(8), dimension(:), allocatable        :: psi
  complex(8)                                   :: emax, emin, alpha
  real(8)                                      :: t
  real(8)                                      :: cpu_start, cpu_end
  character(len=80)                            :: fmt,filename,filename2,directmk
  CHARACTER(len=20)                            :: iter,step,node,pha,direct,typ
  logical                                      :: dirtest

  n=7
  allocate(H(n,n),theta0(n,n),theat1(n,n),theta2(n,n),tot(n,n),psi(n))
  
  H(1,1)=200;H(1,2)=-96;H(1,3)=5;H(1,4)=-4.4;H(1,5)=4.7;H(1,6)=-12.6;H(1,7)=-6.2
  H(2,1)=-96;H(2,2)=320;H(2,3)=33.1;H(2,4)=6.8;H(2,5)=4.5;H(2,6)=7.4;H(2,7)=-0.3
  H(3,1)=5;H(3,2)=33.1;H(3,3)=0;H(3,4)=-51.1;H(3,5)=0.8;H(3,6)=-8.4;H(3,7)=7.6
  H(4,1)=-4.4;H(4,2)=6.8;H(4,3)=-51.1;H(4,4)=110;H(4,5)=-76.6;H(4,6)=-14.2;H(4,7)=-67
  H(5,1)=4.7;H(5,2)=4.5;H(5,3)=0.8;H(5,4)=-76.6;H(5,5)=270;H(5,6)=78.3;H(5,7)=-0.1
  H(6,1)=12.6;H(6,2)=7.4;H(6,3)=-8.4;H(6,4)=-14.2;H(6,5)=78.3;H(6,6)=420;H(6,7)=38.3
  H(7,1)=-6.2;H(7,2)=-0.3;H(7,3)=7.6;H(7,4)=-67;H(7,5)=-0.1;H(7,6)=-38.3;H(7,7)=230
  
  psi=cmplx(0,0);psi(1,1)=cmplx(1,0)


  alpha=(emax-emini)*t/2
  M=int(1.2*alpha+30)
  theta0=cmplx(0,0)
 do i=1,n
  theta0(i,i)=1
   do j=1,n
      theta1(i,j)=2*H(i,j)-(emax+emin)/(emax-emin)*theta0(i,j)
   enddo
 enddo

 do i=1,n
   do j=1,n
      tot(i,j)=BESJN(0,alpha)theta0(i,j)+2*cmplx(0,1)*BESJN(1,alpha)*theta1(i,j)
   enddo
 enddo

 do i=2,M
  do j=1,n
    do k=1,n
      theta2(j,k)=-2*(2*H*theta1(j,k)-(emax+emin)*theta1(j,k))/(emax-emin)-theta0(j,k)
    enddo
  enddo
  do j=1,n
    do k=1,n
      tot(j,k)=tot(j,k)+2*cmplx(0,1)^i*BESJN(i,alpha)*theta2(j,k)
    enddo
  enddo
  theta0=theta1
  theta2=theta1
 end do
 
 psi=exp(cmplx(0,-1)*(emax+emin)*t/2)*matmul(tot*psi)


end program