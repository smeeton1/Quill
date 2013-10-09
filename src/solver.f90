module solver
use eigenSolve
use special
!implicit none
 contains

 subroutine rk45(SO,Rho,t1,t2,h,err)
   integer, parameter :: EP = selected_real_kind(15)
   complex(EP), dimension(:,:), intent(in)     :: SO
   complex(EP), dimension(:), intent(inout)    :: Rho
   real(EP),intent(inout)                      :: h
   real(EP),intent(in)                         :: t1,t2,err

   integer                                     :: n,m,j,i
   real(EP)                                    :: diff, t 
   complex(EP), dimension(:),allocatable       :: f1,f2,f3,f4,f5,f6
   complex(EP), dimension(:),allocatable       :: rho1,rho2

   
   n=size(Rho)
   diff=1.0

   allocate(f1(n),f2(n),f3(n),f4(n),f5(n),f6(n),rho1(n),rho2(n))
   m=ceiling((t2-t1)/h)
   t=t1
   write(*,*)t,h,m
   do i=0,m-1
    do j=1,50
      f1 = h*matmul(SO,Rho)
      f2 = h*matmul(SO,(Rho + f1*h/4))!*(t+h/4)
      f3 = h*matmul(SO,(Rho + (3*f1/32 + 9*f2/32)*h))!*(t+3*h/8)
      f4 = h*matmul(SO,(Rho + (1932*f1/2192 - 7200*f2/2197 + 7296*f3/2197)*h))!*(t+12*h/13)
      f5 = h*matmul(SO,(Rho + (438*f1/216 - 8*f2 + 3680*f3/513 - 845*f4/4104)*h))!*(t+h)
      f6 = h*matmul(SO,(Rho - (8*f1/27 + 2*f2 - 3544*f3/2565 + 1859*f4/4104 - 11*f5/40)*h))!*(t+h/2) 
     
      rho1 = Rho + ((25/216)*f1 + (1408/2565)*f3 + (2197/4104)*f4 -(1/5)*f5)*h
     
      rho2 = Rho + ((16/135)*f1 + (6656/12825)*f3 + (28561/56430)*f4 - (9/50)*f5 + (2/55)*f6)*h
      !write(*,*)f1

      diff = maxval(real(CONJG(rho1)*rho1-CONJG(rho2)*rho2))

      if(diff.lt.err) then
	goto 99
      else
	h=(0.840896*(h*err/abs(diff))**(1/4))*h
	m=(ceiling(0.840896*(m*err/abs(diff))**(1/4)))*m
      endif

    end do
   

   99 Rho = rho1
   t = t + h
   write(*,*)sum(real(conjg(Rho)*Rho)), t,h,t+h,m
   write(*,*)f1
   enddo

 end subroutine



 ! dose the Chebyshev aproximation requires the time the hamiltonian matrix and the wave function.
! works for when H has only real eigenvalues
 subroutine CBM(H,psi,t)
   integer, parameter :: EP = selected_real_kind(15)
   complex(EP), dimension(:,:), intent(in)     :: H
   complex(EP), dimension(:), intent(inout)    :: psi
   real(EP),intent(in)                         :: t
   
   integer                                    :: i, j, k, n, M, ND, em, en
   real(EP), dimension(:,:), allocatable      :: theta0, theta1, theta2, muthold
   complex(EP), dimension(:,:), allocatable   :: tot
   complex(EP), dimension(:), allocatable     :: eighold!, CB1, CB2
   real(EP)                                   :: emax, emin, alpha
   integer                                    :: info
   
   n=size(H,1)
   allocate(theta0(n,n),theta1(n,n),theta2(n,n),tot(n,n),muthold(n,n),eighold(n))
   !allocate(CB1(3))
   call eigen(H,eighold)

   ND=4
   em = maxloc(real(abs(eighold)),1)
   en = minloc(real(abs(eighold)),1)
   emax = eighold(em)
   emin = eighold(en)

   !write(*,*)eighold
   !   write(*,*)'emax=',emax
   !   write(*,*)'emin=',emin



   !write(*,*)'run for T=',t
     alpha=(emax-emin)*t/2
     M=CEILING(1.2*real(alpha)+30)
     !allocate(CB2(M+1))
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

     !call WBsJA(alpha,0,0,ND,CB1)
     !call wbsja(alpha,1,1,ND,CB2)
     do i=1,n
       do j=1,n
         tot(i,j)=BesJN(0,alpha)*theta0(i,j)+2*cmplx(0,1)*BesJN(1,alpha)*theta1(i,j)
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
           !call wbsja(alpha,0,i,ND,CB2)
           tot(j,k)=tot(j,k)+2*(cmplx(0,1)**i)*BesJN(i,alpha)*theta2(j,k)
         enddo
       enddo
       theta0=theta1
       theta1=theta2
     end do
   psi=exp(cmplx(0,-1)*(emax+emin)*t/2)*matmul(tot,psi)

end subroutine

subroutine CBM_comp(H,psi,t)
   integer, parameter :: EP = selected_real_kind(15)
   complex(EP), dimension(:,:), intent(in)     :: H
   complex(EP), dimension(:), intent(inout)    :: psi
   real(EP),intent(in)                         :: t
   
   integer                                    :: i, j, k, n, M, ND, em, en
   real(EP), dimension(:,:), allocatable      :: theta0, theta1, theta2, muthold
   complex(EP), dimension(:,:), allocatable   :: tot
   complex(EP), dimension(:), allocatable     :: eighold!, CB1, CB2
   complex(EP)                                :: emax, emin, alpha, hold
   integer                                    :: info
   
   n=size(H,1)
   allocate(theta0(n,n),theta1(n,n),theta2(n,n),tot(n,n),muthold(n,n),eighold(n))
   !allocate(CB1(3))
   call eigen(H,eighold)

   ND=4
   em = maxloc(real(abs(eighold)),1)
   en = minloc(real(abs(eighold)),1)
   emax = eighold(em)
   emin = eighold(en)

   !write(*,*)eighold
   !   write(*,*)'emax=',emax
   !   write(*,*)'emin=',emin



   !write(*,*)'run for T=',t
     alpha=(emax-emin)*t/2
     M=CEILING(1.2*real(alpha)+30)
     !allocate(CB2(M+1))
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

     
     hold = besseljv_complex(1, alpha)
     do i=1,n
       do j=1,n
         tot(i,j)=besseljv_complex(0,alpha)*theta0(i,j)+2*cmplx(0,1)*hold*theta1(i,j)
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
           !call wbsja(alpha,0,i,ND,CB2)
           tot(j,k)=tot(j,k)+2*(cmplx(0,1)**i)*besseljv_complex(i, alpha)*theta2(j,k)
         enddo
       enddo
       theta0=theta1
       theta1=theta2
     end do
   psi=exp(cmplx(0,-1)*(emax+emin)*t/2)*matmul(tot,psi)

end subroutine


end module