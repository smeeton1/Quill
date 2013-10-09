

program SOtest
  use supperoperator
  implicit none
complex(8), dimension(:,:), allocatable     ::H,D1,SO
integer ::n

n=3
allocate(SO(n*n,n*n),D1(n,n),H(n,n))

H(1,1)=1;H(1,2)=2;H(1,3)=3
H(2,1)=4;H(2,2)=5;H(2,3)=6
H(3,1)=7;H(3,2)=8;H(3,3)=9

D1(1,1)=1;D1(1,2)=1;D1(1,3)=1
D1(2,1)=1;D1(2,2)=1;D1(2,3)=1
D1(3,1)=1;D1(3,2)=1;D1(3,3)=1

 SO(:,:)=cmplx(0,0)
call B_SO_for_H(H,SO)
  open(4,file='SO_H',STATUS='replace',ACTION='write')
  write(4,*)real(SO)
  close(4)
! call A_D_to_SO(D1,SO)
!   open(4,file='SO_D',STATUS='replace',ACTION='write')
!   write(4,*)real(SO)
!   close(4)

end program