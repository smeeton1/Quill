! small program to perform a directed quantum walk on 3 node graphs.

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  use Densutil
  implicit none

 
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n,t
  real(EP)                                     :: alpha
  complex(EP), dimension(:,:), allocatable     :: D
  complex(EP), dimension(:), allocatable       :: rho
  character(len=80)                            :: fmt,filename,dirname
  logical                                      :: v
  
  
  n=3
  t=10
  alpha=0.8
  allocate(D(n,n),rho(n*n))
  v=.false.
  
  ! graph 1
  filename='graph_1'
  D(1,1)=0.0;D(1,2)=1./2.;D(1,3)=1./2.
  D(2,1)=1./2.;D(2,2)=0.0;D(2,3)=1./2.
  D(3,1)=1./2.;D(3,2)=1./2.;D(3,3)=0.0
  rho(:)=cmplx(0,0);rho(1)=cmplx(1,0)
  write(*,*)filename
  call Dir_Gra_Run(D, rho, t, alpha, filename, v)
  
  ! graph 2
  filename='graph_2'
  D(1,1)=0.0;D(1,2)=0.0;D(1,3)=1.
  D(2,1)=1./2.;D(2,2)=0.0;D(2,3)=1./2.
  D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
  rho(:)=cmplx(0,0);rho(1)=cmplx(1,0)
  write(*,*)filename
  call Dir_Gra_Run(D, rho, t, alpha, filename, v)
  
  ! graph 3
  filename='graph_3'
  D(1,1)=0.0;D(1,2)=1./2.;D(1,3)=1./2.
  D(2,1)=0.0;D(2,2)=0.0;D(2,3)=1.
  D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
  rho(:)=cmplx(0,0);rho(1)=cmplx(1,0)
  write(*,*)filename
  call Dir_Gra_Run(D, rho, t, alpha, filename, v)
  
  ! graph 4
  filename='graph_4'
  D(1,1)=0.0;D(1,2)=1.;D(1,3)=0.0
  D(2,1)=0.0;D(2,2)=0.0;D(2,3)=1.
  D(3,1)=1.;D(3,2)=0.0;D(3,3)=0.0
  rho(:)=cmplx(0,0);rho(1)=cmplx(1,0)
  write(*,*)filename
  call Dir_Gra_Run(D, rho, t, alpha, filename, v)
  
  ! graph 5
  filename='graph_5'
  D(1,1)=0.0;D(1,2)=1./2.;D(1,3)=1./2.
  D(2,1)=0.0;D(2,2)=0.0;D(2,3)=0.0
  D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
  rho(:)=cmplx(0,0);rho(1)=cmplx(1,0)
  write(*,*)filename
  call Dir_Gra_Run(D, rho, t, alpha, filename, v)
  
  

end program
  