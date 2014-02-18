! small program to perform a directed quantum walk on 3 node graphs.

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  use Densutil
  use func
  implicit none

 
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n,t,l
  real(EP)                                     :: alpha
  complex(EP), dimension(:,:), allocatable     :: D
  complex(EP), dimension(:), allocatable       :: rho
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,AALPHA,filenamebase,dirnamebase
  logical                                      :: v, r
  
  
  

 CALL GET_COMMAND_ARGUMENT(1,filenamebase)
 CALL GET_COMMAND_ARGUMENT(2,dirnamebase)
 CALL GET_COMMAND_ARGUMENT(3,nn)
 CALL GET_COMMAND_ARGUMENT(4,time)
 read(nn,*)n
 read(time,*)t
 
 allocate(D(n,n),rho(n*n))
 k=1
!   n=3
!   t=100
!   alpha=0.1


  v=.false.
  r=.true.
! write(*,*)size(D,1)
  do j=1,k
    write(dirname,'(a,I5.5)')trim(dirnamebase),j
    call make_dir(dirname)
    call Mat_Gen(D,0.2)
    write(filename,'(3a)')trim(dirname),'/','D'
    call write_Mat(filename,D) 
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
!  call write_Mat('D',D) 
    do l=1,10
      alpha=0.1*real(l)
      write(filename,'(3a,F3.1)')trim(dirname),'/',trim(filenamebase),alpha
      call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
    enddo
  enddo
  


  
!   ! graph 1
!   filename='graph_1_1'
!   D(1,1)=0.0;D(1,2)=1.;D(1,3)=4.
!   D(2,1)=1./2.;D(2,2)=0.0;D(2,3)=1./2.
!   D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
!   rho(:)=cmplx(0,0);rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!   write(*,*)filename
!   call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
!   
!   ! graph 2
!   filename='graph_2'
!   D(1,1)=0.0;D(1,2)=0.0;D(1,3)=1.
!   D(2,1)=1./2.;D(2,2)=0.0;D(2,3)=1./2.
!   D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
!   rho(:)=cmplx(0,0);rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!   write(*,*)filename
!   call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
!   
!   ! graph 3
!   filename='graph_3'
!   D(1,1)=0.0;D(1,2)=1./2.;D(1,3)=1./2.
!   D(2,1)=0.0;D(2,2)=0.0;D(2,3)=1.
!   D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
!   rho(:)=cmplx(0,0);rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!   write(*,*)filename
!   call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
!   
!   ! graph 4
!   filename='graph_4'
!   D(1,1)=0.0;D(1,2)=1.;D(1,3)=0.0
!   D(2,1)=0.0;D(2,2)=0.0;D(2,3)=1.
!   D(3,1)=1.;D(3,2)=0.0;D(3,3)=0.0
!   rho(:)=cmplx(0,0);rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!   write(*,*)filename
!   call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
!   
!   ! graph 5
!   filename='graph_5'
!   D(1,1)=0.0;D(1,2)=1./2.;D(1,3)=1./2.
!   D(2,1)=0.0;D(2,2)=0.0;D(2,3)=0.0
!   D(3,1)=0.0;D(3,2)=1.;D(3,3)=0.0
!   rho(:)=cmplx(0,0);rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!   write(*,*)filename
!   call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
  
  

end program
  