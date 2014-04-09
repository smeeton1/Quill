! small program to perform a directed quantum walk on 3 node graphs.

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  use Densutil
  use func
  use class
  implicit none

 
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n,t,l
  real(EP)                                     :: alpha, err
  complex(EP), dimension(:,:), allocatable     :: D,D2
  complex(EP), dimension(:), allocatable       :: rho,p,rho2
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,AALPHA,filenamebase,dirnamebase
  logical                                      :: v, r, work
  
  
  

 CALL GET_COMMAND_ARGUMENT(1,filenamebase)
 CALL GET_COMMAND_ARGUMENT(2,dirnamebase)
 CALL GET_COMMAND_ARGUMENT(3,nn)
 CALL GET_COMMAND_ARGUMENT(4,time)
 read(nn,*)n
 read(time,*)t
 err=0.000001
 allocate(D(n,n),rho(n*n),p(n),D2(3,3),rho2(9))
 k=100
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
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
    write(filename,'(3a)')trim(dirname),'/',trim(filenamebase)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive quantum non row normalized'
    close(4)
    do l=1,10
        rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
         rho(i+(i-1)*n)=cmplx(1./n,0.0)
      enddo
      alpha=0.1*real(l)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'******************************************************************************'
      write(4,*)'alpha = ',alpha
      close(4)
      call Dir_Gra_Con(D, rho, err, alpha, filename, r)
    enddo
        
    call row_norm(D)
    
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
    call extract_pointerS(rho, p)
    
    
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive quantum row normalized'
    close(4)
    do l=1,10
        rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
        do i=1,n
          rho(i+(i-1)*n)=cmplx(1./n,0.0)
        enddo
      alpha=0.1*real(l)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'******************************************************************************'
      write(4,*)'alpha = ',alpha
      close(4)
      call Dir_Gra_Con(D, rho, err, alpha, filename, r)
    enddo
    
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive class'
    close(4)
    do l=1,10
      alpha=0.1*real(l)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'******************************************************************************'
      write(4,*)'alpha = ',alpha
      close(4)
      call pagerank_it(D, p, err, alpha, filename)
    enddo
       
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'eigenvalue class'
    close(4)
    do l=1,10
      alpha=0.1*real(l)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'******************************************************************************'
      write(4,*)'alpha = ',alpha
      close(4)
      call pagerank_ei(D, p, alpha, filename, work)
    enddo
   enddo
   
   !3 node graphs
   write(dirnamebase,'(3a)')trim('results'),'/',trim('3graphs')
   call make_dir(dirnamebase)
   
   !graph 1
   write(dirname,'(3a)')trim(dirnamebase),'/',trim('graph_1')
   call make_dir(dirname)
   D2(1,1)=0;D2(1,2)=1;D2(1,3)=0
   D2(1,1)=0;D2(1,2)=0;D2(1,3)=1
   D2(1,1)=1;D2(1,2)=0;D2(1,3)=0
   write(filename,'(3a)')trim(dirname),'/','D'
   call write_Mat(filename,D2) 
   do l=1,10
      rho2(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
          rho2(i+(i-1)*3)=cmplx(1./3.,0.0)
      enddo
      alpha=0.1*real(l)
      write(filename,'(3a,I5.5)')trim(dirname),'/',trim(filenamebase),l
      call Dir_Gra_Run(D2, rho2,100, alpha, filename, v, r)
    enddo
    
    
    !graph 2
   write(dirname,'(3a)')trim(dirnamebase),'/',trim('graph_2')
   call make_dir(dirname)
   D2(1,1)=0;D2(1,2)=1;D2(1,3)=0
   D2(1,1)=1;D2(1,2)=0;D2(1,3)=1
   D2(1,1)=1;D2(1,2)=1;D2(1,3)=0
   write(filename,'(3a)')trim(dirname),'/','D'
   call write_Mat(filename,D2) 
   do l=1,10
      rho2(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
          rho2(i+(i-1)*3)=cmplx(1./3.,0.0)
      enddo
      alpha=0.1*real(l)
      write(filename,'(3a,I5.5)')trim(dirname),'/',trim(filenamebase),l
      call Dir_Gra_Run(D2, rho2,100, alpha, filename, v, r)
    enddo
    
    !graph 3
   write(dirname,'(3a)')trim(dirnamebase),'/',trim('graph_3')
   call make_dir(dirname)
   D2(1,1)=0;D2(1,2)=1;D2(1,3)=1
   D2(1,1)=1;D2(1,2)=0;D2(1,3)=1
   D2(1,1)=0;D2(1,2)=0;D2(1,3)=0
   write(filename,'(3a)')trim(dirname),'/','D'
   call write_Mat(filename,D2) 
   do l=1,10
      rho2(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
          rho2(i+(i-1)*3)=cmplx(1./3.,0.0)
      enddo
      alpha=0.1*real(l)
      write(filename,'(3a,I5.5)')trim(dirname),'/',trim(filenamebase),l
      call Dir_Gra_Run(D2, rho2,100, alpha, filename, v, r)
    enddo
    
    !graph 4
   write(dirname,'(3a)')trim(dirnamebase),'/',trim('graph_4')
   call make_dir(dirname)
   D2(1,1)=0;D2(1,2)=1;D2(1,3)=0
   D2(1,1)=1;D2(1,2)=0;D2(1,3)=1
   D2(1,1)=1;D2(1,2)=0;D2(1,3)=0
   write(filename,'(3a)')trim(dirname),'/','D'
   call write_Mat(filename,D2) 
   do l=1,10
      rho2(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
          rho2(i+(i-1)*3)=cmplx(1./3.,0.0)
      enddo
      alpha=0.1*real(l)
      write(filename,'(3a,I5.5)')trim(dirname),'/',trim(filenamebase),l
      call Dir_Gra_Run(D2, rho2,100, alpha, filename, v, r)
    enddo
    
    !graph 5
   write(dirname,'(3a)')trim(dirnamebase),'/',trim('graph_5')
   call make_dir(dirname)
   D2(1,1)=0;D2(1,2)=1;D2(1,3)=0
   D2(1,1)=0;D2(1,2)=0;D2(1,3)=0
   D2(1,1)=1;D2(1,2)=1;D2(1,3)=0
   write(filename,'(3a)')trim(dirname),'/','D'
   call write_Mat(filename,D2) 
   do l=1,10
      rho2(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
      do i=1,n
          rho2(i+(i-1)*3)=cmplx(1./3.,0.0)
      enddo
      alpha=0.1*real(l)
      write(filename,'(3a,I5.5)')trim(dirname),'/',trim(filenamebase),l
      call Dir_Gra_Run(D2, rho2,100, alpha, filename, v, r)
    enddo
   
   
   
   

end program
  