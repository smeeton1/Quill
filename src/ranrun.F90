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
  complex(EP), dimension(:), allocatable       :: rho,p
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,AALPHA,filenamebase,dirnamebase
  logical                                      :: v, r, work
  
  
  

 CALL GET_COMMAND_ARGUMENT(1,filenamebase)
 CALL GET_COMMAND_ARGUMENT(2,dirnamebase)
 CALL GET_COMMAND_ARGUMENT(3,nn)
 CALL GET_COMMAND_ARGUMENT(4,time)
 read(nn,*)n
 read(time,*)k
 err=0.000001
 allocate(D(n,n),rho(n*n),p(n),D2(n,n))

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
!       open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
!       write(4,*)' '
!       write(4,*)'--------------------------------------------------------------------------------'
!       write(4,*)'-D'
!       close(4)
!       call Dir_Gra_Con(-D, rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose D'
      close(4)
      call Dir_Gra_Con(Transpose(D), rho, err, alpha, filename, r)
!       open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
!       write(4,*)' '
!       write(4,*)'--------------------------------------------------------------------------------'
!       write(4,*)'Transpose -D'
!       close(4)
!       call Dir_Gra_Con(Transpose(-D), rho, err, alpha, filename, r)
    enddo
    D2=D    
    call row_norm(D)
    
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
    call extract_pointerS(rho, p)
    
    call col_norm(D2)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive quantum normalized'
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
      call Dir_Gra_Con(D2, rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'row D'
      close(4)
      call Dir_Gra_Con(D, rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose D'
      close(4)
      call Dir_Gra_Con(Transpose(D2), rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose row D'
      close(4)
      call Dir_Gra_Con(Transpose(D), rho, err, alpha, filename, r)
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
       call pagerank_it(D2, p, err, alpha, filename)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'row D'
      close(4)
     call Dir_Gra_Con(D, rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose D'
      close(4)
      call Dir_Gra_Con(Transpose(D2), rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose row D'
      close(4)
      call Dir_Gra_Con(Transpose(D), rho, err, alpha, filename, r)
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
      call pagerank_ei(D2, p, alpha, filename, work)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'row D'
      close(4)
      call Dir_Gra_Con(D, rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose D'
      close(4)
      call Dir_Gra_Con(Transpose(D2), rho, err, alpha, filename, r)
      open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
      write(4,*)' '
      write(4,*)'--------------------------------------------------------------------------------'
      write(4,*)'Transpose row D'
      close(4)
      call Dir_Gra_Con(Transpose(D), rho, err, alpha, filename, r)
    enddo
   enddo

   
   
   
   

end program
  