! small program to perform a directed quantum walk on 3 node graphs.

program photest
  use eigenSolve
  use supperoperator
  use solver
  use until
  use Densutil
  use func
  use class
  use directwalk
  implicit none

 
  integer, parameter :: EP = selected_real_kind(15)
  integer                                      :: i, j,k, n,t,l
  real(EP)                                     :: alpha, err, beta
  complex(EP)                                  :: interact
  complex(EP), dimension(:,:), allocatable     :: D,D2
  complex(EP), dimension(:), allocatable       :: rho, p, rho2, rhod,TE1,TE2,TE3
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,AALPHA,filenamebase,dirnamebase,filename2
  logical                                      :: v, r, work
  
  
  
 CALL GET_COMMAND_ARGUMENT(1,infile)
 CALL GET_COMMAND_ARGUMENT(2,filenamebase)
 CALL GET_COMMAND_ARGUMENT(3,dirnamebase)
 CALL GET_COMMAND_ARGUMENT(4,nn)
 CALL GET_COMMAND_ARGUMENT(5,time)
 read(nn,*)n
 read(time,*)k
 
!  allocate(TE1(4),TE2(4),TE3(16))
!  TE1(1)=1;TE1(2)=1;TE1(3)=0;TE1(4)=0;
!  TE2(1)=2;TE2(2)=3;TE2(3)=2;TE2(4)=2;
!  call k_product_vform(TE1,TE2,TE3)
!  write(filename,'(a)')'TE3'
!  call write_Vec_real(filename,TE3)
!  goto 99
 
 
 err=0.000001
!   n=3
   t=5
   alpha=0.8
  allocate(D(n,n),rho(n*n),p(n),D2(n,n),rho2(n*n),rhod(16*16))

  v=.false.
  r=.true.
! write(*,*)size(D,1)
  call Read_MatR(infile,D)

  do j=1,k
    write(dirname,'(a,I5.5)')trim(dirnamebase),j 
    call make_dir(dirname)
    if(j.ne.1)then
     call Rand_Weight(D)
    endif
    write(filename,'(3a)')trim(dirname),'/','D'
    call write_Mat(filename,D)
    do l=1,11
     alpha=real(l-1)*0.1
     rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
     do i=1,n
       rho(i+(i-1)*n)=cmplx(1./n,0.0)
     enddo
!    rho(1)=cmplx(1.0,0)
     write(filename,'(3a,F4.2)')trim(dirname),'/',trim(filenamebase),alpha

     call Dir_Gra_Run(D, rho, t, alpha, filename, v, r)
        
    enddo
     
    write(dirname,'(a,a,I5.5)')trim(dirnamebase),'interact',j 
    call make_dir(dirname)
    do l=1,11
     alpha=real(l-1)*0.1
     rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
     do i=1,n
       rho(i+(i-1)*n)=cmplx(1./n,0.0)
     enddo
!    rho(1)=cmplx(1.0,0)
     write(filename,'(4a,F4.2)')trim(dirname),'/',trim(filenamebase),'interact',alpha
 
     rho2=rho
     interact=cmplx(0.0,0.0)
 
     call Dir_2walker(D, rho,rho2, t, alpha, filename, r, interact)
    enddo
      
    D2=D

    call row_norm(D)

!     
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
    call extract_pointerS(rho, p)
    
    call col_norm(D2)
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive quantum row normalized'
    close(4)

    call Dir_Gra_Con(D, rho, err, alpha, filename, r)

    
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'itarive class'
    close(4)

    call pagerank_it(transpose(D), p, err, alpha, filename)

       
    open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'eigenvalue class'
    close(4)

     call pagerank_ei(transpose(D), p, alpha, filename, work)
     
         open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
    write(4,*)' '
    write(4,*)'------------------------------------------------------------------------------'
    write(4,*)'Just gamma matrix'
    close(4)
    
    do i=1,n
     p(i)=cmplx(1/sqrt(real(n)),0.0)
    enddo
     
    call Dir_justgamma(transpose(D), p, t, filename, r)
      
   enddo


 err=0.01
 n=16
 do l=1,11
     alpha=real(l-1)*0.1
     do k=1,11
       beta=real(l-1)*0.1
       rhod(:)=cmplx(0,0);
       do i=1,n
        rhod(i+(i-1)*n)=cmplx(1./n,0.0)
       enddo
       write(filename,'(4a,2F4.2)')trim(dirname),'/',trim(filenamebase),'dis',alpha,beta
  
       call Dis_Dir_8line(rhod, err, alpha, beta, filename, r)
       
    enddo
 enddo
 
 99 continue
 
end program
  