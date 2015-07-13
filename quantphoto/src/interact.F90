program interact
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
  integer                                      :: i, j, n,t,l
  real(EP)                                     :: alpha,k
  complex(EP)                                  :: interac
  complex(EP), dimension(:,:), allocatable     :: D
  complex(EP), dimension(:), allocatable       :: rho,  rho2
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,filenamebase,dirnamebase,filename2
  
  
  
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
 

!   n=3
   t=50
   alpha=0.8
  allocate(D(n,n),rho(n*n),rho2(n*n))
  interac=cmplx(k,0.0)
! write(*,*)size(D,1)
  call Read_MatR(infile,D)
   
  write(dirname,'(a,a,F4.2)')trim(dirnamebase),'interact',real(interac)
  call make_dir(dirname)

    
    do l=1,11
     alpha=real(l-1)*0.1
     rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
!      do i=1,n
!        rho(i+(i-1)*n)=cmplx(1./n,0.0)
!      enddo
    rho(1)=cmplx(1.0,0)
     write(filename,'(4a,F4.2)')trim(dirname),'/',trim(filenamebase),'interact',alpha
 
     rho2=rho

     call Dir_2walker(D, rho,rho2, t, alpha, filename, interac)
    enddo
      




 
 99 continue
 
end program
  