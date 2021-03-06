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
  integer                                      :: i, j,k, n,t,l,start
  real(EP)                                     :: alpha, err, beta
  complex(EP)                                  :: interact
  complex(EP), dimension(:,:), allocatable     :: D,D2
  complex(EP), dimension(:), allocatable       :: rho, p, rho2, rhod,TE1,TE2,TE3
  character(len=80)                            :: fmt,filename,dirname,infile,time,nn,AALPHA,filenamebase,dirnamebase,filename2,init
  logical                                      :: v, r, work,suppos
  
  
  
 CALL GET_COMMAND_ARGUMENT(1,infile)
 CALL GET_COMMAND_ARGUMENT(2,filenamebase)
 CALL GET_COMMAND_ARGUMENT(3,dirnamebase)
 CALL GET_COMMAND_ARGUMENT(4,nn)
 CALL GET_COMMAND_ARGUMENT(5,time)
 CALL GET_COMMAND_ARGUMENT(5,init)
 read(nn,*)n
 read(time,*)k
 if(init.eq.'S')then 
 suppos=.true.
 else
 read(time,*)start
 suppos=.false. 
 endif
 
 
!  allocate(TE1(9),TE2(9),TE3(81))
!  TE1(1)=1./3.;TE1(2)=0;TE1(3)=0;TE1(4)=0;TE1(5)=1./3.;TE1(6)=0;TE1(7)=0;TE1(8)=0;TE1(9)=1./3.;
!  TE2(1)=1./3.;TE2(2)=0;TE2(3)=0;TE2(4)=0.;TE2(5)=1./3.;TE2(6)=0;TE2(7)=0;TE2(8)=0;TE2(9)=1./3.;
!  call k_product_vform(TE1,TE2,TE3)
!  call par_traceA_vec(TE3,TE1)
!  write(filename,'(a)')'TE3'
!  call write_Vec_real(filename,TE1)
!  goto 99
 
 
 err=0.00001
!   n=3
   t=200
   alpha=0.8
  allocate(D(n,n),rho(n*n),p(n),D2(n,n),rho2(n*n))

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
    write(filename,'(3a,F4.2)')trim(dirname),'/',trim(filenamebase)
    do l=1,1!1
     alpha=real(l-1)*0.1
     rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
     if(suppos)then
       do i=1,n
          rho(i+(i-1)*n)=cmplx(1./n,0.0)
       enddo
     else
       rho(start+(start-1)*n)=cmplx(1.,0.0)
     endif
!      rho(1)=cmplx(1.0,0)
!     write(filename,'(3a,F4.2)')trim(dirname),'/',trim(filenamebase),alpha
     open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
     write(4,*)' '
     write(4,*)'------------------------------------------------------------------------------'
     write(4,*)'alpha',alpha
     close(4)
     call Dir_Gra_Con(D, rho, err, alpha, filename, r)
!     call Dir_Gra_Run(D, rho, t, 0.1, alpha, filename, v, r)
        
    enddo
    if(suppos)then

!      rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
! !      do i=1,n
! !        rho(i+(i-1)*n)=cmplx(1./n,0.0)
! !      enddo
!      rho(1)=cmplx(1.0,0)
!      write(filename,'(3a,F4.2)')trim(dirname),'/',trim(filenamebase),1.0
      alpha=0.8
!      call Dir_Gra_Run(D, rho, t, 0.1, alpha, filename, v, r)
      
!     D2=D
! 
    call row_norm(D)
! 
! !     
    rho(:)=cmplx(0,0);!rho(1)=cmplx(1./3.,0);rho(5)=cmplx(1./3.,0);rho(9)=cmplx(1./3.,0)
    do i=1,n
      rho(i+(i-1)*n)=cmplx(1./n,0.0)
    enddo
    call extract_pointerS(rho, p)
!     
!     call col_norm(D2)
!     open(4,file=filename,STATUS='unknown',ACCESS='append',ACTION='write')
!     write(4,*)' '
!     write(4,*)'------------------------------------------------------------------------------'
!     write(4,*)'itarive quantum row normalized'
!     close(4)
! 
!     call Dir_Gra_Con(D, rho, err, alpha, filename, r)
! 
!   

     write(filename,'(4a,F4.2)')trim(dirname),'/',trim(filenamebase),'class',alpha
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
     
    endif

    enddo  
    
!     D2=D
! !     do i=1,n
! !      p(i)=cmplx(1/sqrt(real(n)),0.0)
! !     enddo
!     p(n)=cmplx(1.0,0.0)
!     call class_lap(D2)
!     write(filename,'(3a)')trim(dirname),'/','justgamma' 
!     call Dir_justgamma(D2, p, t, 0.1, filename, r, .true.)
    
!     call row_norm(D)
!     write(filename,'(3a)')trim(dirname),'/','justgammarownorm' 
!     call Dir_justgamma(transpose(D), p, t, filename, r)
!     
!     call L_make_1(D)
!     write(filename,'(3a)')trim(dirname),'/','justgammarownormL' 
!     call Dir_justgamma(transpose(D), p, t, filename, r)
!     
!     call L_make_1(D2)
!     write(filename,'(3a)')trim(dirname),'/','justgammaL' 
!     call Dir_justgamma(transpose(D2), p, t, filename, r)

 
 99 continue
 
end program
  