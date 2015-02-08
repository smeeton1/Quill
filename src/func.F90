module func
use eigenSolve
use until
!use special
implicit none


 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp

 contains
 
 
 ! computes the first moment for a time series
 subroutine first_moment(in, moment)
  complex(kdp), dimension(:), intent(in)    :: in
  complex(kdp),intent(out)                  :: moment
  complex(kdp), dimension(:), allocatable   :: out 
  integer(kdp)				    :: N, i
  complex(kdp)                              :: top, bottom
  
  N=size(in)
  allocate(out(N))
  call fftw_1D(in, out)
  
  do i=1,N
    top=top+out(i)*CMPLX(i,0)*CMPLX(i,0)
  enddo
  
  do i=1,N
    bottom=bottom+out(i)
  enddo
   
  moment=sqrt(top/bottom)
  
  deallocate(out)
 
 end subroutine
 
 
  ! computes the second moment for a time series
 subroutine second_moment(in, moment)
  complex(kdp), dimension(:), intent(in)    :: in
  complex(kdp),intent(out)                  :: moment
  complex(kdp), dimension(:), allocatable   :: out 
  integer(kdp)				    :: N, i
  complex(kdp)                              :: top, bottom
  
  N=size(in)
  allocate(out(N))
  call fftw_1D(in, out)
  
  do i=1,N
    top=top+out(i)*CMPLX(i,0)
  enddo
  
  do i=1,N
    bottom=bottom+out(i)
  enddo
   
  moment=top/bottom
  
  deallocate(out)
 
 end subroutine
 
 subroutine write_moments(filename,in)
  complex(kdp), dimension(:,:), intent(in)  :: in
  complex(kdp)                              :: moment1, moment2
  integer(kdp)				    :: n,n2, i
  character(80), intent(in)                 :: filename
  character(90)                             :: filename2
 
  n=size(in,1)
  n2=size(in,2)
  write(filename2,'(a,a)')trim(filename),'_moments'
  open(5,file=filename2,STATUS='replace',ACCESS='append',ACTION='write')
  write(5,'(a,a,a,a,a)')'Node',' ','first moment',' ','second moment'
  do i=1,n2
   call first_moment(in(1:n,i), moment1)
   call second_moment(in(1:n,i), moment2)
   write(5,'(I3,2F8.5)')i,real(moment1),real(moment2)
  enddo
  close(5)
 end subroutine
 
 subroutine Rand_Weight(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    do j=1,m
      if(D(i,j).ne.0)then
       D(i,j)=5*rand()
      endif
    enddo
  enddo
  
  
 end subroutine
 
 subroutine Mat_Gen(D,porb)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  real,intent(in)                            :: porb
  integer				     :: i,j,n,m
  n=size(D,1)
  m=size(D,2)
  
  D(:,:)=cmplx(0,0);
  do i=1,n
    do j=1,m
      if((rand().gt.porb).and.(i.ne.j))then
       D(i,j)=5*rand()
      endif
    enddo
  enddo
  
  
 end subroutine
 
 subroutine row_norm(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  complex(kdp)  			     :: tot
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    tot=cmplx(0.0,0.0)
    do j=1,m
       tot = tot + D(i,j)
    enddo
    if(tot.ne.cmplx(0.0,0.0))then
     do j=1,m
        D(i,J) = D(i,j)/tot
     enddo
    endif
  enddo
    
 end subroutine
 
  subroutine col_norm(D)
  complex(kdp),dimension(:,:),intent(inout)  :: D
  integer				     :: i,j,n,m
  complex(kdp)  			     :: tot
  n=size(D,1)
  m=size(D,2)
  
  do i=1,n
    tot=cmplx(0.0,0.0)
    do j=1,m
       tot = tot + D(j,i)
    enddo
    if(tot.ne.cmplx(0.0,0.0))then
     do j=1,m
        D(j,i) = D(j,i)/tot
     enddo
    endif
  enddo
    
 end subroutine
 
 subroutine L_make_1(D)
  complex(kdp), dimension(:,:), intent(inout)  :: D
  integer                                      :: i,j,n
  complex(kdp)                                 :: row, col

  n=size(D,1)
  row=cmplx(0.0,0.0)
  col=cmplx(0.0,0.0)
  do i=1,n
    do j=1,n
       row=row+D(i,j)
       col=col+D(j,i)
    enddo
    D(i,i)=D(i,i)+row-col
    row=cmplx(0.0,0.0)
    col=cmplx(0.0,0.0)
  enddo

 end subroutine
 
  subroutine class_lap(D)
  complex(kdp), dimension(:,:), intent(inout)  :: D
  integer                                      :: i,j,n
  complex(kdp)                                 :: row

  n=size(D,1)
  row=cmplx(0.0,0.0)
  do i=1,n
    do j=1,n
       row=row+D(j,i)
    enddo
    D(i,i)=D(i,i)-row
    row=cmplx(0.0,0.0)
  enddo

 end subroutine
 
 subroutine class_lap_row(D)
  complex(kdp), dimension(:,:), intent(inout)  :: D
  integer                                      :: i,j,n
  complex(kdp)                                 :: row

  n=size(D,1)
  row=cmplx(0.0,0.0)
  do i=1,n
    do j=1,n
       row=row+D(i,j)
    enddo
    D(i,i)=D(i,i)-row
    row=cmplx(0.0,0.0)
  enddo

 end subroutine
 
 subroutine vec_to_mat(A,B)
  complex(kdp),dimension(:,:),intent(inout)  :: B
  complex(kdp),dimension(:),intent(in)       :: A
  integer				     :: i,j,n
 
  n=size(B,1)
 
   do i =0,n-1
    do j =1,n
      B(i+1,j)=A(j+i*n)
    enddo
  enddo
 
 end subroutine
 
 subroutine mat_to_vec(A,B)
  complex(kdp),dimension(:),intent(inout)    :: B
  complex(kdp),dimension(:,:),intent(in)     :: A
  integer				     :: i,j,n,m
 
  n=size(A,1) 
  
  do i =0,n-1
    do j =1,n
      B(j+i*n)=A(i+1,j)
    enddo
  enddo
 
 end subroutine
 
 !dosethe kronicer product of A and B both must be square matriecs for the same size. C must be n*n 
 subroutine k_product(A,B,C)
  complex(kdp),dimension(:,:),intent(inout)  :: C
  complex(kdp),dimension(:,:),intent(in)     :: A,B
  integer				     :: i,j,k,l,n,m
  
  n=size(A,1)
  m=size(C,1)
  
  do i=1,n
   do j=1,n
    do k=1,n
     do l=1,n
      C(n*(i-1)+k,n*(j-1)+l)=A(i,j)*B(k,l)
     enddo
    enddo
   enddo
  enddo
  
 
 end subroutine
 
 subroutine k_product_vform(A,B,C)
  complex(kdp),dimension(:),intent(inout)    :: C
  complex(kdp),dimension(:),intent(in)       :: A,B
  complex(kdp),dimension(:,:),allocatable    :: D,E,F
  integer				     :: i,j,k,l,n,m,p,q
  
  n=size(A,1)
  m=size(C,1)
  p=int(sqrt(real(n)))
  q=int(sqrt(real(m)))
  
  allocate(D(p,p),E(p,p),F(q,q))
  
  call vec_to_mat(A,D)
  call vec_to_mat(B,E)
  call k_product(D,E,F)
  call mat_to_vec(F,C)
 
 end subroutine
 
 subroutine k_sum(A,B,C)
  complex(kdp),dimension(:,:),intent(inout)  :: C
  complex(kdp),dimension(:,:),intent(in)     :: A,B
  complex(kdp),dimension(:,:), allocatable   :: In,C1,C2
  integer				     :: i,j,k,l,n,m
  
  n=size(A,1)
  m=size(C,1)
  allocate(In(n,n),C1(m,m),C2(m,m))
  
  In(:,:)=cmplx(0,0)
  C1(:,:)=cmplx(0,0)
  C2(:,:)=cmplx(0,0)
  
  do i=1,n
    In(i,i)=cmplx(1,0)
  enddo
  
  call k_product(A,In,C1)
  call k_product(In,B,C2)
  
  C=C1+C2
 
  deallocate(In,C1,C2) 
 
 end subroutine
 
 subroutine par_traceA(A,B)
  complex(kdp),dimension(:,:),intent(inout)  :: B
  complex(kdp),dimension(:,:),intent(in)     :: A
  integer				     :: i,j,k,l,n,m
  
  n=size(B,1)
  m=size(A,1)
  B(:,:)=cmplx(0,0)
  
  do i=0,n-1
   do j=0,n-1
    do k=1,n
      B(i+1,j+1)=B(i+1,j+1)+A(i*n+k,j*n+k)
    enddo
   enddo
  enddo
  
 
 end subroutine
 
 subroutine par_traceA_vec(A,B)
  complex(kdp),dimension(:),intent(inout)    :: B
  complex(kdp),dimension(:),intent(in)       :: A
  complex(kdp),dimension(:,:),allocatable    :: C , D
  integer				     :: i,j,k,l,n,m,p,q
  
  m=size(B,1)
  n=size(A,1)
  p=int(sqrt(real(n)))
  q=int(sqrt(real(m)))
  allocate(C(p,p),D(q,q))


  call vec_to_mat(A,C)
  call par_traceA(C,D)
  call mat_to_vec(D,B)


  deallocate(C,D)
  
 end subroutine
 
 subroutine par_traceB(A,B)
  complex(kdp),dimension(:,:),intent(inout)  :: B
  complex(kdp),dimension(:,:),intent(in)     :: A
  integer				     :: i,j,k,l,n,m
  
  n=size(B,1)
  m=size(A,1)
  B(:,:)=cmplx(0,0)
  
  do i=1,n
   do j=1,n
    do k=0,n-1
      if(i.eq.j)then
       B(i,j)=B(i,j)+ A(i+k*n,j+(k*n))
      else if(i.gt.j)then
       B(i,j)=B(i,j)+ A((i-1)*n+k*n,(j-1)*n+(k*n))
      else
       B(i,j)=B(i,j)+ A((i-1)*n+k*n,j*n+(k*n))
      endif
    enddo
   enddo
  enddo
 
 end subroutine
 
 subroutine par_traceB_vec(A,B)
  complex(kdp),dimension(:),intent(inout)    :: B
  complex(kdp),dimension(:),intent(in)       :: A
  complex(kdp),dimension(:,:),allocatable    :: C , D
  integer				     :: i,j,k,l,n,m,p,q
  
  m=size(B,1)
  n=size(A,1)
  p=int(sqrt(real(n)))
  q=int(sqrt(real(m)))
  allocate(C(p,p),D(q,q))


  call vec_to_mat(A,C)
  call par_traceB(C,D)
  call mat_to_vec(D,B)


  deallocate(C,D)
  
 end subroutine
 
 complex(kdp) function VonNueE(A) result(ent)
  complex(kdp),dimension(:,:),intent(in)     :: A
  complex(8), dimension(:),allocatable       :: v
  integer				     :: i,n
  
  n=size(A,1)
  allocate(v(n))
  call eigen(A,v)
  ent=0
  do i=1,n
    if(v(i).ne.cmplx(0.0,0.0))then
     ent=ent + v(i)*log(v(i))
    endif
  enddo
  deallocate(v)
 
 end function
 
 complex(kdp) function VonNueE_vec(A) result(ent)
  complex(kdp),dimension(:),intent(in)       :: A
  complex(kdp),dimension(:,:),allocatable    :: B
  integer				     :: i,j,n,m
  
  n=size(A,1)
  m=int(sqrt(real(n)))
  allocate(B(m,m))
  call vec_to_mat(A,B)
  ent = VonNueE(B)
  deallocate(B)
 
 end function
 

 
 complex(kdp) function RelativE_vec(A) result(ent)
  complex(kdp),dimension(:),intent(in)       :: A
  complex(kdp),dimension(:,:),allocatable    :: B
  integer				     :: i,j,n,m
  
  n=size(A,1)
  m=int(sqrt(real(n)))
  allocate(B(m,m))
  call vec_to_mat(A,B)
  ent=cmplx(0.0,0.0)
  do i=1,m
   if(B(i,i).ne.cmplx(0.0,0.0))then
    ent=ent + B(i,i)*log(B(i,i))
   endif
  enddo
  ent = ent - VonNueE(B)
  deallocate(B)
  
 end function
 
end module