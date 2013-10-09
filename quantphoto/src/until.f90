module until
 use supperoperator
 implicit none
 contains



! routine to take a transition matrix and make the full super operator
! this will generate the H matrix so it dose not need to be put in
 subroutine L_make_DG(SO,D,alpha)
  integer, parameter :: EP = selected_real_kind(15)
  complex(8), dimension(:,:), intent(in)     :: D
  complex(8), dimension(:,:), intent(inout)  :: SO
  real(EP),intent(in)                        :: alpha
  integer                                    :: i,j,k,l,n
  complex(EP), dimension(:,:), allocatable   :: H


  n=size(D,1)
  allocate(H(n,n))
   
  H(:,:)=cmplx(0,0)


  do i=1,n
    do j=1,n
     if(D(i,j).ne.0)then
       H(i,j)=1
       H(j,i)=1
     end if
    enddo
  enddo
  
  call B_SO_for_H((1-alpha)*H,SO)
  call Stand_to_SO(alpha*D,SO)

  deallocate(H)

 end subroutine






end module