!code to build the Superoperator SO from the hamiltonian H and using damping terms. 
! Both function require that SO, H, D are squar matrixs and SO is n^2 in size withn being the size of both H and D

module supperoperator
use until
use solver
use Densutil
use func
implicit none

 Integer, parameter   :: kdp = selected_real_kind(15)
 private              :: kdp
 contains
!make SO for the Hamiltonian only 
subroutine B_SO_for_H(H,SO)
 complex(kdp), dimension(:,:), intent(in)     :: H
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n


 n=size(H,1)
 SO(:,:)=cmplx(0,0)
 

!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1 
   do j=0,n-1
    do k=1,n
     do l=1,n
        if(i.eq.j)then
	    if(k.eq.l)then
             SO(i*n+k,j*n+l)=cmplx(0,1)*(H(i+1,i+1)-H(k,k))
            elseif(k.lt.l)then
             SO(i*n+k,j*n+l)=-cmplx(0,1)*H(l,k)
            else
	     SO(i*n+k,j*n+l)=-cmplx(0,1)*H(l,k)
	    endif 
        elseif((i.ne.j).and.(k.eq.l))then
            if((i.ne.1).and.(j.eq.1))then
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            elseif((i.eq.1).and.(j.ne.1))then
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            else
              SO(i*n+k,j*n+l)=cmplx(0,1)*H(i+1,j+1)
            endif
        endif
     enddo
    enddo
   enddo
 enddo



end subroutine


!adding in dephasing for the operator style <O|rho|O>
subroutine A_D_to_SO(D,SO)
 complex(kdp), dimension(:,:), intent(in)     :: D
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n


 n=size(D,1)
!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1
   do j=0,n-1
    do k=1,n
     do l=1,n
      if((i.eq.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,i+1)
      elseif((i.lt.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,j+1)
      elseif((i.gt.j).and.(k.eq.l))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,j+1)
      endif
     enddo
    enddo
   enddo
 enddo

end subroutine


! depasing for enery conserving 
subroutine De_en_cojn(D,SO)
 complex(kdp), dimension(:,:), intent(in)     :: D
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n,q

 
 n=size(D,1)
!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1
   do j=0,n-1
    do k=1,n
     do l=1,n
      if((i.eq.j).and.(k.eq.l).and.((i+1).ne.k))then
        do  q=1,n
	    SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-0.5*D(i+1,q)-0.5*D(k,q)
        enddo
      endif
     enddo
    enddo
   enddo
 enddo
 

end subroutine

! depasing for sink 
subroutine De_sink(sink,conect,strength,SO)
 real(kdp), intent(in)                        :: strength
 integer, intent(in)                          :: sink
 integer, dimension(:),intent(in)             :: conect
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,n,m,s
 
 m=size(SO,1)
 n=size(conect)
 s=int(sqrt(real(m)))
 
 do i=1,n
  SO(1+(sink-1)*(s+1),1+(conect(i)-1)*(s+1))=SO(1+(sink-1)*(s+1),1+(conect(i)-1)*(s+1))+strength
  SO(1+(conect(i)-1)*(s+1),1+(conect(i)-1)*(s+1))=SO(1+(conect(i)-1)*(s+1),1+(conect(i)-1)*(s+1))-strength
 enddo


end subroutine


! adds the dephasing for when it has the form: Ot rho O -1/2{Ot O , rho}
subroutine Stand_to_SO(D,SO)
 complex(kdp), dimension(:,:), intent(in)     :: D
 complex(kdp), dimension(:,:), intent(inout)  :: SO
 integer                                      :: i,j,k,l,q,n


 n=size(D,1)
!SO can be seen as being made up of n smaller matrixes of n by n size
!i and j move along the blocks and k,l move with in the blocks
 do i=0,n-1
   do j=0,n-1
    do k=1,n
     do l=1,n
      if((i.eq.j).and.(k.eq.l).and.((i+1).eq.k))then
        do  q=1,n
          if(q.ne.(i+1))then
	    SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-D(i+1,q)
          endif
        enddo
      elseif((i.eq.j).and.(k.eq.l).and.((i+1).ne.k))then
        do  q=1,n
	    SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)-0.5*D(i+1,q)-0.5*D(k,q)
        enddo
      elseif((i.lt.j).and.(k.eq.(i+1)).and.(l.eq.(j+1)))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)+D(l,k)
      elseif((i.gt.j).and.(k.eq.(i+1)).and.(l.eq.(j+1)))then
       SO(i*n+k,j*n+l)=SO(i*n+k,j*n+l)+D(l,k)
      endif
     enddo
    enddo
   enddo
 enddo

end subroutine
 
end module