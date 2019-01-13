!!=========================================================
!! blockinterface_average
!!
!! purpose: dealing with block interface 
!!          with a averaged procedure
!!
!! AUTHOR: Jiamin XU                                       
!! DATE:   2012.11.28
!!=========================================================
subroutine UpdateBufferCoordinate2d
	use blk_var
	use index_var
	use mpi_var
	use glbindex_var
	use flag_var
	implicit none
	
	integer :: i,j,k
	integer :: nBuffer
	integer :: dim_f
	
	integer :: n1, n2, n3
	integer :: n11,n12,n13
	integer :: nn
	integer :: is_bc,  ie_bc,  js_bc,  je_bc,  ks_bc,  ke_bc
	integer :: is_bc_t,ie_bc_t,js_bc_t,je_bc_t,ks_bc_t,ke_bc_t
	integer :: blk_t,subface_t,face_t
	!!index
	integer :: l1, m1
	integer :: ll0,mm0
	integer :: ll, mm
	integer :: l2, m2
	integer :: face,ori
	
	real*8,dimension(:,:,:),allocatable :: ua,ua2
	
	integer :: m00
	integer :: itag
	integer :: ireq(  100, nblock)
	integer :: ista(1,100, nblock)
	
	!!
	integer :: i0,j0,k0

	!!load data to buffer blk
	do m0 = 1,blk_loop
		do ksub = 1,blk(m0)%num_subface
			blk_t = blk(m0)%bc(ksub)%blk_t
			if (blk_t .ge. 0) then
				face = blk(m0)%bc(ksub)%face
				if     (face .eq. 1)then	!!face i-
					is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
					js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
					
					k = 1
					do j = js_bc,je_bc
						do nBuffer = 1, bufferLength
							blk(m0)%bc(ksub)%buffer_send(nBuffer,1,is_bc,j,k) = blk(m0)%x(is_bc+nBuffer,j,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,2,is_bc,j,k) = blk(m0)%y(is_bc+nBuffer,j,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,3,is_bc,j,k) = blk(m0)%z(is_bc+nBuffer,j,k)
						end do 
					end do
				else if(face .eq. 2)then	!!face i+
					is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
					js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
					
					k = 1
					do j = js_bc,je_bc
						do nBuffer = 1, bufferLength
							blk(m0)%bc(ksub)%buffer_send(nBuffer,1,ie_bc,j,k) = blk(m0)%x(ie_bc-nBuffer,j,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,2,ie_bc,j,k) = blk(m0)%y(ie_bc-nBuffer,j,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,3,ie_bc,j,k) = blk(m0)%z(ie_bc-nBuffer,j,k)
						end do
					end do					
				else if(face .eq. 3)then	!!face j-
					is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
					js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
					
					k = 1
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							blk(m0)%bc(ksub)%buffer_send(nBuffer,1,i,js_bc,k) = blk(m0)%x(i,js_bc+nBuffer,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,2,i,js_bc,k) = blk(m0)%y(i,js_bc+nBuffer,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,3,i,js_bc,k) = blk(m0)%z(i,js_bc+nBuffer,k)
						end do 
					end do					
				else if(face .eq. 4)then	!!face j+
					is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
					js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
					
					k = 1
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							blk(m0)%bc(ksub)%buffer_send(nBuffer,1,i,je_bc,k) = blk(m0)%x(i,je_bc-nBuffer,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,2,i,je_bc,k) = blk(m0)%y(i,je_bc-nBuffer,k)
							blk(m0)%bc(ksub)%buffer_send(nBuffer,3,i,je_bc,k) = blk(m0)%z(i,je_bc-nBuffer,k)
						end do
					end do
				end if
				!!*						
			end if
		end do
	end do
	!!*end load data to buffer blk
	
	!!send buffer data
	do m0 = 1,blk_loop                                 !!start: m0
		call m02m(m,m0,myid,myn)                       !!start: m
		do ksub = 1,blk(m0)%num_subface                !!start: ksub
			blk_t     = blk(m0)%bc(ksub)%blk_t         !!end:   m
			subface_t = blk(m0)%bc(ksub)%subface_t     !!end:   ksub
			dest_id   = blk(m0)%bc(ksub)%rank_t        !!end:   rank
			if (blk_t .ge. 0) then
				dim_f     = bufferLength*5*blk(m0)%bc(ksub)%dim_f
				itag      = 100000*blk(m0)%bc(ksub)%f_no + m
				if     (dest_id .ne. myid)then
					call MPI_isend(blk(m0)%bc(ksub)%buffer_send,dim_f,       &
					               MPI_REAL8,dest_id,itag,MPI_COMM_WORLD,ireq(ksub,m),ierr)
				else if(dest_id .eq. myid)then  
					is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
					js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
					ks_bc = blk(m0)%bc(ksub)%ks; ke_bc = blk(m0)%bc(ksub)%ke
					m00 = blk_t - myid*myn
					do k = ks_bc,ke_bc
						do j = js_bc,je_bc
							do i = is_bc,ie_bc
								do nBuffer = 1, bufferLength
									blk(m00)%bc(subface_t)%buffer_recv(nBuffer,1,i,j,k) = blk(m0)%bc(ksub)%buffer_send(nBuffer,1,i,j,k)
									blk(m00)%bc(subface_t)%buffer_recv(nBuffer,2,i,j,k) = blk(m0)%bc(ksub)%buffer_send(nBuffer,2,i,j,k)
									blk(m00)%bc(subface_t)%buffer_recv(nBuffer,3,i,j,k) = blk(m0)%bc(ksub)%buffer_send(nBuffer,3,i,j,k)
								end do 							
							end do
						end do
					end do
				end if
			end if
		end do
	end do
	!!*	
	
	!!receive data from buffer block
	do m0 = 1,blk_loop
		call m02m(m,m0,myid,myn)
		do ksub = 1,blk(m0)%num_subface
			blk_t     = blk(m0)%bc(ksub)%blk_t
			subface_t = blk(m0)%bc(ksub)%subface_t
			source_id = blk(m0)%bc(ksub)%rank_t
			if (blk_t .ge. 0) then
				if (source_id .ne. myid)then
					dim_f     = bufferLength*5*blk(m0)%bc(ksub)%dim_f
					itag      = 100000*subface_t + blk_t
					m00       = blk_t - source_id*myn
					call MPI_irecv(blk(m0)%bc(ksub)%buffer_recv,dim_f,       &
					               MPI_REAL8,source_id,itag,MPI_COMM_WORLD,ireq(subface_t,blk_t),ierr)
				end if
			end if
		end do
	end do
	!!*   
    
	!!mpi wait
	!!send buffer data
	do m0 = 1,blk_loop
		call m02m(m,m0,myid,myn)
		do ksub = 1,blk(m0)%num_subface
			blk_t     = blk(m0)%bc(ksub)%blk_t
			subface_t = blk(m0)%bc(ksub)%subface_t
			dest_id   = blk(m0)%bc(ksub)%rank_t
			if (blk_t .ge. 0) then
				if (dest_id .ne. myid)then
					call MPI_wait(ireq(ksub,m),ista(1,ksub,m),ierr)	
				end if			
			end if
		end do
	end do
	!!*
  
	!!receive data from buffer block
	do m0 = 1,blk_loop
		call m02m(m,m0,myid,myn)
		do ksub = 1,blk(m0)%num_subface
			blk_t     = blk(m0)%bc(ksub)%blk_t
			subface_t = blk(m0)%bc(ksub)%subface_t
			dest_id   = blk(m0)%bc(ksub)%rank_t
			if (blk_t .ge. 0) then
				if (dest_id .ne. myid)then
					call MPI_wait(ireq(subface_t,blk_t),ista(1,subface_t,blk_t),ierr)	
				end if				
			end if
		end do
	end do
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*	

	!!average procedure
	!!load buffer data to blk
	do m0 = 1,blk_loop
		call m02m(m,m0,myid,myn)
		do ksub=1,blk(m0)%num_subface
			!!current blk and subface info
			blk_t     = blk(m0)%bc(ksub)%blk_t
			subface_t = blk(m0)%bc(ksub)%subface_t
			is_bc     = blk(m0)%bc(ksub)%is
			ie_bc     = blk(m0)%bc(ksub)%ie
			js_bc     = blk(m0)%bc(ksub)%js
			je_bc     = blk(m0)%bc(ksub)%je
			
			n1 = ie_bc - is_bc + 1
			n2 = je_bc - js_bc + 1
			
			if (blk_t .ge. 0) then
				!!neighb subface dimension
				face_t  = blk0(blk_t)%bc0(subface_t)%face
				
				is_bc_t = blk0(blk_t)%bc0(subface_t)%is
				ie_bc_t = blk0(blk_t)%bc0(subface_t)%ie
				js_bc_t = blk0(blk_t)%bc0(subface_t)%js
				je_bc_t = blk0(blk_t)%bc0(subface_t)%je

				n11 = ie_bc_t - is_bc_t + 1
				n12 = je_bc_t - js_bc_t + 1
				nn  = max(n11,n12)
				
				!!load buffer_recv to ua
				allocate(ua(bufferLength,3,nn))
				if      (face_t .eq. 1)then	!!boundary  i-
					k = 1
					do ll0 = 1,nn
						ua(:,1,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t,js_bc_t+ll0-1,k)
						ua(:,2,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t,js_bc_t+ll0-1,k)
						ua(:,3,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t,js_bc_t+ll0-1,k)
					end do
				else if (face_t .eq. 2)then	!!boundary  i+
					k = 1
					do ll0 = 1,nn
						ua(:,1,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,1,ie_bc_t,js_bc_t+ll0-1,k)
						ua(:,2,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,2,ie_bc_t,js_bc_t+ll0-1,k)
						ua(:,3,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,3,ie_bc_t,js_bc_t+ll0-1,k)
					end do
				else if (face_t .eq. 3)then	!!boundary  j-
					k = 1
					do ll0 = 1,nn
						ua(:,1,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+ll0-1,js_bc_t,k)
						ua(:,2,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+ll0-1,js_bc_t,k)
						ua(:,3,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+ll0-1,js_bc_t,k)
					end do
				else if (face_t .eq. 4)then	!!boundary  j+
					k = 1
					do ll0 = 1,nn	
						ua(:,1,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+ll0-1,je_bc_t,k)
						ua(:,2,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+ll0-1,je_bc_t,k)
						ua(:,3,ll0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+ll0-1,je_bc_t,k)
					end do
				end if
				!!*
				!!end update ua				
				
				ori = blk(m0)%bc(ksub)%ori
				allocate(ua2(bufferLength,3,nn))
				if      (ori .eq. 2)then
					do ll0 = 1,nn
						ua2(:,:,ll0) = ua(:,:,ll0)
					end do
				else if (ori .eq. 3)then
					do ll0 = 1,nn
						ua2(:,:,ll0) = ua(:,:,nn-ll0+1)
					end do
				end if
				!!*
				
				face = blk(m0)%bc(ksub)%face
				!!i- face
				if (face .eq. 1) then !!boundary  i- 
					i0 = is_bc
					k0 = 1
					do j=1,nn
						j0 = js_bc+j-1
						do nBuffer = 1, bufferLength
							blk(m0)%x(i0-nBuffer,j0,k0) = ua2(nBuffer,1,j)
							blk(m0)%y(i0-nBuffer,j0,k0) = ua2(nBuffer,2,j)
							blk(m0)%z(i0-nBuffer,j0,k0) = ua2(nBuffer,3,j)
						end do 
					end do
				end if
				!!*
				
				!!i+ face
				if (face .eq. 2) then !!boundary  i+
					i0 = ie_bc
					k0 = 1
					do j=1,nn
						j0 = js_bc+j-1
						do nBuffer = 1, bufferLength
							blk(m0)%x(i0+nBuffer,j0,k0) = ua2(nBuffer,1,j)
							blk(m0)%y(i0+nBuffer,j0,k0) = ua2(nBuffer,2,j)
							blk(m0)%z(i0+nBuffer,j0,k0) = ua2(nBuffer,3,j)
						end do 
					end do
				end if
				!!*
				
				!!j- face
				if (face .eq. 3) then	!!boundary  j-
					j0 = js_bc
					k0 = 1
					do i=1,nn
						i0 = is_bc+i-1
						do nBuffer = 1, bufferLength
							blk(m0)%x(i0,j0-nBuffer,k0) = ua2(nBuffer,1,i)
							blk(m0)%y(i0,j0-nBuffer,k0) = ua2(nBuffer,2,i)
							blk(m0)%z(i0,j0-nBuffer,k0) = ua2(nBuffer,3,i)
						end do 
					end do
				end if
				!!*
				
				!!j+ face
				if (face .eq. 4) then	!!boundary  j+
					j0 = je_bc
					k0 = 1
					do i=1,nn
						i0 = is_bc+i-1
						do nBuffer = 1, bufferLength
							blk(m0)%x(i0,j0+nBuffer,k0) = ua2(nBuffer,1,i) 
							blk(m0)%y(i0,j0+nBuffer,k0) = ua2(nBuffer,2,i)
							blk(m0)%z(i0,j0+nBuffer,k0) = ua2(nBuffer,3,i)
						end do
					end do
				end if
				!!*
				deallocate(ua,ua2)
			endif
		enddo
	enddo
	
	return
end subroutine
