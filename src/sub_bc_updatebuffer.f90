!!=========================================================
!! blockinterface_average
!!
!! purpose: dealing with block interface 
!!          with a averaged procedure
!!
!! AUTHOR: Jiamin XU                                       
!! DATE:   2012.11.28
!!=========================================================
subroutine UpdateBuffer
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
	integer :: ind
	
	real*8,dimension(:,:,:,:),allocatable :: ua,ua2
	
	integer :: m00
	integer :: itag
	integer :: ireq(  100, nblock)
	integer :: ista(1,100, nblock)
	
	!!
	integer :: i0,j0,k0
	real*8  :: q1(5),q2(5)

	!!load data to buffer blk
	do m0 = 1,blk_loop
		do ksub = 1,blk(m0)%num_subface
			blk_t = blk(m0)%bc(ksub)%blk_t
			if (blk_t .ge. 0) then
				face = blk(m0)%bc(ksub)%face

				is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
				js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
				ks_bc = blk(m0)%bc(ksub)%ks; ke_bc = blk(m0)%bc(ksub)%ke
				if     (face .eq. 1)then	!!face i-
					do k = ks_bc,ke_bc
					do j = js_bc,je_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,is_bc,j,k) = blk(m0)%q(ind,is_bc+nBuffer,j,k)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,is_bc,j,k) = blk(m0)%nut(is_bc+nBuffer,j,k)
							end if
						end do 
					end do
					end do
				else if(face .eq. 2)then	!!face i+
					do k = ks_bc,ke_bc
					do j = js_bc,je_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,ie_bc,j,k) = blk(m0)%q(ind,ie_bc-nBuffer,j,k)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,ie_bc,j,k) = blk(m0)%nut(ie_bc-nBuffer,j,k)
							end if
						end do
					end do		
					end do			
				else if(face .eq. 3)then	!!face j-
					do k = ks_bc,ke_bc
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,i,js_bc,k) = blk(m0)%q(ind,i,js_bc+nBuffer,k)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,i,js_bc,k) = blk(m0)%nut(i,js_bc+nBuffer,k)
							end if
						end do 
					end do
					end do					
				else if(face .eq. 4)then	!!face j+
					do k = ks_bc,ke_bc
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,i,je_bc,k) = blk(m0)%q(ind,i,je_bc-nBuffer,k)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,i,je_bc,k) = blk(m0)%nut(i,je_bc-nBuffer,k)
							end if
						end do
					end do
					end do
				else if(face .eq. 5)then	!!face k-
					do j = js_bc,je_bc
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,i,j,ks_bc) = blk(m0)%q(ind,i,j,ks_bc+nBuffer)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,i,j,ks_bc) = blk(m0)%nut(i,j,ks_bc+nBuffer)
							end if
						end do
					end do
					end do
				else if(face .eq. 6)then	!!face k-
					do j = js_bc,je_bc
					do i = is_bc,ie_bc
						do nBuffer = 1, bufferLength
							do ind = 1, 5
								blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,i,j,ke_bc) = blk(m0)%q(ind,i,j,ke_bc-nBuffer)
							end do
							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%bc(ksub)%buffer_send(nBuffer,6,i,j,ke_bc) = blk(m0)%nut(i,j,ke_bc-nBuffer)
							end if
						end do
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
				dim_f     = bufferLength*9*blk(m0)%bc(ksub)%dim_f
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
						do ind = 1,9
							blk(m00)%bc(subface_t)%buffer_recv(nBuffer,ind,i,j,k) = blk(m0)%bc(ksub)%buffer_send(nBuffer,ind,i,j,k)
						end do
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
					dim_f     = bufferLength*9*blk(m0)%bc(ksub)%dim_f
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
			ks_bc     = blk(m0)%bc(ksub)%ks
			ke_bc     = blk(m0)%bc(ksub)%ke
			
			n1 = ie_bc - is_bc + 1
			n2 = je_bc - js_bc + 1
			n3 = ke_bc - ks_bc + 1
			
			if (blk_t .ge. 0) then
				!!neighb subface dimension
				face_t  = blk0(blk_t)%bc0(subface_t)%face
				
				is_bc_t = blk0(blk_t)%bc0(subface_t)%is
				ie_bc_t = blk0(blk_t)%bc0(subface_t)%ie
				js_bc_t = blk0(blk_t)%bc0(subface_t)%js
				je_bc_t = blk0(blk_t)%bc0(subface_t)%je
				ks_bc_t = blk0(blk_t)%bc0(subface_t)%ks
				ke_bc_t = blk0(blk_t)%bc0(subface_t)%ke

				n11 = ie_bc_t - is_bc_t + 1
				n12 = je_bc_t - js_bc_t + 1
				n13 = ke_bc_t - ks_bc_t + 1

				if      (face_t .eq. 1)then !!boundary  i-
					l1 = n13	!!k ==> l
					m1 = n12    !!j ==> m
				else if (face_t .eq. 2)then !!boundary  i+
					l1 = n12	!!j ==> l
					m1 = n13	!!k ==> m
				else if (face_t .eq. 3)then !!boundary  j-
					l1 = n11	!!i ==> l
					m1 = n13	!!k ==> m
				else if (face_t .eq. 4)then !!boundary  j+
					l1 = n13	!!k ==> l
					m1 = n11	!!i ==> m
				else if (face_t .eq. 5)then !!boundary  k-
					l1 = n12	!!j ==> l
					m1 = n11	!!i ==> m
				else if (face_t .eq. 6)then !!boundary  k+
					l1 = n11	!!i ==> l
					m1 = n12	!!j ==> m
				end if
				!!*
			
				!!load buffer_recv to ua
				allocate(ua(bufferLength,9,l1,m1))
				if      (face_t .eq. 1)then	!!boundary  i-
					do mm0 = 1,m1
					do ll0 = 1,l1
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,is_bc_t,js_bc_t+mm0-1,ks_bc_t+ll0-1)
					end do
					end do
				else if (face_t .eq. 2)then	!!boundary  i+
					do mm0=1,m1
					do ll0=1,l1
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,ie_bc_t,js_bc_t+ll0-1,ks_bc_t+mm0-1)
					end do
					end do
				else if (face_t .eq. 3)then	!!boundary  j-
					do mm0=1,m1
					do ll0=1,l1
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,is_bc_t+ll0-1,js_bc_t,ks_bc_t+mm0-1)
					end do
					end do
				else if (face_t .eq. 4)then	!!boundary  j+
					do mm0=1,m1
					do ll0=1,l1		
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,is_bc_t+mm0-1,je_bc_t,ks_bc_t+ll0-1)
					end do
					end do
				else if (face_t .eq. 5) then
					do mm0=1,m1
					do ll0=1,l1		
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,is_bc_t+mm0-1,js_bc_t+ll0-1,ks_bc_t)
					end do
					end do
				else if (face_t .eq. 6) then
					do mm0=1,m1
					do ll0=1,l1		
						ua(:,1,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,1,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
						ua(:,2,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,2,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
						ua(:,3,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,3,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
						ua(:,4,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,4,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
						ua(:,5,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,5,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
						ua(:,6,ll0,mm0) = blk(m0)%bc(ksub)%buffer_recv(:,6,is_bc_t+ll0-1,js_bc_t+mm0-1,ke_bc_t)
					end do
					end do
				end if
				!!*
				!!end update ua				
				
				ori = blk(m0)%bc(ksub)%ori
				if      (ori .eq. 1)then
					l2 = l1
					m2 = m2
					allocate(ua2(bufferLength,9,l2,m2))
					do mm = 1,m2
					do ll = 1,l2
						ua2(:,:,ll,mm) = ua(:,:,ll,m1-mm+1)
					end do
					end do
				else if (ori .eq. 2)then
					l2 = m1
					m2 = l1
					allocate(ua2(bufferLength,9,l2,m2))
					do mm = 1,m2
					do ll = 1,l2
						ua2(:,:,ll,mm) = ua(:,:,mm,ll)
					end do
					end do
				else if (ori .eq. 3)then
					l2 = l1
					m2 = m1
					allocate(ua2(bufferLength,9,l2,m2))
					do mm = 1,m2
					do ll = 1,l2
						ua2(:,:,ll,mm) = ua(:,:,l1-ll+1,mm)
					end do
					end do
				else if (ori .eq. 4)then
					l2 = m1
					m2 = l1
					allocate(ua2(bufferLength,9,l2,m2))
					do mm = 1,m2
					do ll = 1,l2
						ua2(:,:,ll,mm) = ua(:,:,l1-mm+1,m1-ll+1)
					end do
					end do
				end if
				!!*
				
				face = blk(m0)%bc(ksub)%face
				!!i- face
				if (face .eq. 1) then !!boundary  i- 
					do k=1,m2
					do j=1,l2
						i0 = is_bc
						j0 = js_bc+k-1
						k0 = ks_bc+j-1	
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0-nBuffer,j0,k0) = ua2(nBuffer,1,j,k)
							blk(m0)%q(2,i0-nBuffer,j0,k0) = ua2(nBuffer,2,j,k)
							blk(m0)%q(3,i0-nBuffer,j0,k0) = ua2(nBuffer,3,j,k)
							blk(m0)%q(4,i0-nBuffer,j0,k0) = ua2(nBuffer,4,j,k)
							blk(m0)%q(5,i0-nBuffer,j0,k0) = ua2(nBuffer,5,j,k)

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0-nBuffer,j0,k0) = ua2(nBuffer,6,j,k)
							end if
						end do 
					end do
					end do
				end if
				!!*
				
				!!i+ face
				if (face .eq. 2) then !!boundary  i+
					do k=1,m2
					do j=1,l2
						i0 = ie_bc
						j0 = js_bc+j-1
						k0 = ks_bc+k-1
							
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0+nBuffer,j0,k0) = ua2(nBuffer,1,j,k)
							blk(m0)%q(2,i0+nBuffer,j0,k0) = ua2(nBuffer,2,j,k)
							blk(m0)%q(3,i0+nBuffer,j0,k0) = ua2(nBuffer,3,j,k)
							blk(m0)%q(4,i0+nBuffer,j0,k0) = ua2(nBuffer,4,j,k)
							blk(m0)%q(5,i0+nBuffer,j0,k0) = ua2(nBuffer,5,j,k) 

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0+nBuffer,j0,k0) = ua2(nBuffer,6,j,k)
							end if
						end do 
					end do
					end do
				end if
				!!*
				
				!!j- face
				if (face .eq. 3) then	!!boundary  j-
					do k=1,m2
					do i=1,l2
						i0 = is_bc+i-1
						j0 = js_bc
						k0 = ks_bc+k-1
							
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0,j0-nBuffer,k0) = ua2(nBuffer,1,i,k)
							blk(m0)%q(2,i0,j0-nBuffer,k0) = ua2(nBuffer,2,i,k)
							blk(m0)%q(3,i0,j0-nBuffer,k0) = ua2(nBuffer,3,i,k)
							blk(m0)%q(4,i0,j0-nBuffer,k0) = ua2(nBuffer,4,i,k)
							blk(m0)%q(5,i0,j0-nBuffer,k0) = ua2(nBuffer,5,i,k)

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0,j0-nBuffer,k0) = ua2(nBuffer,6,i,k)
							end if
						end do 
					end do
					end do
				end if
				!!*
				
				!!j+ face
				if (face .eq. 4) then	!!boundary  j+
					do k=1,m2
					do i=1,l2
						i0 = is_bc+k-1
						j0 = je_bc
						k0 = ks_bc+i-1
							
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0,j0+nBuffer,k0) = ua2(nBuffer,1,i,k) 
							blk(m0)%q(2,i0,j0+nBuffer,k0) = ua2(nBuffer,2,i,k)
							blk(m0)%q(3,i0,j0+nBuffer,k0) = ua2(nBuffer,3,i,k)
							blk(m0)%q(4,i0,j0+nBuffer,k0) = ua2(nBuffer,4,i,k)
							blk(m0)%q(5,i0,j0+nBuffer,k0) = ua2(nBuffer,5,i,k)

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0,j0+nBuffer,k0) = ua2(nBuffer,6,i,k)
							end if
						end do
					end do
					end do
				end if
				!!*

				!!k- face
				if (face .eq. 5) then	!!boundary  j+
					do j=1,m2
					do i=1,l2
						i0 = is_bc+j-1
						j0 = js_bc+i-1
						k0 = ks_bc
							
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0,j0,k0-nBuffer) = ua2(nBuffer,1,i,j) 
							blk(m0)%q(2,i0,j0,k0-nBuffer) = ua2(nBuffer,2,i,j)
							blk(m0)%q(3,i0,j0,k0-nBuffer) = ua2(nBuffer,3,i,j)
							blk(m0)%q(4,i0,j0,k0-nBuffer) = ua2(nBuffer,4,i,j)
							blk(m0)%q(5,i0,j0,k0-nBuffer) = ua2(nBuffer,5,i,j)

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0,j0,k0-nBuffer) = ua2(nBuffer,6,i,j)
							end if
						end do
					end do
					end do
				end if
				!!*

				!!k+ face
				if (face .eq. 6) then	!!boundary  j+
					do j=1,m2
					do i=1,l2
						i0 = is_bc+i-1
						j0 = js_bc+j-1
						k0 = ke_bc
								
						do nBuffer = 1, bufferLength
							blk(m0)%q(1,i0,j0,k0+nBuffer) = ua2(nBuffer,1,i,j) 
							blk(m0)%q(2,i0,j0,k0+nBuffer) = ua2(nBuffer,2,i,j)
							blk(m0)%q(3,i0,j0,k0+nBuffer) = ua2(nBuffer,3,i,j)
							blk(m0)%q(4,i0,j0,k0+nBuffer) = ua2(nBuffer,4,i,j)
							blk(m0)%q(5,i0,j0,k0+nBuffer) = ua2(nBuffer,5,i,j)

							if (iflag_turbulence .eq. iflag_sa) then
								blk(m0)%nut(i0,j0,k0+nBuffer) = ua2(nBuffer,6,i,j)
							end if
						end do
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
