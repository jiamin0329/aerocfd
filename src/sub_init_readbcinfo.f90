!!=========================================================
!!	read bcinfo
!!
!!	purpose: 1. read block interface condition info
!!           2. read physical boundary condition info
!!           3. allocate buffer(send/recv) block 
!!              for block interface
!!           4. bc info is load into both blk and blk0
!!                                                               
!!	author: jiamin xu
!!	date:   2012.11.22
!!=========================================================
subroutine init_readbcinfo
	use glbindex_var
	use blk_var
	use index_var
	use mpi_var
	use output_var
	use flag_var
	implicit none
	!!***************************************************************
	integer                          :: num_subface
	integer,dimension(:),allocatable :: f_no
	integer,dimension(:),allocatable :: face
	integer,dimension(:),allocatable :: is0,ie0,js0,je0,ks0,ke0
	integer,dimension(:),allocatable :: blk_t
	integer,dimension(:),allocatable :: subface_t
	integer,dimension(:),allocatable :: rank_c,rank_t
	integer,dimension(:),allocatable :: ori
	integer,dimension(:),allocatable :: fwh
	!!***************************************************************
	integer                          :: m_t,ksub_t
	integer                          :: is_s,ie_s,js_s,je_s,ks_s,ke_s
	integer                          :: is_r,ie_r,js_r,je_r,ks_r,ke_r
	!!***************************************************************	
	!! subface info
	integer :: idim_subface, jdim_subface, kdim_subface
	real*8,dimension(:,:,:),allocatable :: temp_subface
	integer :: i,j,k
	!!reading bcinfo file
	if (myid .eq. root) then
		write(*,*)'***********************************************************'
		write(*,*)"reading bcinfo file:"
		write(*,*)'***********************************************************'
		open(99,file = bcinfo)
		read(99,*)
		read(99,*)
		read(99,*)
	end if
	!!*
	
	do m = 1,nblock
		!!get number of subfaces in the block
		if(myid .eq. root) then
			!!reading bc info
			read(99,*) 
			read(99,*) 
			read(99,*) num_subface
			read(99,*)
		end if
		call MPI_BCAST(num_subface, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		
		!!global index
		blk0(m)%num_subface = num_subface
		allocate(blk0(m)%bc0(num_subface))

		allocate(f_no          (num_subface))
		allocate(face          (num_subface))
		allocate(is0           (num_subface))
		allocate(ie0           (num_subface))
		allocate(js0           (num_subface))
		allocate(je0           (num_subface))
		allocate(ks0           (num_subface))
		allocate(ke0           (num_subface))
		allocate(blk_t         (num_subface))
		allocate(subface_t     (num_subface))
		allocate(ori           (num_subface))
		allocate(fwh           (num_subface))
		allocate(rank_c        (num_subface))
		allocate(rank_t        (num_subface))

		!! Reading bcinfo from inputfile, from root thread
		if(myid .eq. root) then
			do ksub = 1,num_subface
				if      (iflag_dimension .eq. iflag_2d)then
					read(99,*) f_no(ksub),face(ksub),                  &
					           is0(ksub),ie0(ksub),js0(ksub),je0(ksub),&
					           blk_t(ksub),subface_t(ksub),ori(ksub),  &
					           fwh(ksub)
					
					ks0(ksub) = 1
					ke0(ksub) = 1           
					call m2node(rank_c(ksub),          m,myn,lastn,nblock,numprocs)
					call m2node(rank_t(ksub),blk_t(ksub),myn,lastn,nblock,numprocs)  
				else if (iflag_dimension .eq. iflag_3d)then
					read(99,*) f_no(ksub),face(ksub),                                        &
					           is0(ksub),ie0(ksub),js0(ksub),je0(ksub),ks0(ksub),ke0(ksub),  &
					           blk_t(ksub),subface_t(ksub),ori(ksub),fwh(ksub)
					           
					call m2node(rank_c(ksub),          m,myn,lastn,nblock,numprocs)
					call m2node(rank_t(ksub),blk_t(ksub),myn,lastn,nblock,numprocs)  
				end if     
			end do
			read(99,*)
		end if
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!!*

		call MPI_BCAST(f_no,      num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(face,      num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(is0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(ie0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(js0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(je0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(ks0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(ke0,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk_t,     num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(subface_t, num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(ori,       num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(rank_c,    num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(rank_t,    num_subface,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		!!id of thread to receive
		call m2node(dest_id,m,myn,lastn,nblock,numprocs)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!!*
		
		if (myid .eq. 0)then
			if (dest_id .eq. 0)then
				m0 = m
!!				write(*,25) m0,blk_loop,myid,num_subface
!!25			format("(native)block",I5 ," of",I5 ," in rank",I4 ," total subface =",I4)
				allocate(blk(m0)%bc(num_subface))
				blk(m0)%num_subface = num_subface
				do ksub = 1,num_subface
					blk(m0)%bc(ksub)%f_no      = f_no          (ksub)
					blk(m0)%bc(ksub)%face      = face          (ksub)
					blk(m0)%bc(ksub)%is        = is0           (ksub)
					blk(m0)%bc(ksub)%ie        = ie0           (ksub)
					blk(m0)%bc(ksub)%js        = js0           (ksub)
					blk(m0)%bc(ksub)%je        = je0           (ksub) 
					blk(m0)%bc(ksub)%ks        = ks0           (ksub)
					blk(m0)%bc(ksub)%ke        = ke0           (ksub)
					blk(m0)%bc(ksub)%dim_f     =(ke0(ksub) - ks0(ksub) + 1) &
					                           *(je0(ksub) - js0(ksub) + 1) &
					                           *(ie0(ksub) - is0(ksub) + 1)
					blk(m0)%bc(ksub)%blk_t     = blk_t         (ksub)
					blk(m0)%bc(ksub)%subface_t = subface_t     (ksub)
					blk(m0)%bc(ksub)%ori       = ori           (ksub)
					blk(m0)%bc(ksub)%fwh       = fwh           (ksub)
					blk(m0)%bc(ksub)%rank_c    = rank_c        (ksub)
					blk(m0)%bc(ksub)%rank_t    = rank_t        (ksub)
				end do
			else if(dest_id .gt. 0)then
				call MPI_send(f_no      ,num_subface, MPI_INTEGER,dest_id,89,MPI_COMM_WORLD,ierr)
				call MPI_send(face      ,num_subface, MPI_INTEGER,dest_id,90,MPI_COMM_WORLD,ierr)
				call MPI_send(is0       ,num_subface, MPI_INTEGER,dest_id,91,MPI_COMM_WORLD,ierr)
				call MPI_send(ie0       ,num_subface, MPI_INTEGER,dest_id,92,MPI_COMM_WORLD,ierr)
				call MPI_send(js0       ,num_subface, MPI_INTEGER,dest_id,93,MPI_COMM_WORLD,ierr)
				call MPI_send(je0       ,num_subface, MPI_INTEGER,dest_id,94,MPI_COMM_WORLD,ierr)
				call MPI_send(ks0       ,num_subface, MPI_INTEGER,dest_id,95,MPI_COMM_WORLD,ierr)
				call MPI_send(ke0       ,num_subface, MPI_INTEGER,dest_id,96,MPI_COMM_WORLD,ierr)
				call MPI_send(blk_t     ,num_subface, MPI_INTEGER,dest_id,97,MPI_COMM_WORLD,ierr)
				call MPI_send(subface_t ,num_subface, MPI_INTEGER,dest_id,98,MPI_COMM_WORLD,ierr)
				call MPI_send(ori       ,num_subface, MPI_INTEGER,dest_id,99,MPI_COMM_WORLD,ierr)
				call MPI_send(fwh       ,num_subface, MPI_INTEGER,dest_id,86,MPI_COMM_WORLD,ierr)
				call MPI_send(rank_c    ,num_subface, MPI_INTEGER,dest_id,88,MPI_COMM_WORLD,ierr)
				call MPI_send(rank_t    ,num_subface, MPI_INTEGER,dest_id,87,MPI_COMM_WORLD,ierr)     		
			end if
		end if

		source_id = 0
		
		if(myid .gt. 0 .and. myid .eq. dest_id)then
			m0 = m - myid*myn
!!			write(*,30) m0,blk_loop,myid,num_subface
!!30		format("(native)block",I5 ," of",I5 ," in rank",I4 ," total subface =",I4)
			blk(m0)%num_subface = num_subface
			allocate(blk(m0)%bc(num_subface))
			call MPI_recv(rank_t    ,num_subface,MPI_INTEGER,source_id,87,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(rank_c    ,num_subface,MPI_INTEGER,source_id,88,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(f_no      ,num_subface,MPI_INTEGER,source_id,89,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(face      ,num_subface,MPI_INTEGER,source_id,90,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(is0       ,num_subface,MPI_INTEGER,source_id,91,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(ie0       ,num_subface,MPI_INTEGER,source_id,92,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(js0       ,num_subface,MPI_INTEGER,source_id,93,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(je0       ,num_subface,MPI_INTEGER,source_id,94,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(ks0       ,num_subface,MPI_INTEGER,source_id,95,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(ke0       ,num_subface,MPI_INTEGER,source_id,96,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(blk_t     ,num_subface,MPI_INTEGER,source_id,97,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(subface_t ,num_subface,MPI_INTEGER,source_id,98,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(ori       ,num_subface,MPI_INTEGER,source_id,99,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(fwh       ,num_subface,MPI_INTEGER,source_id,86,MPI_COMM_WORLD,status,ierr)
			
			do ksub = 1,num_subface
				blk(m0)%bc(ksub)%f_no      = f_no (ksub)
				blk(m0)%bc(ksub)%face      = face (ksub)
				blk(m0)%bc(ksub)%is        = is0  (ksub)
				blk(m0)%bc(ksub)%ie        = ie0  (ksub)
				blk(m0)%bc(ksub)%js        = js0  (ksub)
				blk(m0)%bc(ksub)%je        = je0  (ksub) 
				blk(m0)%bc(ksub)%ks        = ks0  (ksub)
				blk(m0)%bc(ksub)%ke        = ke0  (ksub)
				blk(m0)%bc(ksub)%dim_f     =(ke0 (ksub) - ks0 (ksub) + 1) &
				                           *(je0 (ksub) - js0 (ksub) + 1) &
				                           *(ie0 (ksub) - is0 (ksub) + 1)
				blk(m0)%bc(ksub)%blk_t     = blk_t (ksub)
				blk(m0)%bc(ksub)%subface_t = subface_t(ksub)
				blk(m0)%bc(ksub)%ori       = ori (ksub)
				blk(m0)%bc(ksub)%fwh       = fwh (ksub)
				blk(m0)%bc(ksub)%rank_c    = rank_c(ksub)
				blk(m0)%bc(ksub)%rank_t    = rank_t(ksub)
			end do
			!!*
		end if
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		do ksub = 1,num_subface
			blk0(m)%bc0(ksub)%f_no     = f_no(ksub)
			blk0(m)%bc0(ksub)%face     = face(ksub)
			blk0(m)%bc0(ksub)%is       = is0 (ksub)
			blk0(m)%bc0(ksub)%ie       = ie0 (ksub)
			blk0(m)%bc0(ksub)%js       = js0 (ksub)
			blk0(m)%bc0(ksub)%je       = je0 (ksub)
			blk0(m)%bc0(ksub)%ks       = ks0 (ksub)
			blk0(m)%bc0(ksub)%ke       = ke0 (ksub)
			blk0(m)%bc0(ksub)%dim_f    =(ke0 (ksub) - ks0 (ksub) + 1) &
				                       *(je0 (ksub) - js0 (ksub) + 1) &
				                       *(ie0 (ksub) - is0 (ksub) + 1)
			blk0(m)%bc0(ksub)%blk_t    = blk_t(ksub)
			blk0(m)%bc0(ksub)%subface_t= subface_t(ksub)
			blk0(m)%bc0(ksub)%ori      = ori      (ksub)
			blk0(m)%bc0(ksub)%fwh      = fwh      (ksub)
			blk0(m)%bc0(ksub)%rank_c   = rank_c   (ksub)
			blk0(m)%bc0(ksub)%rank_t   = rank_t   (ksub)
		end do

		deallocate(f_no     )
		deallocate(face     )
		deallocate(is0      )
		deallocate(ie0      )
		deallocate(js0      )
		deallocate(je0      )
		deallocate(ks0      )
		deallocate(ke0      )
		deallocate(blk_t    )
		deallocate(subface_t)
		deallocate(ori      )
		deallocate(fwh      )
		deallocate(rank_c   )
		deallocate(rank_t   )
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	end do

	if(myid .eq. 0)then  
		write(*,*)'***********************************************************'
		write(*,*)'finish reading bcinfo!'
		write(*,*)'***********************************************************'
		close(99)
!		pause
	end if 
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
	return
end subroutine


subroutine init_readbcinfo2
	use glbindex_var
	use blk_var
	use index_var
	use mpi_var
	use output_var
	use flag_var
	implicit none
	!! subface info
	integer :: idim_subface, jdim_subface, kdim_subface
	real*8,dimension(:,:,:),allocatable :: temp_subface
	integer :: i,j,k

	integer :: m_t,ksub_t
	integer :: is_s,ie_s,js_s,je_s,ks_s,ke_s
	integer :: is_r,ie_r,js_r,je_r,ks_r,ke_r

	do m = 1,nblock
		do ksub = 1,blk0(m)%num_subface
			idim_subface = blk0(m)%bc0(ksub)%ie - blk0(m)%bc0(ksub)%is + 1
			jdim_subface = blk0(m)%bc0(ksub)%je - blk0(m)%bc0(ksub)%js + 1
			kdim_subface = blk0(m)%bc0(ksub)%ke - blk0(m)%bc0(ksub)%ks + 1	
			
			print *, "********", myid, m, ksub, "face:", blk0(m)%bc0(ksub)%face, idim_subface, jdim_subface, kdim_subface
			if      (blk0(m)%bc0(ksub)%face .eq. 1 .or. blk0(m)%bc0(ksub)%face .eq. 2) then
				allocate(temp_subface(3, jdim_subface, kdim_subface))
				allocate(blk0(m)%bc0(ksub)%x(jdim_subface,kdim_subface))
				allocate(blk0(m)%bc0(ksub)%y(jdim_subface,kdim_subface))
				allocate(blk0(m)%bc0(ksub)%z(jdim_subface,kdim_subface))
			else if (blk0(m)%bc0(ksub)%face .eq. 3 .or. blk0(m)%bc0(ksub)%face .eq. 4) then
				allocate(temp_subface(3, idim_subface, kdim_subface))
				allocate(blk0(m)%bc0(ksub)%x(idim_subface,kdim_subface))
				allocate(blk0(m)%bc0(ksub)%y(idim_subface,kdim_subface))
				allocate(blk0(m)%bc0(ksub)%z(idim_subface,kdim_subface))
			else if	(blk0(m)%bc0(ksub)%face .eq. 5 .or. blk0(m)%bc0(ksub)%face .eq. 6) then
				allocate(temp_subface(3, idim_subface, jdim_subface))
				allocate(blk0(m)%bc0(ksub)%x(idim_subface,jdim_subface))
				allocate(blk0(m)%bc0(ksub)%y(idim_subface,jdim_subface))
				allocate(blk0(m)%bc0(ksub)%z(idim_subface,jdim_subface))
			end if 
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

			call m2node(dest_id,m,myn,lastn,nblock,numprocs)
			if      (dest_id .eq. myid) then
				m0 = m - myid*myn
				!!print *, "********", myid, m, blk0(m)%ni, blk0(m)%nj, blk0(m)%nk
				!! i- face
				if (blk0(m)%bc0(ksub)%face .eq. 1) then
				do k = blk0(m)%bc0(ksub)%ks, blk0(m)%bc0(ksub)%ke
				do j = blk0(m)%bc0(ksub)%js, blk0(m)%bc0(ksub)%je
					temp_subface(1,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%x(blk0(m)%is0,j,k)
					temp_subface(2,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%y(blk0(m)%is0,j,k)
					temp_subface(3,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%z(blk0(m)%is0,j,k)
				end do
				end do
				end if

				!! i+ face
				if (blk0(m)%bc0(ksub)%face .eq. 2) then
				do k = blk0(m)%bc0(ksub)%ks, blk0(m)%bc0(ksub)%ke
				do j = blk0(m)%bc0(ksub)%js, blk0(m)%bc0(ksub)%je
					temp_subface(1,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%x(blk0(m)%ie0,j,k)
					temp_subface(2,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%y(blk0(m)%ie0,j,k)
					temp_subface(3,j-blk0(m)%bc0(ksub)%js+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%z(blk0(m)%ie0,j,k)
				end do
				end do
				end if

				!! j- face
				if (blk0(m)%bc0(ksub)%face .eq. 3) then
				do k = blk0(m)%bc0(ksub)%ks, blk0(m)%bc0(ksub)%ke
				do i = blk0(m)%bc0(ksub)%is, blk0(m)%bc0(ksub)%ie
					temp_subface(1,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%x(i,blk0(m)%js0,k)
					temp_subface(2,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%y(i,blk0(m)%js0,k)
					temp_subface(3,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%z(i,blk0(m)%js0,k)
				end do
				end do
				end if

				!! j+ face
				if (blk0(m)%bc0(ksub)%face .eq. 4) then
				do k = blk0(m)%bc0(ksub)%ks, blk0(m)%bc0(ksub)%ke
				do i = blk0(m)%bc0(ksub)%is, blk0(m)%bc0(ksub)%ie
					temp_subface(1,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%x(i,blk0(m)%je0,k)
					temp_subface(2,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%y(i,blk0(m)%je0,k)
					temp_subface(3,i-blk0(m)%bc0(ksub)%is+1,k-blk0(m)%bc0(ksub)%ks+1) = blk(m-myid*myn)%z(i,blk0(m)%je0,k)
				end do
				end do
				end if

				!! k- face
				if (blk0(m)%bc0(ksub)%face .eq. 5) then
				do j = blk0(m)%bc0(ksub)%js, blk0(m)%bc0(ksub)%je
				do i = blk0(m)%bc0(ksub)%is, blk0(m)%bc0(ksub)%ie
					temp_subface(1,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%x(i,j,blk0(m)%ks0)
					temp_subface(2,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%y(i,j,blk0(m)%ks0)
					temp_subface(3,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%z(i,j,blk0(m)%ks0)
				end do
				end do
				end if

				!! k+ face
				if (blk0(m)%bc0(ksub)%face .eq. 6) then
				do j = blk0(m)%bc0(ksub)%js, blk0(m)%bc0(ksub)%je
				do i = blk0(m)%bc0(ksub)%is, blk0(m)%bc0(ksub)%ie
					temp_subface(1,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%x(i,j,blk0(m)%ke0)
					temp_subface(2,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%y(i,j,blk0(m)%ke0)
					temp_subface(3,i-blk0(m)%bc0(ksub)%is+1,j-blk0(m)%bc0(ksub)%js+1) = blk(m-myid*myn)%z(i,j,blk0(m)%ke0)
				end do
				end do
				end if
			end if
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

			call m2node(dest_id,m,myn,lastn,nblock,numprocs)
			call MPI_BCAST(temp_subface, idim_subface*jdim_subface*kdim_subface, MPI_REAL8, dest_id, MPI_COMM_WORLD,ierr)
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
			do k = blk0(m)%bc0(ksub)%ks-blk0(m)%bc0(ksub)%ks+1, blk0(m)%bc0(ksub)%ke-blk0(m)%bc0(ksub)%ks+1
			do j = blk0(m)%bc0(ksub)%js-blk0(m)%bc0(ksub)%js+1, blk0(m)%bc0(ksub)%je-blk0(m)%bc0(ksub)%js+1
			do i = blk0(m)%bc0(ksub)%is-blk0(m)%bc0(ksub)%is+1, blk0(m)%bc0(ksub)%ie-blk0(m)%bc0(ksub)%is+1
				if 		(blk0(m)%bc0(ksub)%face .eq. 1) then
					blk0(m)%bc0(ksub)%x(j,k) = temp_subface(1,j,k)
					blk0(m)%bc0(ksub)%y(j,k) = temp_subface(2,j,k)
					blk0(m)%bc0(ksub)%z(j,k) = temp_subface(3,j,k)
				else if (blk0(m)%bc0(ksub)%face .eq. 2) then
					blk0(m)%bc0(ksub)%x(j,k) = temp_subface(1,j,k)
					blk0(m)%bc0(ksub)%y(j,k) = temp_subface(2,j,k)
					blk0(m)%bc0(ksub)%z(j,k) = temp_subface(3,j,k)
				else if (blk0(m)%bc0(ksub)%face .eq. 3) then
					blk0(m)%bc0(ksub)%x(i,k) = temp_subface(1,i,k)
					blk0(m)%bc0(ksub)%y(i,k) = temp_subface(2,i,k)
					blk0(m)%bc0(ksub)%z(i,k) = temp_subface(3,i,k)
				else if	(blk0(m)%bc0(ksub)%face .eq. 4) then
					blk0(m)%bc0(ksub)%x(i,k) = temp_subface(1,i,k)
					blk0(m)%bc0(ksub)%y(i,k) = temp_subface(2,i,k)
					blk0(m)%bc0(ksub)%z(i,k) = temp_subface(3,i,k)
				else if	(blk0(m)%bc0(ksub)%face .eq. 5) then
					blk0(m)%bc0(ksub)%x(i,j) = temp_subface(1,i,j)
					blk0(m)%bc0(ksub)%y(i,j) = temp_subface(2,i,j)
					blk0(m)%bc0(ksub)%z(i,j) = temp_subface(3,i,j)
				else if	(blk0(m)%bc0(ksub)%face .eq. 6) then
					blk0(m)%bc0(ksub)%x(i,j) = temp_subface(1,i,j)
					blk0(m)%bc0(ksub)%y(i,j) = temp_subface(2,i,j)
					blk0(m)%bc0(ksub)%z(i,j) = temp_subface(3,i,j)
				end if
			end do
			end do
			end do

			deallocate(temp_subface)
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		end do
		!!*		
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	end do !! end main loop m, global block loop
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*

	do m0 = 1,blk_loop
		do ksub = 1,blk(m0)%num_subface
			if(blk(m0)%bc(ksub)%blk_t .gt. 0)then
				is_s = blk(m0)%bc(ksub)%is
				ie_s = blk(m0)%bc(ksub)%ie
				js_s = blk(m0)%bc(ksub)%js
				je_s = blk(m0)%bc(ksub)%je
				ks_s = blk(m0)%bc(ksub)%ks
				ke_s = blk(m0)%bc(ksub)%ke
			
				allocate(blk(m0)%bc(ksub)%buffer_send    (bufferLength,9,is_s:ie_s,js_s:je_s,ks_s:ke_s))
				!!allocate(blk(m0)%bc(ksub)%buffer_send_tur(bufferLength,  is_s:ie_s,js_s:je_s,ks_s:ke_s))
			
				!! for overlap method average procedure
				allocate(blk(m0)%bc(ksub)%buffer(9,is_s:ie_s,js_s:je_s,ks_s:ke_s))

				m_t    = blk(m0)%bc(ksub)%blk_t
				ksub_t = blk(m0)%bc(ksub)%subface_t
			
				is_r = blk0(m_t)%bc0(ksub_t)%is
				ie_r = blk0(m_t)%bc0(ksub_t)%ie
				js_r = blk0(m_t)%bc0(ksub_t)%js
				je_r = blk0(m_t)%bc0(ksub_t)%je
				ks_r = blk0(m_t)%bc0(ksub_t)%ks
				ke_r = blk0(m_t)%bc0(ksub_t)%ke
			
				allocate(blk(m0)%bc(ksub)%buffer_recv    (bufferLength,9,is_r:ie_r,js_r:je_r,ks_r:ke_r))
				!!allocate(blk(m0)%bc(ksub)%buffer_recv_tur(bufferLength,  is_r:ie_r,js_r:je_r,ks_r:ke_r))
			
				blk(m0)%bc(ksub)%buffer_send     = 0.d0
				blk(m0)%bc(ksub)%buffer_recv     = 0.d0
				!!blk(m0)%bc(ksub)%buffer_send_tur = 0.d0
				!!blk(m0)%bc(ksub)%buffer_recv_tur = 0.d0
			end if
		end do
	end do
	!!*

	return
end subroutine
