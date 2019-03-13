!!=========================================================================== 
!! memory allocation
!!
!! purpose: 1. read total block number and comupte myn&lastn(blk_loop)
!!          2. read block dimension and allocate memory in each processor
!!          3. root rank is used to read grid file and then x,y,z are sent 
!!             to the proper processor
!!                                                                     
!! author: jiamin xu                                    
!! date:   2012.11.22                      
!!============================================================================ 
subroutine init_allocatememory
	use ns_const
	use blk_var
	use mpi_var
	use index_var
	use glbindex_var
	use output_var
	use flag_var
	implicit none
	
	integer :: i,j,k
	integer :: temp = 0
	integer :: idim,jdim,kdim
	integer,allocatable,dimension(:)     :: tempni,tempnj,tempnk  !!temp array for block dimension
	real*8, allocatable,dimension(:,:,:) :: tempx, tempy, tempz   !!temp array for coordinate
	integer,allocatable,dimension(:)     :: global_is,  global_ie,  global_js,  global_je,  global_ks,  global_ke
	integer,allocatable,dimension(:)     :: global_is0, global_ie0, global_js0, global_je0, global_ks0, global_ke0
	integer,allocatable,dimension(:)     :: global_is1, global_ie1, global_js1, global_je1, global_ks1, global_ke1
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	integer :: subface

	if(myid .eq. root)then
		open(99,file = grid)
		read(99,*) nblock
	end if

	!!read block dimension
	allocate(tempni(nblock))
	allocate(tempnj(nblock))
	allocate(tempnk(nblock))
	allocate(global_is (nblock))
	allocate(global_ie (nblock))
	allocate(global_js (nblock))
	allocate(global_je (nblock))
	allocate(global_ks (nblock))
	allocate(global_ke (nblock))
	allocate(global_is0(nblock))
	allocate(global_ie0(nblock))
	allocate(global_js0(nblock))
	allocate(global_je0(nblock))
	allocate(global_ks0(nblock))
	allocate(global_ke0(nblock))
	allocate(global_is1(nblock))
	allocate(global_ie1(nblock))
	allocate(global_js1(nblock))
	allocate(global_je1(nblock))
	allocate(global_ks1(nblock))
	allocate(global_ke1(nblock))
	!!*
	
	if (myid .eq. root)then
		do m = 1,nblock
			if     (iflag_dimension .eq. iflag_2d) then
				read(99,*) tempni(m),tempnj(m)
				tempnk(m) = 1
			else if(iflag_dimension .eq. iflag_3d) then
				read(99,*) tempni(m),tempnj(m),tempnk(m)
			end if
			
			global_is0(m) = 1
			global_ie0(m) = tempni(m) 
			global_js0(m) = 1
			global_je0(m) = tempnj(m) 
			global_ks0(m) = 1
			global_ke0(m) = tempnk(m)

			global_is(m) = global_is0(m)
			global_ie(m) = global_ie0(m)
			global_js(m) = global_js0(m)
			global_je(m) = global_je0(m)
			global_ks(m) = global_ks0(m)
			global_ke(m) = global_ke0(m)

			do ksub = 1, blk0(m)%num_subface
				if (blk0(m)%bc0(ksub)%blk_t .gt. 0) then 
					if (blk0(m)%bc0(ksub)%face .eq. 1 ) then
						global_is(m) = global_is0(m) - bufferLength
					end if
					if (blk0(m)%bc0(ksub)%face .eq. 2 ) then
						global_ie(m) = global_ie0(m) + bufferLength
					end if
					if (blk0(m)%bc0(ksub)%face .eq. 3 ) then
						global_js(m) = global_js0(m) - bufferLength
					end if
					if (blk0(m)%bc0(ksub)%face .eq. 4 ) then
						global_je(m) = global_je0(m) + bufferLength
					end if
					if (blk0(m)%bc0(ksub)%face .eq. 5 ) then
						global_ks(m) = global_ks0(m) - bufferLength
					end if
					if (blk0(m)%bc0(ksub)%face .eq. 6 ) then
						global_ke(m) = global_ke0(m) + bufferLength
					end if
				end if
			end do

			if (global_is0(m) .eq. global_is(m)) then
				global_is1(m) = global_is0(m) + 1
			else 
				global_is1(m) = global_is0(m)
			end if

			if (global_ie0(m) .eq. global_ie(m)) then
				global_ie1(m) = global_ie0(m) - 1
			else 
				global_ie1(m) = global_ie0(m)
			end if

			if (global_js0(m) .eq. global_js(m)) then
				global_js1(m) = global_js0(m) + 1
			else 
				global_js1(m) = global_js0(m)
			end if

			if (global_je0(m) .eq. global_je(m)) then
				global_je1(m) = global_je0(m) - 1
			else 
				global_je1(m) = global_je0(m)
			end if

			if (global_ks0(m) .eq. global_ks(m) .and. iflag_dimension .eq. iflag_3d) then
				global_ks1(m) = global_ks0(m) + 1
			else 
				global_ks1(m) = global_ks0(m)
			end if

			if (global_ke0(m) .eq. global_ke(m) .and. iflag_dimension .eq. iflag_3d) then
				global_ke1(m) = global_ke0(m) - 1
			else 
				global_ke1(m) = global_ke0(m)
			end if

			blk0(m)%ni = global_ie0(m)-global_is0(m)+1
			blk0(m)%nj = global_je0(m)-global_js0(m)+1
			blk0(m)%nk = global_ke0(m)-global_ks0(m)+1

			blk0(m)%is = global_is(m)
			blk0(m)%ie = global_ie(m)
			blk0(m)%js = global_js(m)
			blk0(m)%je = global_je(m)
			blk0(m)%ks = global_ks(m)
			blk0(m)%ke = global_ke(m)

			blk0(m)%is0 = global_is0(m)
			blk0(m)%ie0 = global_ie0(m)
			blk0(m)%js0 = global_js0(m)
			blk0(m)%je0 = global_je0(m)
			blk0(m)%ks0 = global_ks0(m)
			blk0(m)%ke0 = global_ke0(m)

			blk0(m)%is1 = global_is1(m)
			blk0(m)%ie1 = global_ie1(m)
			blk0(m)%js1 = global_js1(m)
			blk0(m)%je1 = global_je1(m)
			blk0(m)%ks1 = global_ks1(m)
			blk0(m)%ke1 = global_ke1(m)
		end do
	end if
	
	call MPI_BCAST(tempni,    nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tempnj,    nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tempnk,    nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_is, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ie, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_js, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_je, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ks, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ke, nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_is0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ie0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_js0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_je0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ks0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ke0,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_is1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ie1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_js1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_je1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ks1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(global_ke1,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	!!*
	
	!!tempni/tempnj/tempnk ---> nx,ny,nz
	if      (myid .lt. numprocs-1) then
		do m0 = 1,myn
			blk(m0)%ni = global_ie(m0+myid*myn)-global_is(m0+myid*myn)+1
			blk(m0)%nj = global_je(m0+myid*myn)-global_js(m0+myid*myn)+1
			blk(m0)%nk = global_ke(m0+myid*myn)-global_ks(m0+myid*myn)+1
			temp = temp + blk(m0)%ni*blk(m0)%nj*blk(m0)%nk

			blk(m0)%is  = global_is (m0+myid*myn)
			blk(m0)%ie  = global_ie (m0+myid*myn)
			blk(m0)%js  = global_js (m0+myid*myn)
			blk(m0)%je  = global_je (m0+myid*myn)
			blk(m0)%ks  = global_ks (m0+myid*myn)
			blk(m0)%ke  = global_ke (m0+myid*myn)

			blk(m0)%is0 = global_is0(m0+myid*myn)
			blk(m0)%ie0 = global_ie0(m0+myid*myn)
			blk(m0)%js0 = global_js0(m0+myid*myn)
			blk(m0)%je0 = global_je0(m0+myid*myn)
			blk(m0)%ks0 = global_ks0(m0+myid*myn)
			blk(m0)%ke0 = global_ke0(m0+myid*myn)

			blk(m0)%is1 = global_is1(m0+myid*myn)
			blk(m0)%ie1 = global_ie1(m0+myid*myn)
			blk(m0)%js1 = global_js1(m0+myid*myn)
			blk(m0)%je1 = global_je1(m0+myid*myn)
			blk(m0)%ks1 = global_ks1(m0+myid*myn)
			blk(m0)%ke1 = global_ke1(m0+myid*myn)
		end do
	else if (myid .eq. numprocs-1) then
		do m0 = 1,lastn
			blk(m0)%ni = global_ie(nblock-lastn+m0)-global_is(nblock-lastn+m0)+1
			blk(m0)%nj = global_je(nblock-lastn+m0)-global_js(nblock-lastn+m0)+1
			blk(m0)%nk = global_ke(nblock-lastn+m0)-global_ks(nblock-lastn+m0)+1
			temp = temp + blk(m0)%ni*blk(m0)%nj*blk(m0)%nk

			blk(m0)%is  = global_is (nblock-lastn+m0)
			blk(m0)%ie  = global_ie (nblock-lastn+m0)
			blk(m0)%js  = global_js (nblock-lastn+m0)
			blk(m0)%je  = global_je (nblock-lastn+m0)
			blk(m0)%ks  = global_ks (nblock-lastn+m0)
			blk(m0)%ke  = global_ke (nblock-lastn+m0)

			blk(m0)%is0 = global_is0(nblock-lastn+m0)
			blk(m0)%ie0 = global_ie0(nblock-lastn+m0)
			blk(m0)%js0 = global_js0(nblock-lastn+m0)
			blk(m0)%je0 = global_je0(nblock-lastn+m0)
			blk(m0)%ks0 = global_ks0(nblock-lastn+m0)
			blk(m0)%ke0 = global_ke0(nblock-lastn+m0)

			blk(m0)%is1 = global_is1(nblock-lastn+m0)
			blk(m0)%ie1 = global_ie1(nblock-lastn+m0)
			blk(m0)%js1 = global_js1(nblock-lastn+m0)
			blk(m0)%je1 = global_je1(nblock-lastn+m0)
			blk(m0)%ks1 = global_ks1(nblock-lastn+m0)
			blk(m0)%ke1 = global_ke1(nblock-lastn+m0)
		end do
	end if
	
	write(*,5), temp,myid
5 format(I8," nodes in processor",I3)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	if(myid .eq. root)then
		totalnodes = 0
		do m = 1,nblock
			totalnodes = totalnodes + blk0(m)%ni*blk0(m)%nj*blk0(m)%nk
		end do
		write (*,*)"***************"  
		print *, 'total nodes = ',totalnodes
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*

	do m0 = 1,blk_loop
		!!block dimension
		is = blk(m0)%is; ie = blk(m0)%ie
		js = blk(m0)%js; je = blk(m0)%je
		ks = blk(m0)%ks; ke = blk(m0)%ke
		!!write(*,*)is,ie,js,je,ks,ke
		!!pause
		!!allocate memory for each block
		!!coordinates
		allocate (blk(m0)%x(is:ie,js:je,ks:ke))
		allocate (blk(m0)%y(is:ie,js:je,ks:ke))
		allocate (blk(m0)%z(is:ie,js:je,ks:ke))
		
		!!conservative variables
		allocate (blk(m0)%q      (5,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%q0(-1:0,5,is:ie,js:je,ks:ke))
		!!*
		
		!!inviscid term
		!!allocate (blk(m0)%e      (5,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%f      (5,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%g      (5,is:ie,js:je,ks:ke))
		!!inviscid flux
		!!allocate (blk(m0)%dedxi  (5,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%dfdeta (5,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%dgdzeta(5,is:ie,js:je,ks:ke))
		!!*
		
		!!primitive variables
		allocate (blk(m0)%pri_v(7,is:ie,js:je,ks:ke))
		!!*
		
		!!viscous flux
		!!velocity gradients       
		allocate (blk(m0)%dudx (9,is:ie,js:je,ks:ke))
		!!stress tensor
		allocate (blk(m0)%tao  (6,is:ie,js:je,ks:ke))
		
		allocate (blk(m0)%amu  (3,is:ie,js:je,ks:ke))
		allocate (blk(m0)%vor    (is:ie,js:je,ks:ke))
		
		!!right hand side terms
		allocate (blk(m0)%rhsv (5,is:ie,js:je,ks:ke))
		allocate (blk(m0)%rhsi (5,is:ie,js:je,ks:ke))
		allocate (blk(m0)%rhs  (5,is:ie,js:je,ks:ke))
		!!*
		
		!!non-physical time step
		allocate (blk(m0)%dt        (is:ie,js:je,ks:ke))
		!!*
		
		!!turbulence model
		allocate (blk(m0)%length    (is:ie,js:je,ks:ke))
		!!sa model
		allocate (blk(m0)%nut       (is:ie,js:je,ks:ke))
		allocate (blk(m0)%sa_rhs    (is:ie,js:je,ks:ke))
		!!kwsst model
		!!allocate (blk(m0)%k         (is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%omg       (is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%kw_rhs1   (is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%kw_rhs2   (is:ie,js:je,ks:ke))
		!!*end turbulence model
		
		!!jacobian parameters
		allocate (blk(m0)%dxidx   (9,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%dxidxh  (9,is:ie,js:je,ks:ke))
		!!allocate (blk(m0)%dxidxdxi(9,is:ie,js:je,ks:ke))
		allocate (blk(m0)%inv_j     (is:ie,js:je,ks:ke))
		
		!!allocate (blk(m0)%dxidx2   (3,is:ie,js:je,ks:ke))
		allocate (blk(m0)%aalpha   (3,is:ie,js:je,ks:ke))
		allocate (blk(m0)%bbeta    (3,is:ie,js:je,ks:ke))
		allocate (blk(m0)%ggamma   (3,is:ie,js:je,ks:ke))
		!!dist and spacing
		allocate (blk(m0)%dist      (is:ie,js:je,ks:ke))
		allocate (blk(m0)%spci      (is:ie,js:je,ks:ke))
		allocate (blk(m0)%spcj      (is:ie,js:je,ks:ke))
		allocate (blk(m0)%spck      (is:ie,js:je,ks:ke))
		allocate (blk(m0)%vol       (is:ie,js:je,ks:ke))
		allocate (blk(m0)%spacing   (is:ie,js:je,ks:ke))
	end do
	
	if(myid .eq. 0)then
		print *,'***********************************************************'
		print *,"finish allocating memory!"
		print *,'***********************************************************'
		!!pause
	end if
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	!!**********************************************************************!!
	!!                      reading grid                                    !!
	!!**********************************************************************!!
	if(myid .eq. 0)then
		write (*,*)'***********************************************************'
		write (*,*)'reading grid......'
		write (*,*)'***********************************************************' 
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	do m = 1,nblock
		is = global_is(m); ie = global_ie(m)
		js = global_js(m); je = global_je(m)
		ks = global_ks(m); ke = global_ke(m)

		is0 = global_is0(m); ie0 = global_ie0(m)
		js0 = global_js0(m); je0 = global_je0(m)
		ks0 = global_ks0(m); ke0 = global_ke0(m)

		idim = ie0-is0+1
		jdim = je0-js0+1
		kdim = ke0-ks0+1
		
		allocate(tempx(is0:ie0 ,js0:je0, ks0:ke0))
		allocate(tempy(is0:ie0 ,js0:je0, ks0:ke0))
		allocate(tempz(is0:ie0 ,js0:je0, ks0:ke0))

		if (myid .eq. 0)then  		
			if     (iflag_dimension .eq. iflag_2d)then
				read(99,*)(((tempx(i,j,k),i=is0,ie0), j=js0,je0), k=ks0,ke0), &
						  (((tempy(i,j,k),i=is0,ie0), j=js0,je0), k=ks0,ke0)
				tempz = 0.d0

				do k = ks0,ke0
				do j = js0,je0
				do i = is0,ie0
					tempx(i,j,k) = tempx(i,j,k)*scale_x/lref
					tempy(i,j,k) = tempy(i,j,k)*scale_y/lref
				end do
				end do
				end do  
			else if(iflag_dimension .eq. iflag_3d)then
				read(99,*)(((tempx(i,j,k),i=is0,ie0), j=js0,je0), k=ks0,ke0), &
						  (((tempy(i,j,k),i=is0,ie0), j=js0,je0), k=ks0,ke0), &
						  (((tempz(i,j,k),i=is0,ie0), j=js0,je0), k=ks0,ke0)
			          	
				do k = ks0,ke0
				do j = js0,je0
				do i = is0,ie0
					tempx(i,j,k) = tempx(i,j,k)*scale_x/lref
					tempy(i,j,k) = tempy(i,j,k)*scale_y/lref
					tempz(i,j,k) = tempz(i,j,k)*scale_z/lref
				end do
				end do
				end do
				!!*
			end if
		end if
		
		!!id of thread to receive
		call m2node(dest_id,m,myn,lastn,nblock,numprocs)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!!*
		
		if (myid .eq. 0)then	
			if     (dest_id .eq. 0)then
				m0 = m
!!				write(*,15) m,nblock,m0,blk_loop,myid		
!!15			format("reading grid: block",I5 ," of",I5," to block",I5," of",I5," in rank",I4)
				is0 = global_is0(m0); ie0 = global_ie0(m0)
				js0 = global_js0(m0); je0 = global_je0(m0)
				ks0 = global_ks0(m0); ke0 = global_ke0(m0)
				do k = ks0,ke0
					do j = js0,je0
						do i = is0,ie0
							blk(m0)%x(i,j,k) = tempx(i,j,k)
							blk(m0)%y(i,j,k) = tempy(i,j,k)
							blk(m0)%z(i,j,k) = tempz(i,j,k)
						end do 
					end do
				end do
			else if(dest_id .gt. 0)then
				call MPI_send(tempx,idim*jdim*kdim,MPI_REAL8,dest_id,97,MPI_COMM_WORLD,ierr)
				call MPI_send(tempy,idim*jdim*kdim,MPI_REAL8,dest_id,98,MPI_COMM_WORLD,ierr)
				call MPI_send(tempz,idim*jdim*kdim,MPI_REAL8,dest_id,99,MPI_COMM_WORLD,ierr)	      		
			end if
		end if
		
		source_id = 0
		
		if(myid .gt. 0 .and. myid .eq. dest_id)then
			m0 = m - myid*myn
!!			write(*,20) m,nblock,m0,blk_loop,myid
!!20		format("reading grid: block",I5 ," of",I5," to block",I5," of",I5," in rank",I4)
			call MPI_recv(tempx,idim*jdim*kdim,MPI_REAL8,source_id,97,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(tempy,idim*jdim*kdim,MPI_REAL8,source_id,98,MPI_COMM_WORLD,status,ierr)
			call MPI_recv(tempz,idim*jdim*kdim,MPI_REAL8,source_id,99,MPI_COMM_WORLD,status,ierr)

			is0 = global_is0(m0 + myid*myn); ie0 = global_ie0(m0 + myid*myn)
			js0 = global_js0(m0 + myid*myn); je0 = global_je0(m0 + myid*myn)
			ks0 = global_ks0(m0 + myid*myn); ke0 = global_ke0(m0 + myid*myn)
			do k = ks0,ke0
			do j = js0,je0
			do i = is0,ie0
				blk(m0)%x(i,j,k) = tempx(i,j,k)
				blk(m0)%y(i,j,k) = tempy(i,j,k)
				blk(m0)%z(i,j,k) = tempz(i,j,k)
			end do 
			end do
			end do
		end if
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		deallocate(tempx,tempy,tempz)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	end do
	
	do m = 1,nblock
		call MPI_BCAST(blk0(m)%ni,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%nj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		call MPI_BCAST(blk0(m)%is,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ie,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%js,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%je,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ks,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ke,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


		call MPI_BCAST(blk0(m)%is0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ie0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%js0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%je0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ks0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(blk0(m)%ke0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	end do
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
	if(myid .eq. 0)then  
		write(*,*)'***********************************************************'
		write(*,*)'finish reading grid!'
		write(*,*)'***********************************************************'
		close(99)
	end if
	
	deallocate(tempni,tempnj,tempnk)
	deallocate(global_is,  global_ie,  global_js,  global_je,  global_ks,  global_ke)
	deallocate(global_is0, global_ie0, global_js0, global_je0, global_ks0, global_ke0)
	deallocate(global_is1, global_ie1, global_js1, global_je1, global_ks1, global_ke1)
	
	return
end subroutine


