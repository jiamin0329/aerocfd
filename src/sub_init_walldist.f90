!!=========================================================!!
!! get wall distance                                       !!
!!                                                         !!
!! purpose: to compute the wall distance in each node      !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.11.22                                      !!
!!=========================================================!!
subroutine init_walldist
	use blk_var
	use mpi_var
	use glbindex_var
	use flag_var
	use index_var
	implicit none
		
	integer :: i,j,k                   !!index of point in box
	real*8  :: xx,yy,zz                !!coordinate of point in the box
	integer :: i0,j0,k0                !!index of point in box
	integer :: is0,ie0,js0,je0,ks0,ke0 !!box dimension
	integer :: is,ie,js,je,ks,ke !!box dimension
	real*8  :: xx0,yy0,zz0             !!coordinate of point in the box
	real*8  :: temp
	real*8  :: walldist         
	character(len = 180):: test
	integer :: idim,jdim,kdim
	real*8,allocatable,dimension(:,:,:) :: tempdist
	integer :: glb_m
	
	
	do m0 = 1, blk_loop
		!!block dimension
		is = blk(m0)%is0;ie = blk(m0)%ie0
		js = blk(m0)%js0;je = blk(m0)%je0
		ks = blk(m0)%ks0;ke = blk(m0)%ke0
		!!*
		print *, "*****************", myid, m0, "of", blk_loop
		!!loop over current block 
		do k = ks,ke
		do j = js,je
		do i = is,ie
			xx = blk(m0)%x(i,j,k)
			yy = blk(m0)%y(i,j,k)
			zz = blk(m0)%z(i,j,k)	

			walldist = 1.d9
			do glb_m = 1,nblock
				do ksub = 1,blk0(glb_m)%num_subface
					if (blk0(glb_m)%bc0(ksub)%blk_t .eq. bc_wall)then
						is0 = blk0(glb_m)%bc0(ksub)%is; ie0 = blk0(glb_m)%bc0(ksub)%ie
						js0 = blk0(glb_m)%bc0(ksub)%js; je0 = blk0(glb_m)%bc0(ksub)%je
						ks0 = blk0(glb_m)%bc0(ksub)%ks; ke0 = blk0(glb_m)%bc0(ksub)%ke

						!!print *, myid, m0, glb_m, ksub, is0, ie0, js0, je0, ks0, ke0
						do k0 = ks0-ks0+1,ke0-ks0+1
						do j0 = js0-js0+1,je0-js0+1
						do i0 = is0-is0+1,ie0-is0+1
							!!print *, blk0(glb_m)%bc0(ksub)%face, i0, j0, k0
							if 		(blk0(glb_m)%bc0(ksub)%face .eq. 1 .or. blk0(glb_m)%bc0(ksub)%face .eq. 2) then
								xx0 = blk0(glb_m)%bc0(ksub)%x(j0,k0)
								yy0 = blk0(glb_m)%bc0(ksub)%y(j0,k0)
								zz0 = blk0(glb_m)%bc0(ksub)%z(j0,k0)
							else if (blk0(glb_m)%bc0(ksub)%face .eq. 3 .or. blk0(glb_m)%bc0(ksub)%face .eq. 4) then
								xx0 = blk0(glb_m)%bc0(ksub)%x(i0,k0)
								yy0 = blk0(glb_m)%bc0(ksub)%y(i0,k0)
								zz0 = blk0(glb_m)%bc0(ksub)%z(i0,k0)
							else if (blk0(glb_m)%bc0(ksub)%face .eq. 5 .or. blk0(glb_m)%bc0(ksub)%face .eq. 6) then
								xx0 = blk0(glb_m)%bc0(ksub)%x(i0,j0)
								yy0 = blk0(glb_m)%bc0(ksub)%y(i0,j0)
								zz0 = blk0(glb_m)%bc0(ksub)%z(i0,j0)
							end if
							!!if (i .eq. 10 .and. j .eq. 10) then
							!!	print *, myid, m0, glb_m, ksub, xx0, yy0, zz0
							!!end if

							temp = sqrt((xx-xx0)**2 + (yy-yy0)**2 + (zz-zz0)**2)													
							if (temp .le. walldist)then 
								walldist = temp
							end if
						end do
						end do
						end do
					end if
				end do
			end do

			blk(m0)%dist(i,j,k) = walldist
		end do
		end do
		end do
	end do

	!!*********************************************************************!!
	!!geom file output
	write(test,"('result/dist_'I5.5'.dat')"),myid
	if (myid .eq. root)then
		write(*,*) "writing dist_debug parameters......"
	end if
	
	open (99,file = test)
	
	do m0 = 1,blk_loop
		is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
		js0 = blk(m0)%js0;je0 = blk(m0)%je0
		ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0
		
		write(99,*) "variables = x,y,z,dist"
		write(99,*) "zone i= ", ie0-is0+1, "j= ", je0-js0+1, "k= ", ke0-ks0+1, "datapacking=point"
		
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),blk(m0)%dist(i,j,k)
		end do
		end do
		end do
	end do
	close(99)
	!!*********************************************************************!!
	return
end subroutine


