!!=========================================================!!
!! cp&cf output                                            !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.12.05                                      !!
!!=========================================================!!
subroutine surface_output
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
  implicit none
  
  integer :: idim,jdim,kdim
  real*8  :: qinf
  real*8  :: p,p0
  real*8  :: prscoeff
  integer :: i,j,k
  integer :: is, ie, js, je, ks, ke !!block dimension
  character(len = 180) :: surf
  
  p0   = d_0*t_0/gamma/ma/ma
  qinf = 0.5d0*d_0*u_0**2
  
  if(myid .eq. root)then
		write(*,*) "writing surface data......"
	end if
  write(surf,"('result/surfcp_'I6.6'_'I4.4'.dat')"),timestep,myid
  open(99, file = surf)

  do m0 = 1,blk_loop
  	do ksub = 1,blk(m0)%num_subface
  		if(blk(m0)%bc(ksub)%blk_t .eq. bc_wall)then
  			is = blk(m0)%bc(ksub)%is
  			ie = blk(m0)%bc(ksub)%ie
  			js = blk(m0)%bc(ksub)%js
  			je = blk(m0)%bc(ksub)%je
  			ks = blk(m0)%bc(ksub)%ks
  			ke = blk(m0)%bc(ksub)%ke
  			
  			idim = ie-is+1
  			jdim = je-js+1
  			kdim = ke-ks+1
  			
  			write(99,*) "variables = x,y,z,cp"
				write(99,*) "zone i= ",idim, "j= ",jdim, "k= ",kdim,"datapacking=point"

				if     (blk(m0)%bc(ksub)%face .eq. 1)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				else if(blk(m0)%bc(ksub)%face .eq. 2)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				else if(blk(m0)%bc(ksub)%face .eq. 3)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				else if(blk(m0)%bc(ksub)%face .eq. 4)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				else if(blk(m0)%bc(ksub)%face .eq. 5)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				else if(blk(m0)%bc(ksub)%face .eq. 6)then
					do k = ks,ke
						do j = js,je
							do i = is,ie
								p = blk(m0)%pri_v(5,i,j,k)
								prscoeff = (p-p0)/qinf
								write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),prscoeff
							end do
						end do
					end do
				end if	
			end if
		end do
	end do
	close(99)
	!!*
	
	return
end subroutine
