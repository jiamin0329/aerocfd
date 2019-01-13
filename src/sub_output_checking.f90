!!=========================================================!!
!! checking flow field                                     !!
!!                                                         !!
!! coordination: x,y,z                                     !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!! geom variables: spci,spcj,spck                          !!
!!                 spacing,dist,vol                        !!
!!                                                         !!
!! Author: jiamin xu                                       !!
!! Date:   2012.12.05                                      !!
!!=========================================================!!
subroutine checking_output
	use blk_var
	use mpi_var
	use index_var
	implicit none
	
	integer :: i,j,k
	integer :: temp = 0
	character(len = 180) :: flow
	integer :: is, ie, js, je, ks, ke !!block dimension
	!!*********************************************************************!!
	!!geom file output
	
	do m0 = 1,blk_loop
		is = 1; ie = blk(m0)%ni
		js = 1; je = blk(m0)%nj
		ks = 1; ke = blk(m0)%nk
		
		do k = ks,ke
			do j = js,je
				do i = is,ie
					if(blk(m0)%pri_v(6,i,j,k) .lt. 0)then
						write(*,*) "Negative Temperature!!!"
                        write(*,*) i,j,k
						temp = 1
						goto 5
					end if
				end do
			end do
		end do
	end do
	!!*
	
5 if(temp .eq. 1) then
		write(flow,"('result/checking_'I4.4'.dat')"),myid
		open (99,file = flow)
		
		do m0 = 1,blk_loop
			is = 1; ie = blk(m0)%ni
			js = 1; je = blk(m0)%nj
			ks = 1; ke = blk(m0)%nk
			
			write(99,*) "variables = x,y,z,d,u,v,w,p"
			write(99,*) "zone i= ",ie, "j= ",je, "k= ",ke,"datapacking=point"
			
			do k = ks,ke
				do j = js,je
					do i = is,ie
						write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k), &
						            blk(m0)%pri_v(1,i,j,k),                             &
						            blk(m0)%pri_v(2,i,j,k),                             &
						            blk(m0)%pri_v(3,i,j,k),                             &
						            blk(m0)%pri_v(4,i,j,k),                             &
						            blk(m0)%pri_v(5,i,j,k)
					end do
				end do
			end do
		end do
		close(99)
		read(*,*)
	end if
	!!*********************************************************************!!

	return
end subroutine
