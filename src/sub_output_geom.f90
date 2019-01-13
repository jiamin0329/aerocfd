!!=========================================================!!
!! geomtry output                                          !!
!!                                                         !!
!! coordination: x,y,z                                     !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!! geom variables: spci,spcj,spck                          !!
!!                 spacing,dist,vol                        !!
!!                                                         !!
!! Author: jiamin xu                                       !!
!! Date:   2012.12.05                                      !!
!!=========================================================!!
subroutine geom_output
	use blk_var
	use mpi_var
	use index_var
	implicit none
	
	integer :: i,j,k
	character(len = 180) :: geom
	integer :: is, ie, js, je, ks, ke !!block dimension
	integer :: idim, jdim, kdim
	!!*********************************************************************!!
	!!geom file output
	write(geom,"('result/geom_'I2.2'.dat')"),myid
	if (myid .eq. root) then
		write(*,*) "writing geom parameters......"
	end if
	open (99,file = geom)
	
	do m0 = 1,blk_loop
		is = blk(m0)%is0;ie = blk(m0)%ie0
		js = blk(m0)%js0;je = blk(m0)%je0
		ks = blk(m0)%ks0;ke = blk(m0)%ke0

		idim = ie-is+1
		jdim = je-js+1
		kdim = ke-ks+1

		write(99,*) "variables = x,y,z,spacing,vol"
		write(99,*) "zone i= ",idim, "j= ",jdim, "k= ",kdim,"datapacking=point"
		
		do k = ks,ke
		do j = js,je
		do i = is,ie
			write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k),blk(m0)%spacing(i,j,k),blk(m0)%vol(i,j,k)
		end do
		end do
		end do
	end do
	
	close(99)
	!!*********************************************************************!!

	return
end subroutine
