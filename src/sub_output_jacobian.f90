!!=========================================================!!
!! jacobian metrix output                                  !!
!!                                                         !!
!! coordination: x,y,z                                     !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!! geom variables: spci,spcj,spck                          !!
!!                 spacing,dist,vol                        !!
!!                                                         !!
!! Author: jiamin xu                                       !!
!! Date:   2012.12.05                                      !!
!!=========================================================!!
subroutine jacobian_output
	use blk_var
	use mpi_var
	use index_var
	implicit none
	
	integer :: i,j,k
	integer :: is, ie, js, je, ks, ke !!block dimension
	integer :: idim, jdim, kdim
	character(len = 180) :: geom
	!!*********************************************************************!!
	!!geom file output
	write(geom,"('result/jacobian_'I2.2'.dat')"),myid
	if (myid .eq. root)then
		write(*,*) "writing jacobian parameters......"
	end if
	
	open (99,file = geom)
	
	do m0 = 1,blk_loop
		is = blk(m0)%is0; ie = blk(m0)%ie0
		js = blk(m0)%js0; je = blk(m0)%je0
		ks = blk(m0)%ks0; ke = blk(m0)%ke0
		idim = ie-is+1
		jdim = je-js+1
		kdim = ke-ks+1

		write(99,*) "variables = x,y,z,dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz,inv_j,"
		write(99,*) "            alpha1, alpha2, alpha3, beta1, beta2, beta3, gamma1, gamma2, gamma3"
		write(99,*) "zone i= ",idim, "j= ",jdim, "k= ",kdim,"datapacking=point"
		
		do k = ks,ke
		do j = js,je
		do i = is,ie
			write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k), &
					    blk(m0)%dxidx(1,i,j,k),blk(m0)%dxidx(2,i,j,k),blk(m0)%dxidx(3,i,j,k),&
					    blk(m0)%dxidx(4,i,j,k),blk(m0)%dxidx(5,i,j,k),blk(m0)%dxidx(6,i,j,k),&
					    blk(m0)%dxidx(7,i,j,k),blk(m0)%dxidx(8,i,j,k),blk(m0)%dxidx(9,i,j,k),&
						blk(m0)%inv_j(i,j,k), &
						blk(m0)%aalpha(1,i,j,k),blk(m0)%aalpha(2,i,j,k),blk(m0)%aalpha(3,i,j,k), &
						blk(m0)%bbeta (1,i,j,k),blk(m0)%bbeta (2,i,j,k),blk(m0)%bbeta (3,i,j,k), &
						blk(m0)%ggamma(1,i,j,k),blk(m0)%ggamma(3,i,j,k),blk(m0)%ggamma(3,i,j,k)
		end do
		end do
		end do
	end do
	close(99)
	!!*********************************************************************!!

	return
end subroutine
