!!=========================================================!!
!!  residual calculation                                   !!
!!                                                         !!
!!                                                         !!
!!  author: jiamin xu                                      !!
!!  date:   2012.12.24                                     !!
!!=========================================================!!
subroutine residual
	use blk_var
	use flag_var
	use ns_const
	use mpi_var
	use index_var
	implicit none
    
	integer :: i,j,k
	integer :: is1, ie1, js1, je1, ks1, ke1
	real*8  :: resi_rms(5)
	real*8  :: resi_tur(2)
	real*8  :: temp(5),temp_tur(2)
	
	character(len = 180) :: resi_out
	
	resi_rms = 0.d0
	resi_tur = 0.d0
	!!main equation residual
	do m0 = 1,blk_loop
		is1 = blk(m0)%is1; ie1 = blk(m0)%ie1
		js1 = blk(m0)%js1; je1 = blk(m0)%je1
		ks1 = blk(m0)%ks1; ke1 = blk(m0)%ke1
			
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			resi_rms(1) = resi_rms(1) + (blk(m0)%rhs(1,i,j,k)*blk(m0)%dt(i,j,k)/blk(m0)%inv_j(i,j,k))**2
			resi_rms(2) = resi_rms(2) + (blk(m0)%rhs(2,i,j,k)*blk(m0)%dt(i,j,k)/blk(m0)%inv_j(i,j,k))**2
			resi_rms(3) = resi_rms(3) + (blk(m0)%rhs(3,i,j,k)*blk(m0)%dt(i,j,k)/blk(m0)%inv_j(i,j,k))**2
			resi_rms(4) = resi_rms(4) + (blk(m0)%rhs(4,i,j,k)*blk(m0)%dt(i,j,k)/blk(m0)%inv_j(i,j,k))**2
			resi_rms(5) = resi_rms(5) + (blk(m0)%rhs(5,i,j,k)*blk(m0)%dt(i,j,k)/blk(m0)%inv_j(i,j,k))**2
					
			resi_tur(1) = resi_tur(1) + (blk(m0)%sa_rhs(i,j,k)*blk(m0)%dt(i,j,k))**2
		end do
		end do
		end do
		!!*
	end do

	!!reduction
	call MPI_Reduce(resi_rms,temp,    5,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_Reduce(resi_tur,temp_tur,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	!!*
	temp     = sqrt(temp    /totalnodes)
	temp_tur = sqrt(temp_tur/totalnodes)
	!!output to screen
	if(myid .eq. root)then
		print *, 'residual:'
		write(*,5) temp(1),temp(2),temp(3),temp(5)
5	format("D: ",E11.4," U: ",E11.4," V: ",E11.4," E: ",E11.4)
		write(*,10) temp_tur(1)
10	format("nut: ",E11.4)
	end if
    
	return
end subroutine

