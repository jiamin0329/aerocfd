!!=========================================================
!! physical boundary condition
!!
!! purpose: to set physical boundary conditions
!!
!! AUTHOR: Jiamin XU                                       
!! DATE:   2012.11.28
!!=========================================================
subroutine physical_bc
	use ns_const
	use blk_var
	use index_var
	use mpi_var
	use flag_var
	use sa_var
	implicit none
    
	integer :: face
	integer :: is, ie, js, je, ks, ke
	integer :: is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc

	real*8 :: d,nut
	integer :: i,j,k
	
	do m0 = 1,blk_loop
		do ksub = 1,blk(m0)%num_subface
			!!face type
			face  = blk(m0)%bc(ksub)%face
			!!subface dimension
			is_bc = blk(m0)%bc(ksub)%is; ie_bc = blk(m0)%bc(ksub)%ie
			js_bc = blk(m0)%bc(ksub)%js; je_bc = blk(m0)%bc(ksub)%je
			ks_bc = blk(m0)%bc(ksub)%ks; ke_bc = blk(m0)%bc(ksub)%ke
			!!block dimension
			is = blk(m0)%is; ie = blk(m0)%ie
			js = blk(m0)%js; je = blk(m0)%je
			ks = blk(m0)%ks; ke = blk(m0)%ke
			
			if      (blk(m0)%bc(ksub)%blk_t .eq. bc_user    )then
			!!	call userbc (blk(m0)%q,blk(m0)%x,blk(m0)%y,blk(m0)%z, &
			!!	             dt,timestep,d_0,u_0,p_0,gamma,ma,cv,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_periodic)then
				print *, "periodic bc is not available now!"
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_inflow  )then
				call inflow (blk(m0)%q,blk(m0)%nut,d_0,u_0,p_0,aoa,gamma,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_outflow )then
				call outflow(blk(m0)%q,blk(m0)%nut,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_symm    )then
				call symm(blk(m0)%q,blk(m0)%nut,blk(m0)%dxidx,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_wall    )then
				if      (iflag_solver .eq. iflag_nssolver   )then
					call nonslipwall (blk(m0)%q,blk(m0)%nut,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
				else if (iflag_solver .eq. iflag_eulersolver)then
					call slipwall(blk(m0)%q,blk(m0)%nut,blk(m0)%dxidx,face,is_bc,ie_bc,js_bc,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
				end if
			else if (blk(m0)%bc(ksub)%blk_t .eq. bc_farfield)then
				call farfield (blk(m0)%q,blk(m0)%nut,blk(m0)%dxidx, &
				               d_0,u_0,p_0,c_0,aoa,gamma,face,         &
				               is_bc,ie_bc,js,je_bc,ks_bc,ke_bc,is,ie,js,je,ks,ke)
			end if
		end do
	end do


	do m0 = 1,blk_loop
		is = blk(m0)%is; ie = blk(m0)%ie
		js = blk(m0)%js; je = blk(m0)%je
		ks = blk(m0)%ks; ke = blk(m0)%ke
		
		do k = ks,ke
		do j = js,je
		do i = is,ie
			d    = blk(m0)%pri_v(1,i,j,k)
			nut  = blk(m0)%nut  (i,j,k)
			leix = blk(m0)%nut  (i,j,k)/blk(m0)%amu(1,i,j,k)*d
			leix = max(leix,0.001d0)
			fv1  = leix**3/(leix**3 + cv1**3)
					
			blk(m0)%nut(i,j,k)   = max(nut,0.0001d0)
			blk(m0)%nut(i,j,k)   = min(nut,1.d6)					
			blk(m0)%amu(2,i,j,k) = blk(m0)%nut(i,j,k)*fv1*d
		end do
		end do
		end do	
	end do    

	return
end subroutine





