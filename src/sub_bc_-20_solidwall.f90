!!*************************************************************
!! solid wall boundary condition
!!
!! (1) slip wall for euler solver
!! (2) nonslip wall for ns solver
!!
!! author: Jiamin Xu
!! date  : 2018.08.27
!!*************************************************************
subroutine slipwall(q,tur,dxidx,face,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use flag_var
	implicit none

	integer :: i,j,k                   !!index
	integer :: is, ie, js, je, ks, ke  !!subface dimension
	integer :: is0,ie0,js0,je0,ks0,ke0 !!block dimension
	integer :: face
	
	real*8  :: q      (5,is0:ie0,js0:je0,ks0:ke0)
	real*8  :: dxidx  (9,is0:ie0,js0:je0,ks0:ke0)
	real*8  :: tur    (  is0:ie0,js0:je0,ks0:ke0)
    
	real*8  :: dxidx0,dxidy0,dxidz0,temp
	real*8  :: cu
	
	!!subface i-
	if (face .eq. 1) then
		do k = ks,ke
			do j = js,je
				dxidx0 = dxidx(1,is,j,k)
				dxidy0 = dxidx(2,is,j,k)
				dxidz0 = dxidx(3,is,j,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 =-dxidx0/temp
				dxidy0 =-dxidy0/temp
				dxidz0 =-dxidz0/temp
				
				cu = q(2,is+1,j,k)*dxidx0 + q(3,is+1,j,k)*dxidy0 + q(4,is+1,j,k)*dxidz0
				
				q(1,is,j,k) = q(1,is+1,j,k)
				q(2,is,j,k) = q(2,is+1,j,k) - cu*dxidx0
				q(3,is,j,k) = q(3,is+1,j,k) - cu*dxidy0
				q(4,is,j,k) = q(4,is+1,j,k) - cu*dxidz0
				q(5,is,j,k) = q(5,is+1,j,k)		
				if (iflag_turbulence .ge. 1) then
					tur(is,j,k) = tur(is+1,j,k)
				end if
			end do						
		end do	
	end if
	!!*end subface i-
	
	!!subface i+
	if (face .eq. 2) then
		do k = ks,ke
			do j = js,je
				dxidx0 = dxidx(1,ie,j,k)
				dxidy0 = dxidx(2,ie,j,k)
				dxidz0 = dxidx(3,ie,j,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 = dxidx0/temp
				dxidy0 = dxidy0/temp
				dxidz0 = dxidz0/temp
				
				cu = q(2,ie-1,j,k)*dxidx0 + q(3,ie-1,j,k)*dxidy0 + q(4,ie-1,j,k)*dxidz0
				
				q(1,ie,j,k) = q(1,ie-1,j,k)
				q(2,ie,j,k) = q(2,ie-1,j,k) - cu*dxidx0
				q(3,ie,j,k) = q(3,ie-1,j,k) - cu*dxidy0
				q(4,ie,j,k) = q(4,ie-1,j,k) - cu*dxidz0
				q(5,ie,j,k) = q(5,ie-1,j,k)	

				if (iflag_turbulence .ge. 1) then
					tur(ie,j,k) = tur(ie-1,j,k)
				end if
			end do						
		end do	
	end if
	!!*end subface i+
	
	!!subface j-
	if (face .eq. 3) then
		do k = ks,ke
			do i = is,ie
				dxidx0 = dxidx(4,i,js,k)
				dxidy0 = dxidx(5,i,js,k)
				dxidz0 = dxidx(6,i,js,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 =-dxidx0/temp
				dxidy0 =-dxidy0/temp
				dxidz0 =-dxidz0/temp
				
				cu = q(2,i,js+1,k)*dxidx0 + q(3,i,js+1,k)*dxidy0 + q(4,i,js+1,k)*dxidz0
				
				q(1,i,js,k) = q(1,i,js+1,k)
				q(2,i,js,k) = q(2,i,js+1,k) - cu*dxidx0
				q(3,i,js,k) = q(3,i,js+1,k) - cu*dxidy0
				q(4,i,js,k) = q(4,i,js+1,k)	- cu*dxidz0
				q(5,i,js,k) = q(5,i,js+1,k)

				if (iflag_turbulence .ge. 1) then
					tur(i,js,k) = tur(i,js+1,k)
				end if
			end do								
		end do	
	end if
	!!*end subface j-
	
	!!subface j+
	if (face .eq. 4) then
		do k = ks,ke
			do i = is,ie
				dxidx0 = dxidx(4,i,je,k)
				dxidy0 = dxidx(5,i,je,k)
				dxidz0 = dxidx(6,i,je,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 = dxidx0/temp
				dxidy0 = dxidy0/temp
				dxidz0 = dxidz0/temp
				
				cu = q(2,i,je-1,k)*dxidx0 + q(3,i,je-1,k)*dxidy0 + q(4,i,je-1,k)*dxidz0
				
				q(1,i,je,k) = q(1,i,je-1,k)
				q(2,i,je,k) = q(2,i,je-1,k) - cu*dxidx0
				q(3,i,je,k) = q(3,i,je-1,k) - cu*dxidy0
				q(4,i,je,k) = q(4,i,je-1,k) - cu*dxidz0
				q(5,i,je,k) = q(5,i,je-1,k)

				if (iflag_turbulence .ge. 1) then
					tur(i,je,k) = tur(i,je-1,k)
				end if
			end do							
		end do			
	end if
	!!*end subface j+
	
	!!subface k-
	if (face .eq. 5) then
		do j = js,je
			do i = is,ie
				dxidx0 = dxidx(7,i,j,ks)
				dxidy0 = dxidx(8,i,j,ks)
				dxidz0 = dxidx(9,i,j,ks)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 =-dxidx0/temp
				dxidy0 =-dxidy0/temp
				dxidz0 =-dxidz0/temp
				
				cu = q(2,i,j,ks+1)*dxidx0 + q(3,i,j,ks+1)*dxidy0 + q(4,i,j,ks+1)*dxidz0
				
				q(1,i,j,ks) = q(1,i,j,ks+1)
				q(2,i,j,ks) = q(2,i,j,ks+1) - cu*dxidx0
				q(3,i,j,ks) = q(3,i,j,ks+1) - cu*dxidy0
				q(4,i,j,ks) = q(4,i,j,ks+1) - cu*dxidz0
				q(5,i,j,ks) = q(5,i,j,ks+1)

				if (iflag_turbulence .ge. 1) then
					tur(i,j,ks) = tur(i,j,ks+1)
				end if
			end do							
		end do			
	end if
	!!*end subface k-
	
	!!subface k+
	if (face .eq. 6) then
		do j = js,je
			do i = is,ie
				dxidx0 = dxidx(7,i,j,ke)
				dxidy0 = dxidx(8,i,j,ke)
				dxidz0 = dxidx(9,i,j,ke)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 = dxidx0/temp
				dxidy0 = dxidy0/temp
				dxidz0 = dxidz0/temp
				
				cu = q(2,i,j,ke-1)*dxidx0 + q(3,i,j,ke-1)*dxidy0 + q(4,i,j,ke-1)*dxidz0
				
				q(1,i,j,ke) = q(1,i,j,ke-1)
				q(2,i,j,ke) = q(2,i,j,ke-1) - cu*dxidx0
				q(3,i,j,ke) = q(3,i,j,ke-1) - cu*dxidy0
				q(4,i,j,ke) = q(4,i,j,ke-1) - cu*dxidz0
				q(5,i,j,ke) = q(5,i,j,ke-1)
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ke) = tur(i,j,ke-1)
				end if
			end do							
		end do			
	end if
	!!*end subface k+

	return
end subroutine

subroutine nonslipwall (q,tur,face,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use flag_var
	implicit none
	
	integer :: i,j,k                       !!index
	integer :: is, ie, js, je, ks, ke      !!subface dimension
	integer :: is0,ie0,js0,je0,ks0,ke0     !!block dimension
	integer :: face                        !!face type(1-6)
	real*8  :: q(5,is0:ie0,js0:je0,ks0:ke0)!!conservative variables
	real*8  :: tur(is0:ie0,js0:je0,ks0:ke0)!!conservative variables
	
	!!subface i-
	if (face .eq. 1) then
		do k = ks,ke
			do j = js,je
				q(1,is,j,k) = (11.d0*q(1,is+1,j,k) - 7.d0*q(1,is+2,j,k) + 2.d0*q(1,is+3,j,k))/6.d0
				q(2,is,j,k) = 0.d0
				q(3,is,j,k) = 0.d0
				q(4,is,j,k) = 0.d0
				q(5,is,j,k) = (11.d0*q(5,is+1,j,k) - 7.d0*q(5,is+2,j,k) + 2.d0*q(5,is+3,j,k))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(is,j,k) = 0.d0
				end if
			end do
		end do
	end if
	!!*end subface i-
	
	!!subface i+
	if (face .eq. 2) then
		do k = ks,ke
			do j = js,je
				q(1,ie,j,k) = (11.d0*q(1,ie-1,j,k) - 7.d0*q(1,ie-2,j,k) + 2.d0*q(1,ie-3,j,k))/6.d0
				q(2,ie,j,k) = 0.d0
				q(3,ie,j,k) = 0.d0
				q(4,ie,j,k) = 0.d0
				q(5,ie,j,k) = (11.d0*q(5,ie-1,j,k) - 7.d0*q(5,ie-2,j,k) + 2.d0*q(5,ie-3,j,k))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(ie,j,k) = 0.d0
				end if
			end do
		end do	
	end if
	!!*end subface i+
	
	!!subface j-
	if (face .eq. 3) then
		do k = ks,ke
			do i = is,ie
				q(1,i,js,k) = (11.d0*q(1,i,js+1,k) - 7.d0*q(1,i,js+2,k) + 2.d0*q(1,i,js+3,k))/6.d0
				q(2,i,js,k) = 0.d0
				q(3,i,js,k) = 0.d0
				q(4,i,js,k) = 0.d0
				q(5,i,js,k) = (11.d0*q(5,i,js+1,k) - 7.d0*q(5,i,js+2,k) + 2.d0*q(5,i,js+3,k))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(i,js,k) = 0.d0
				end if
			end do
		end do		
	end if
	!!*end subface j-
	
	!!subface j+
	if (face .eq. 4) then
		do k = ks,ke
			do i = is,ie
				q(1,i,je,k) = (11.d0*q(1,i,je-1,k) - 7.d0*q(1,i,je-2,k) + 2.d0*q(1,i,je-3,k))/6.d0
				q(2,i,je,k) = 0.d0
				q(3,i,je,k) = 0.d0
				q(4,i,je,k) = 0.d0
				q(5,i,je,k) = (11.d0*q(5,i,je-1,k) - 7.d0*q(5,i,je-2,k) + 2.d0*q(5,i,je-3,k))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(i,je,k) = 0.d0
				end if
			end do
		end do	
	end if
	!!*end subface j+

	!!subface k-
	if (face .eq. 5) then
		do j = js,je
			do i = is,ie
				q(1,i,j,ks) = (11.d0*q(1,i,j,ks+1) - 7.d0*q(1,i,j,ks+2) + 2.d0*q(1,i,j,ks+3))/6.d0
				q(2,i,j,ks) = 0.d0
				q(3,i,j,ks) = 0.d0
				q(4,i,j,ks) = 0.d0
				q(5,i,j,ks) = (11.d0*q(5,i,j,ks+1) - 7.d0*q(5,i,j,ks+2) + 2.d0*q(5,i,j,ks+3))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ks) = 0.d0
				end if
			end do
		end do		
	end if
	!!*end subface k-
	
	!!subface k+
	if (face .eq. 6) then
		do j = js,je
			do i = is,ie
				q(1,i,j,ke) = (11.d0*q(1,i,j,ke-1) - 7.d0*q(1,i,j,ke-2) + 2.d0*q(1,i,j,ke-3))/6.d0
				q(2,i,j,ke) = 0.d0                                                                              
				q(3,i,j,ke) = 0.d0                                                                              
				q(4,i,j,ke) = 0.d0                                                                              
				q(5,i,j,ke) = (11.d0*q(5,i,j,ke-1) - 7.d0*q(5,i,j,ke-2) + 2.d0*q(5,i,j,ke-3))/6.d0
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ke) = 0.d0
				end if
			end do
		end do	
	end if
	!!*end subface k+
			
	return
end subroutine