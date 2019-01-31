!!*************************************************************
!! inflow boundary condition
!!
!! author: Jiamin Xu
!! date:   2013.08.27
!!*************************************************************
subroutine inflow (q,tur,d_0,u_0,p_0,aoa,gamma, &
	               face,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use flag_var
	use sa_var
	implicit none
	
	integer :: i,j,k                   !!index
	integer :: is, ie, js, je, ks, ke  !!subface dimension
	integer :: is0,ie0,js0,je0,ks0,ke0 !!block dimension
	integer :: face
	real*8  :: q(5,is0:ie0,js0:je0,ks0:ke0)
    real*8  :: tur(is0:ie0,js0:je0,ks0:ke0)
    
	real*8  :: gamma
	real*8  :: d_0,u_0,p_0
	real*8  :: aoa

    
	!!subface i-
	if (face .eq. 1) then
		do k = ks,ke
			do j = js,je
				q(1,is,j,k) = d_0
				q(2,is,j,k) = d_0*u_0*cos(aoa)
				q(3,is,j,k) = d_0*u_0*sin(aoa)
				q(4,is,j,k) = 0.d0
				q(5,is,j,k) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				if (iflag_turbulence .ge. 1) then
					tur(is,j,k) = nutinf
				end if
				
			end do
		end do									
	end if
	!!*end subface i-
	
	!!subface i+
	if (face .eq. 2) then
		do k = ks,ke
			do j = js,je
				q(1,ie,j,k) = d_0                            
				q(2,ie,j,k) = d_0*u_0*cos(aoa)                   
				q(3,ie,j,k) = d_0*u_0*sin(aoa)  
				q(4,ie,j,k) = 0.d0                    
				q(5,ie,j,k) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				
				if (iflag_turbulence .ge. 1) then
					tur(ie,j,k) = nutinf
				end if
			end do
		end do
	end if
	!!*end subface i+
	
	!!subface j-
	if (face .eq. 3) then
		do k = ks,ke
			do i = is,ie
				q(1,i,js,k) = d_0                            
				q(2,i,js,k) = d_0*u_0*cos(aoa)                   
				q(3,i,js,k) = d_0*u_0*sin(aoa)   
				q(4,i,js,k) = 0.d0                
				q(5,i,js,k) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				
				if (iflag_turbulence .ge. 1) then
					tur(i,js,k) = nutinf
				end if
			end do
		end do	
	end if
	!!*end subface j-
	
	!!subface j+
	if (face .eq. 4) then
		do k = ks,ke
			do i = is,ie
				q(1,i,je,k) = d_0                            
				q(2,i,je,k) = d_0*u_0*cos(aoa)                   
				q(3,i,je,k) = d_0*u_0*sin(aoa)
				q(4,i,je,k) = 0.d0
				q(5,i,je,k) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				
				if (iflag_turbulence .ge. 1) then
					tur(i,je,k) = nutinf
				end if
			end do
		end do			
	end if
	!!*end subface j+    

	!!subface k-
	if (face .eq. 5) then
		do j = js,je
			do i = is,ie
				q(1,i,j,ks) = d_0                            
				q(2,i,j,ks) = d_0*u_0*cos(aoa)                   
				q(3,i,j,ks) = d_0*u_0*sin(aoa)   
				q(4,i,j,ks) = 0.d0                
				q(5,i,j,ks) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ks) = nutinf
				end if
			end do
		end do	
	end if
	!!*end subface k-
	
	!!subface k+
	if (face .eq. 6) then
		do j = js,je
			do i = is,ie
				q(1,i,j,ke) = d_0                            
				q(2,i,j,ke) = d_0*u_0*cos(aoa)                   
				q(3,i,j,ke) = d_0*u_0*sin(aoa)
				q(4,i,j,ke) = 0.d0
				q(5,i,j,ke) = p_0/(gamma-1.d0) + 0.5d0*u_0**2

				
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ke) = nutinf
				end if
			end do
		end do			
	end if
	!!*end subface k+ 
		
	return
end subroutine