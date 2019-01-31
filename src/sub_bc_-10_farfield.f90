!!*************************************************************
!! farfield boundary condition
!!
!! author: Jiamin Xu
!! date:   2013.08.27
!!*************************************************************
subroutine farfield (q,tur,dxidx,d_0,u_0,p_0,c_0,aoa,gamma,face,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use sa_var
	use flag_var
	implicit none
	
	integer :: i,j,k                       !!index
	integer :: is, ie, js, je, ks, ke      !!subface dimension
	integer :: is0,ie0,js0,je0,ks0,ke0     !!block dimension
	integer :: face                        !!face type(1-6)
	real*8  :: q(5,is0:ie0,js0:je0,ks0:ke0)!!conservative variables
	real*8  :: tur(is0:ie0,js0:je0,ks0:ke0)!!conservative variables
	real*8  :: d_0,u_0,p_0,c_0,aoa,gamma
	real*8  :: dxidx  (9,is0:ie0,js0:je0,ks0:ke0)
	
	real*8 :: dxidx0, dxidy0, dxidz0, temp
	real*8 :: d_inf, u_inf, v_inf, w_inf, p_inf, s_inf, c_inf, vn_inf
	real*8 :: d_p,   u_p,   v_p,   w_p,   p_p,   s_p,   c_p,   vn_p
	real*8 :: d_b,   u_b,   v_b,   w_b,   p_b,   s_b,   c_b,   vn_b 
	real*8 :: u_ref, v_ref, w_ref, u_ref_t
	real*8 :: turVal
	real*8 :: rp, rm
	
	!!infinite flow condition
	d_inf = d_0
	u_inf = u_0*cos(aoa)
	v_inf = u_0*sin(aoa)
	w_inf = 0.d0
	p_inf = p_0
	s_inf = p_0/d_0**gamma
	c_inf = c_0
	!!*
	
	!!subface i-
	if (face .eq. 1) then
		do k = ks,ke
			do j = js,je
				dxidx0 = dxidx(1,is,j,k)
				dxidy0 = dxidx(2,is,j,k)
				dxidz0 = dxidx(3,is,j,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 = -dxidx0/temp
				dxidy0 = -dxidy0/temp
				dxidz0 = -dxidz0/temp
				
				!!interior point
				d_p = q(1,is+1,j,k)
				u_p = q(2,is+1,j,k)/d_p
				v_p = q(3,is+1,j,k)/d_p
				w_p = q(4,is+1,j,k)/d_p
				p_p =(q(5,is+1,j,k) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p

				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+

				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)

				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(is+1,j,k)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,is,j,k) = d_b
				q(2,is,j,k) = d_b*u_b
				q(3,is,j,k) = d_b*v_b
				q(4,is,j,k) = d_b*w_b
				q(5,is,j,k) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(is,j,k) = turVal
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
				
				!!interior point
				d_p = q(1,ie-1,j,k)
				u_p = q(2,ie-1,j,k)/d_p
				v_p = q(3,ie-1,j,k)/d_p
				w_p = q(4,ie-1,j,k)/d_p
				p_p =(q(5,ie-1,j,k) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p
				
				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+
				
				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)
				
				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(ie-1,j,k)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,ie+0,j,k) = d_b
				q(2,ie+0,j,k) = d_b*u_b
				q(3,ie+0,j,k) = d_b*v_b
				q(4,ie+0,j,k) = d_b*w_b
				q(5,ie+0,j,k) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(ie,j,k) = turVal
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
				
				dxidx0 = -dxidx0/temp
				dxidy0 = -dxidy0/temp
				dxidz0 = -dxidz0/temp
				
				!!interior point
				d_p = q(1,i,js+1,k)
				u_p = q(2,i,js+1,k)/d_p
				v_p = q(3,i,js+1,k)/d_p
				w_p = q(4,i,js+1,k)/d_p
				p_p =(q(5,i,js+1,k) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p

				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+

				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)

				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(i,js+1,k)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,i,js,k) = d_b
				q(2,i,js,k) = d_b*u_b
				q(3,i,js,k) = d_b*v_b
				q(4,i,js,k) = d_b*w_b
				q(5,i,js,k) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(i,js,k) = turVal
				end if
			end do
		end do	
	end if
	!!end subface j-
	
	!!subface j+
	if (face .eq. 4) then
		do k = ks, ke
			do i = is, ie
				dxidx0 = dxidx(4,i,je,k)
				dxidy0 = dxidx(5,i,je,k)
				dxidz0 = dxidx(6,i,je,k)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 = dxidx0/temp
				dxidy0 = dxidy0/temp
				dxidz0 = dxidz0/temp
				
				!! interior point
				d_p = q(1,i,je-1,k)
				u_p = q(2,i,je-1,k)/d_p
				v_p = q(3,i,je-1,k)/d_p
				w_p = q(4,i,je-1,k)/d_p
				p_p =(q(5,i,je-1,k) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p
				
				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+
				
				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)
				
				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(i,je-1,k)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,i,je,k) = d_b
				q(2,i,je,k) = d_b*u_b
				q(3,i,je,k) = d_b*v_b
				q(4,i,je,k) = d_b*w_b
				q(5,i,je,k) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(i,je,k) = turVal
				end if
			end do
		end do	
	end if
	!!*end subface j+
	
	!!subface k-
	if (face .eq. 5) then
		do j = js, je
			do i = is, ie
				dxidx0 = dxidx(7,i,j,ks)
				dxidy0 = dxidx(8,i,j,ks)
				dxidz0 = dxidx(9,i,j,ks)
				
				temp   = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp   = sqrt(temp)
				
				dxidx0 =-dxidx0/temp
				dxidy0 =-dxidy0/temp
				dxidz0 =-dxidz0/temp
				
				!! interior point
				d_p = q(1,i,j,ks+1)
				u_p = q(2,i,j,ks+1)/d_p
				v_p = q(3,i,j,ks+1)/d_p
				w_p = q(4,i,j,ks+1)/d_p
				p_p =(q(5,i,j,ks+1) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p
				
				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+
				
				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)
				
				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(i,j,ks+1)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,i,j,ks) = d_b
				q(2,i,j,ks) = d_b*u_b
				q(3,i,j,ks) = d_b*v_b
				q(4,i,j,ks) = d_b*w_b
				q(5,i,j,ks) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ks) = turVal
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
				
				!!interior point
				d_p = q(1,i,j,ke-1)
				u_p = q(2,i,j,ke-1)/d_p
				v_p = q(3,i,j,ke-1)/d_p
				w_p = q(4,i,j,ke-1)/d_p
				p_p =(q(5,i,j,ke-1) - 0.5d0*d_p*(u_p**2 + v_p**2 + w_p**2))*(gamma-1.d0)
				c_p = sqrt(gamma*p_p/d_p)
				s_p = p_p/(d_p**gamma)
				
				vn_inf = u_inf*dxidx0 + v_inf*dxidy0 + w_inf*dxidz0   !!* (vn)??
				vn_p   = u_p  *dxidx0 + v_p  *dxidy0 + w_p  *dxidz0   !!* (vn)p
				
				rm     = vn_inf - 2.d0*c_inf/(gamma-1.d0)             !!* r-
				rp     = vn_p   + 2.d0*c_p  /(gamma-1.d0)             !!* r+
				
				vn_b   = 0.5d0*(rp + rm)
				c_b    = 0.25d0*(gamma-1.d0)*(rp-rm)
				
				if(vn_b .lt. 0.d0) then !!*inflow
					s_b   = s_inf
					u_ref = u_inf
					v_ref = v_inf
					w_ref = w_inf
					turVal= nutinf
				end if
				
				if(vn_b .ge. 0.d0) then !!*outflow
					s_b   = s_p
					u_ref = u_p
					v_ref = v_p
					w_ref = w_p
					turVal=tur(i,j,ke-1)
				end if
				
				u_ref_t = u_ref*dxidx0 + v_ref*dxidy0 + w_ref*dxidz0
				
				u_b = u_ref + (vn_b - u_ref_t)*dxidx0
				v_b = v_ref + (vn_b - u_ref_t)*dxidy0
				w_b = w_ref + (vn_b - u_ref_t)*dxidz0
				d_b = (c_p*c_p/(gamma*s_b))**(1.d0/(gamma-1.d0))
				p_b = s_b*d_b**gamma  
				
				q(1,i,j,ke) = d_b
				q(2,i,j,ke) = d_b*u_b
				q(3,i,j,ke) = d_b*v_b
				q(4,i,j,ke) = d_b*w_b
				q(5,i,j,ke) = p_b/(gamma-1.d0) + 0.5d0*d_b*(u_b**2 + v_b**2 + w_b**2)
				if (iflag_turbulence .ge. 1) then
					tur(i,j,ke) = turVal
				end if
			end do
		end do	
	end if
	!!*end subface k+
		
	return
end subroutine