!!=========================================================!!
!! unsteady second order crank nicholson method            !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.12.24                                      !!
!!=========================================================!!
subroutine crank2nd(q,q0,dt,pri_v,rhs,rhsi,rhsv,dxidx,inv_j, &
	                is,ie,js,je,ks,ke, &
					is0,ie0,js0,je0,ks0,ke0, &
					is1,ie1,js1,je1,ks1,ke1, &
					gamma,mm,ssub)
    use flag_var
	implicit none
	!!*********************************************************************************
	integer :: mm, ssub
	integer :: i,j,k
	real*8  :: gamma
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: inv_j (       is:ie,js:je,ks:ke)
	real*8  :: dxidx (     9,is:ie,js:je,ks:ke)
	real*8  :: rhsi  (     5,is:ie,js:je,ks:ke)
	real*8  :: rhsv  (     5,is:ie,js:je,ks:ke)
	real*8  :: rhs   (     5,is:ie,js:je,ks:ke)
	real*8  :: pri_v (     7,is:ie,js:je,ks:ke)
	real*8  :: dt    (       is:ie,js:je,ks:ke)
	real*8  :: q0    (-1:0,5,is:ie,js:je,ks:ke)
	real*8  :: q     (     5,is:ie,js:je,ks:ke)
	!!*********************************************************************************
	integer :: mdim
	real*8  :: d,u,v,w,q2,p,c,h
	real*8  :: dxidx0,dxidy0,dxidz0,temp
	real*8  :: cu
	real*8  :: b1,b2,b3
	real*8  :: a1,a2,a3,a4,a5
	real*8  :: s1,s2,s3,s4,s5
	real*8  :: div,hdv,r
	real*8  :: alph, alph_i
	real*8,dimension(:,:,:),allocatable  :: dtsub
	real*8,allocatable,dimension(:,:,:) :: rou
	real*8,allocatable,dimension(:,:,:,:) :: dq
	character(len = 180) :: debugFile

	allocate (rou  (  is:ie,js:je,ks:ke))
	allocate (dtsub(  is:ie,js:je,ks:ke))
	allocate (dq   (5,is:ie,js:je,ks:ke))
	!!for viscous term in greatest eigenvalue
	!!=========================================================                   
	!! for innter iteration
	!! crank-nicolson method:             dtsub =      0.5d0*dt
	!! euler implicit method:             dtsub =      1.0d0*dt
	!! 3-points backward implicit method: dtsub = 2.d0/3.0d0*dt
	!!=========================================================
	!! for inner iteration parameters
	!! aphi2 = 0.0d0                                 !1st order
	!! aphi2 = 0.5d0                                 !2nd order
	!!=========================================================
	!!get right hand side terms
	!!rhs = rhsi + rhsv
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		rhs(:,i,j,k) = rhsi(:,i,j,k) + rhsv(:,i,j,k)
	end do
	end do
	end do
	!!*
	
	alph   = 0.5d0
	alph_i = 1.d0/(1.d0 + alph)
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		rhs(:,i,j,k) = alph_i*(rhs(:,i,j,k)*dt(i,j,k) - &
			          ((1.d0+     alph)*q (   :,i,j,k)  &
			          -(1.d0+2.d0*alph)*q0( 0,:,i,j,k)  &
			          +(          alph)*q0(-1,:,i,j,k))*inv_j(i,j,k))
	end do
	end do
	end do

	if (iflag_dimension .eq. iflag_2d) then
		k=1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			dtsub(i,j,k) = dt(i,j,k)/(1.d0+alph)
		end do
		end do
	else if (iflag_dimension .eq. iflag_3d) then
		do k = ks1-1,ke1+1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			dtsub(i,j,k) = dt(i,j,k)/(1.d0+alph)
		end do
		end do
		end do
	end if
	
	!!initialize dq
	dq = 0.d0
	
	if (iflag_dimension .eq. iflag_2d) then
		!!==============================================================
		!!                      xi direction                                                      
		!!==============================================================
		do k = ks1,ke1
		do j = js1,je1
			do i = is1-1,ie1+1
				d  = pri_v(1,i,j,k)
				u  = pri_v(2,i,j,k)
				v  = pri_v(3,i,j,k)
				w  = pri_v(4,i,j,k)
				p  = pri_v(5,i,j,k)
				c  = pri_v(7,i,j,k)
					
				dxidx0 = dxidx(1,i,j,k)
				dxidy0 = dxidx(2,i,j,k)
				dxidz0 = dxidx(3,i,j,k)
					
				temp = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp = sqrt(temp)
				cu   = u*dxidx0 + v*dxidy0 + w*dxidz0
					
				rou(i,j,k) = abs(cu) + abs(c*temp)
			end do  
			!!forward sweep
			do i = is1, ie1
				d  = pri_v(1,i-1,j,k)
				u  = pri_v(2,i-1,j,k)
				v  = pri_v(3,i-1,j,k)
				w  = pri_v(4,i-1,j,k)
				p  = pri_v(5,i-1,j,k)
				c  = pri_v(7,i-1,j,k)
					
				dxidx0 = dxidx(1,i-1,j,k)
				dxidy0 = dxidx(2,i-1,j,k)
				dxidz0 = dxidx(3,i-1,j,k)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2)         		
				r  = rou(i-1,j,k)
				h  = (c**2.d0)/(gamma-1.d0) + q2
				cu = dxidx0*u + dxidy0*v + dxidz0*w
					
				s1 = dq(1,i-1,j,k)
				s2 = dq(2,i-1,j,k)
				s3 = dq(3,i-1,j,k)
				s4 = dq(4,i-1,j,k)
				s5 = dq(5,i-1,j,k)
				!!========================================                                                  
				!!sub iteration
				hdv = 0.5d0*dtsub(i,j,k)
				b1  = hdv*(-cu*s1 + dxidx0*s2 + dxidy0*s3 + dxidz0*s4)
				b2  = hdv*( q2*s1      - u*s2      - v*s3      - w*s4 + s5)*(gamma-1.d0)
				b3  = hdv*( cu+r )
				a1  = b1               + b3*s1
				a2  = b1*u + b2*dxidx0 + b3*s2
				a3  = b1*v + b2*dxidy0 + b3*s3
				a4  = b1*w + b2*dxidz0 + b3*s4
				a5  = b1*h + b2*cu     + b3*s5
				!!========================================
				r = rou(i,j,k)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(rhs(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(rhs(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(rhs(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(rhs(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(rhs(5,i,j,k) + a5)
			end do
			!!backward sweep   
			do i = ie1,is1,-1 
				d  = pri_v(1,i+1,j,k)
				u  = pri_v(2,i+1,j,k)
				v  = pri_v(3,i+1,j,k)
				w  = pri_v(4,i+1,j,k)
				p  = pri_v(5,i+1,j,k)
				c  = pri_v(7,i+1,j,k)
					
				dxidx0 = dxidx(1,i+1,j,k)
				dxidy0 = dxidx(2,i+1,j,k)
				dxidz0 = dxidx(3,i+1,j,k)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2) 
				r  = rou(i+1,j,k) 
				h  = (c**2.d0)/(gamma-1.d0) + q2 
				cu = dxidx0*u + dxidy0*v + dxidz0*w 
					
				s1 = dq(1,i+1,j,k) 
				s2 = dq(2,i+1,j,k) 
				s3 = dq(3,i+1,j,k) 
				s4 = dq(4,i+1,j,k) 
				s5 = dq(5,i+1,j,k) 
				!!=================================================                                                   
				!!sub iteration
				hdv = 0.5d0*dtsub(i,j,k)
				b1  = hdv*(-cu*s1 + dxidx0*s2 + dxidy0*s3 + dxidz0*s4) 
				b2  = hdv*( q2*s1      - u*s2      - v*s3      - w*s4 + s5)*(gamma-1.d0)
				b3  = hdv*( cu-r ) 
				a1  = b1               + b3*s1
				a2  = b1*u + b2*dxidx0 + b3*s2
				a3  = b1*v + b2*dxidy0 + b3*s3
				a4  = b1*w + b2*dxidz0 + b3*s4
				a5  = b1*h + b2*cu     + b3*s5
				!!=================================================
				r = rou(i,j,k)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
		end do
		end do
		!!=====================end xi direction==========================
  	
		!!===============================================================
		!!                       eta direction                                                     
		!!===============================================================
		do k = ks1,ke1
		do i = is1,ie1
			do j = js1-1,je1+1
				d  = pri_v(1,i,j,k)
				u  = pri_v(2,i,j,k)
				v  = pri_v(3,i,j,k)
				w  = pri_v(4,i,j,k)
				p  = pri_v(5,i,j,k)
				c  = pri_v(7,i,j,k)
					
				dxidx0 = dxidx(4,i,j,k)
				dxidy0 = dxidx(5,i,j,k)
				dxidz0 = dxidx(6,i,j,k)
					
				temp = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp = sqrt(temp)
					
				cu   = u*dxidx0 + v*dxidy0 + w*dxidz0
				
				rou(i,j,k) = abs(cu) + abs(c*temp)
			end do  
			!! forward sweep  
			do j = js1, je1
				d  = pri_v(1,i,j-1,k)
				u  = pri_v(2,i,j-1,k)
				v  = pri_v(3,i,j-1,k)
				w  = pri_v(4,i,j-1,k)
				p  = pri_v(5,i,j-1,k)
				c  = pri_v(7,i,j-1,k)
					
				dxidx0 = dxidx(4,i,j-1,k)
				dxidy0 = dxidx(5,i,j-1,k)
				dxidz0 = dxidx(6,i,j-1,k)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2)         			
				r  = rou(i,j-1,k)
				h  = (c**2.d0)/(gamma-1.d0) + q2 
				cu = dxidx0*u + dxidy0*v + dxidz0*w  
					
				s1 = dq(1,i,j-1,k) 
				s2 = dq(2,i,j-1,k) 
				s3 = dq(3,i,j-1,k) 
				s4 = dq(4,i,j-1,k) 
				s5 = dq(5,i,j-1,k)
				!!==================================================
				!!sub iteration
				hdv = 0.5d0*dtsub(i,j,k)
				b1  = hdv*(-cu*s1 + dxidx0*s2 + dxidy0*s3 + dxidz0*s4) 
				b2  = hdv*( q2*s1      - u*s2      - v*s3      - w*s4 + s5)*(gamma-1.d0)
				b3  = hdv*( cu+r ) 
				a1  = b1               + b3*s1
				a2  = b1*u + b2*dxidx0 + b3*s2
				a3  = b1*v + b2*dxidy0 + b3*s3
				a4  = b1*w + b2*dxidz0 + b3*s4
				a5  = b1*h + b2*cu     + b3*s5
				!!==================================================
				r = rou(i,j,k)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(dq(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(dq(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(dq(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(dq(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(dq(5,i,j,k) + a5)
			end do  
			!! backward sweep   
			do j = je1,js1,-1
				d  = pri_v(1,i,j+1,k)
				u  = pri_v(2,i,j+1,k)
				v  = pri_v(3,i,j+1,k)
				w  = pri_v(4,i,j+1,k)
				p  = pri_v(5,i,j+1,k)
				c  = pri_v(7,i,j+1,k)
					
				dxidx0 = dxidx(4,i,j+1,k)
				dxidy0 = dxidx(5,i,j+1,k)
				dxidz0 = dxidx(6,i,j+1,k)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2)         				
				r  = rou(i,j+1,k) 
				h  = (c**2.d0)/(gamma-1.d0) + q2
					
				cu = dxidx0*u + dxidy0*v + dxidz0*w
					
				s1 = dq(1,i,j+1,k) 
				s2 = dq(2,i,j+1,k) 
				s3 = dq(3,i,j+1,k) 
				s4 = dq(4,i,j+1,k) 
				s5 = dq(5,i,j+1,k) 
				!!================================================                                                  
				!!sub iteration
				hdv = 0.5d0*dtsub(i,j,k)
				b1  = hdv*(-cu*s1 + dxidx0*s2 + dxidy0*s3 + dxidz0*s4) 
				b2  = hdv*( q2*s1      - u*s2      - v*s3      - w*s4 + s5)*(gamma-1.d0)
				b3  = hdv*( cu-r ) 
				a1  = b1               + b3*s1 
				a2  = b1*u + b2*dxidx0 + b3*s2 
				a3  = b1*v + b2*dxidy0 + b3*s3 
				a4  = b1*w + b2*dxidz0 + b3*s4 
				a5  = b1*h + b2*cu     + b3*s5 
				!!================================================  
				r = rou(i,j,k)                                  
				!!sub iteration           
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
		end do
		end do
		!!=====================end eta direction=========================	
	else if (iflag_dimension .eq. iflag_3d) then
		print *, "2nd crank for 3d is not supported!"
		stop
	end if
	
	!!update q
	!!q(n+1) = q(n) + dq(n)
	do k = ks1, ke1
	do j = js1, je1
	do i = is1, ie1
		q(:,i,j,k) = q(:,i,j,k) + dq(:,i,j,k)/inv_j(i,j,k)
	end do
	end do
	end do
	!!*end update q

	deallocate(rou,dtsub,dq)

	return
end subroutine crank2nd 