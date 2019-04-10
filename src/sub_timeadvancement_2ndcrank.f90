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
					gamma)
    use flag_var
	implicit none
	!!*********************************************************************************
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
	real*8,allocatable,dimension(:,:) :: rou
	real*8,allocatable,dimension(:,:,:,:) :: dq
	mdim = max(ie-is+1,je-js+1,ke-ks+1)
	allocate (rou  (mdim,mdim))
	allocate (dtsub(is:ie,js:je,ks:ke))
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
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		rhs(:,i,j,k) = rhsi(:,i,j,k) + rhsv(:,i,j,k)
	end do
	end do
	end do
	!!*
	
	alph   = 0.5d0
	alph_i = 1.d0/(1.d0 + alph)
	
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		rhs(:,i,j,k) = alph_i*(rhs(:,i,j,k)*dt(i,j,k) - &
			          ((1.d0+     alph)*q (   :,i,j,k)  &
			          -(1.d0+2.d0*alph)*q0( 0,:,i,j,k)  &
			          +(          alph)*q0(-1,:,i,j,k))*inv_j(i,j,k))
	end do
	end do
	end do
	
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		dtsub(i,j,k) = dt(i,j,k)/(1.d0+alph)
	end do
	end do
	end do
	
	!!initialize dq
	dq = 0.d0
	
	if (iflag_dimension .eq. iflag_2d) then
		!!==============================================================
		!!                      xi direction                                                      
		!!==============================================================
		do k = ks0,ke0
			!! get vector r
			do j = js0,je0
			do i = is0,ie0
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
					
				rou(i,j) = abs(cu) + abs(c*temp)
			end do
			end do
			!! end get vector r
			!!============================================================   
			!!                     forward sweep
			!!============================================================
			do j = js0, je0
			do i = is0+1, ie0
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
				r  = rou(i-1,j)
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(rhs(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(rhs(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(rhs(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(rhs(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(rhs(5,i,j,k) + a5)
			end do
			end do
			!!==================end forward sweep=======================
  	
			!!============================================================   
			!!                    backward sweep
			!!============================================================    
			do j = js0,je0
			do i = ie0-1,is0,-1 
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
				r  = rou(i+1,j) 
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
			end do
			!!==================end backward sweep=======================
		end do
		!!=====================end xi direction==========================
  	
		!!===============================================================
		!!                       eta direction                                                     
		!!===============================================================
		do k = ks0,ke0
			!! get vector r
			do j = js0,je0
			do i = is0,ie0
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
				
				rou(i,j) = abs(cu) + abs(c*temp)
			end do
			end do
			!!end get vector r
			!!============================================================   
			!!                     forward sweep
			!!============================================================     
			do j = js0+1, je0
			do i = is0, ie0
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
				r  = rou(i,j-1)
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(dq(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(dq(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(dq(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(dq(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(dq(5,i,j,k) + a5)
			end do
			end do
			!!==================end forward sweep=======================
			
			!!==========================================================   
			!!                    backward sweep
			!!==========================================================    
			do j = je0-1,js0,-1
			do i = is0,ie0
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
				r  = rou(i,j+1) 
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
				r = rou(i,j)                                  
				!!sub iteration           
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
			end do
			!!==================end backward sweep=======================
			!!*
		end do
		!!=====================end eta direction=========================	
	else if (iflag_dimension .eq. iflag_3d) then
		!!==============================================================
		!!                      xi direction                                                      
		!!==============================================================
		do k = ks0,ke0
			!! get vector r
			do j = js0,je0
			do i = is0,ie0
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
					
				rou(i,j) = abs(cu) + abs(c*temp)
			end do
			end do
			!! end get vector r
			!!============================================================   
			!!                     forward sweep
			!!============================================================
			do j = js0, je0
			do i = is0+1, ie0
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
				r  = rou(i-1,j)
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(rhs(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(rhs(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(rhs(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(rhs(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(rhs(5,i,j,k) + a5)
			end do
			end do
			!!==================end forward sweep=======================
  	
			!!============================================================   
			!!                    backward sweep
			!!============================================================    
			do j = js0,je0
			do i = ie0-1,is0,-1 
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
				r  = rou(i+1,j) 
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
			end do
			!!==================end backward sweep=======================
		end do
		!!=====================end xi direction==========================
  	
		!!===============================================================
		!!                       eta direction                                                     
		!!===============================================================
		do k = ks0,ke0
			!! get vector r
			do j = js0,je0
			do i = is0,ie0
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
					
				rou(i,j) = abs(cu) + abs(c*temp)
			end do
			end do
			!!end get vector r
			!!============================================================   
			!!                     forward sweep
			!!============================================================     
			do j = js0+1, je0
			do i = is0, ie0
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
				r  = rou(i,j-1)
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
				r = rou(i,j)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = div*(dq(1,i,j,k) + a1)
				dq(2,i,j,k) = div*(dq(2,i,j,k) + a2)
				dq(3,i,j,k) = div*(dq(3,i,j,k) + a3)
				dq(4,i,j,k) = div*(dq(4,i,j,k) + a4)
				dq(5,i,j,k) = div*(dq(5,i,j,k) + a5)
			end do
			end do
			!!==================end forward sweep=======================
			
			!!==========================================================   
			!!                    backward sweep
			!!==========================================================    
			do j = je0-1,js0,-1
			do i = is0,ie0
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
				r  = rou(i,j+1) 
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
				r = rou(i,j)                                  
				!!sub iteration           
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r)
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
			end do
			!!==================end backward sweep=======================
			!!*
		end do
		!!=====================end eta direction=========================
		
		!!===============================================================
		!!                      zeta direction
		!!===============================================================
		do j = js0,je0
			!!get vector r
			do k = ks0,ke0
			do i = is0,ie0
				d  = pri_v(1,i,j,k)
				u  = pri_v(2,i,j,k)
				v  = pri_v(3,i,j,k)
				w  = pri_v(4,i,j,k)
				p  = pri_v(5,i,j,k)
				c  = pri_v(7,i,j,k)
					
				dxidx0 = dxidx(7,i,j,k)
				dxidy0 = dxidx(8,i,j,k)
				dxidz0 = dxidx(9,i,j,k)
  	              
				temp = dxidx0**2 + dxidy0**2 + dxidz0**2
				temp = sqrt(temp)
				cu   = u*dxidx0 + v*dxidy0 + w*dxidz0
					
				rou(i,k) = abs(cu) + abs(c*temp)
			end do
			end do
			!!============================================================   
			!!                     forward sweep                                             
			!!============================================================     
			do k = ks0+1,ke0
			do i = is0,ie0
				d  = pri_v(1,i,j,k-1)
				u  = pri_v(2,i,j,k-1)
				v  = pri_v(3,i,j,k-1)
				w  = pri_v(4,i,j,k-1)
				p  = pri_v(5,i,j,k-1)
				c  = pri_v(7,i,j,k-1)
					
				dxidx0 = dxidx(7,i,j,k-1)
				dxidy0 = dxidx(8,i,j,k-1)
				dxidz0 = dxidx(9,i,j,k-1)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2)         				
				r  = rou(i,k-1) 
				h  = (c**2.d0)/(gamma-1.d0) + q2 
				cu = dxidx0*u + dxidy0*v + dxidz0*w 
					
				s1 = dq(1,i,j,k-1) 
				s2 = dq(2,i,j,k-1) 
				s3 = dq(3,i,j,k-1) 
				s4 = dq(4,i,j,k-1) 
				s5 = dq(5,i,j,k-1)
				!!================================================                                                  
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
				!!================================================
				r = rou(i,k)
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r) 
				dq(1,i,j,k) = div*(dq(1,i,j,k) + a1) 
				dq(2,i,j,k) = div*(dq(2,i,j,k) + a2) 
				dq(3,i,j,k) = div*(dq(3,i,j,k) + a3) 
				dq(4,i,j,k) = div*(dq(4,i,j,k) + a4) 
				dq(5,i,j,k) = div*(dq(5,i,j,k) + a5) 
			end do
			end do
			!!==================end forward sweep=====================
			
			!!======================================================== 
			!!                   backward sweep
			!!========================================================   
			do k = ke0-1, ks0, -1
			do i = is0, ie0
				d  = pri_v(1,i,j,k+1)
				u  = pri_v(2,i,j,k+1)
				v  = pri_v(3,i,j,k+1)
				w  = pri_v(4,i,j,k+1)
				p  = pri_v(5,i,j,k+1)
				c  = pri_v(7,i,j,k+1)
					
				dxidx0 = dxidx(7,i,j,k+1)
				dxidy0 = dxidx(8,i,j,k+1)
				dxidz0 = dxidx(9,i,j,k+1)
					
				q2 = 0.5d0*(u**2 + v**2 + w**2)        				
				r  = rou(i,k+1)
				h  = (c**2.d0)/(gamma-1.d0) + q2
					
				cu = dxidx0*u + dxidy0*v + dxidz0*w
					
				s1 = dq(1,i,j,k+1)
				s2 = dq(2,i,j,k+1)
				s3 = dq(3,i,j,k+1)
				s4 = dq(4,i,j,k+1)
				s5 = dq(5,i,j,k+1)
				!!================================================                                                 
				!!sub iteration
				hdv = 0.5d0*dtsub(i,j,k)
				b1  = hdv*(- cu*s1 + dxidx0*s2 + dxidy0*s3 + dxidz0*s4)
				b2  = hdv*(  q2*s1      - u*s2      - v*s3      - w*s4 + s5)*(gamma-1.d0)
				b3  = hdv*(  cu-r )
				a1  = b1               + b3*s1
				a2  = b1*u + b2*dxidx0 + b3*s2
				a3  = b1*v + b2*dxidy0 + b3*s3
				a4  = b1*w + b2*dxidz0 + b3*s4
				a5  = b1*h + b2*cu     + b3*s5
				!!================================================
				r = rou(i,k)                                  
				!!sub iteration
				div = 1.d0/(1.d0 + dtsub(i,j,k)*r) 
				dq(1,i,j,k) = dq(1,i,j,k) - div*a1
				dq(2,i,j,k) = dq(2,i,j,k) - div*a2
				dq(3,i,j,k) = dq(3,i,j,k) - div*a3
				dq(4,i,j,k) = dq(4,i,j,k) - div*a4
				dq(5,i,j,k) = dq(5,i,j,k) - div*a5
			end do
			end do
			!!==================end backward sweep=======================
		end do
		!!====================end zeta direction=======================	
	end if
	
	!!update q
	!!q(n+1) = q(n) + dq(n)
	do k = ks0, ke0
	do j = js0, je0
	do i = is0, ie0
		q(:,i,j,k) = q(:,i,j,k) + dq(:,i,j,k)/inv_j(i,j,k)
	end do
	end do
	end do
	!!*end update q

	deallocate(rou,dtsub,dq)

	return
end subroutine crank2nd 