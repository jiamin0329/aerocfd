!!*********************************************************
!!	get local time-step                                    
!!                                                         
!!	ref: CFL3D User's Manual Page.236                      
!!	                                                       
!!  purpose: to get local time-step scale                  
!!���������� for steady simulation                         
!!                                                         
!!  input:                                                 
!!	block dimension:	is,ie,js,je,ks,ke                    
!!	reynolds number:	re                                   
!!	clf number:       cfl                                  
!!	jacobian matrix:  dxidx,dxidx2                         
!!	viscosity coeff:  amu                                  
!!	primitive variable: pri_v                              
!!                                                         
!!	output:                                                
!!	local time-step: dt                                    
!!                                                         
!!                                                         
!!	Author: Jiamin Xu                                      
!!	Date:   2012.12.24                                     
!!*********************************************************
subroutine find_dt(dt,dtmax,dtmin,cflmax,cflmin,pri_v,amu,dxidx, &
                   re,cfl,is,ie,js,je,ks,ke,iflag_time,iflag_steady,iflag_unsteady)
	implicit none
	!!*******************************************************   
	integer :: iflag_time,iflag_steady,iflag_unsteady
	integer :: is,ie,js,je,ks,ke
	real*8  :: cfl,re
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	real*8  :: dxidx2(3,is:ie,js:je,ks:ke)
	real*8  :: amu   (3,is:ie,js:je,ks:ke)
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: dtmax,dtmin,cflmax,cflmin
	real*8  :: dt    (  is:ie,js:je,ks:ke)
	!!*******************************************************  
	integer :: i,j,k 
	real*8  :: dxidx0,  dxidy0,  dxidz0
	real*8  :: detadx0, detady0, detadz0
	real*8  :: dzetadx0,dzetady0,dzetadz0
	real*8  :: tempxi,tempeta,tempzeta
	real*8  :: d0,u0,v0,w0,p0,c0
	real*8  :: uu,vv,ww
	real*8  :: t1,t2,t3
	!!*******************************************************
	dtmin  = 1.d10
	dtmax  = 0.d0
	cflmin = 1.d10
	cflmax = 0.d0
	!!*******************************************************
	if (iflag_time .eq. iflag_steady)then
		do k = ks,ke
		do j = js,je
		do i = is,ie
			dxidx0   = dxidx(1,i,j,k); dxidy0   = dxidx(2,i,j,k); dxidz0   = dxidx(3,i,j,k)
			detadx0  = dxidx(4,i,j,k); detady0  = dxidx(5,i,j,k); detadz0  = dxidx(6,i,j,k)
			dzetadx0 = dxidx(7,i,j,k); dzetady0 = dxidx(8,i,j,k); dzetadz0 = dxidx(9,i,j,k)
				
			tempxi   = sqrt(dxidx0  **2 + dxidy0  **2 + dxidz0  **2)
			tempeta  = sqrt(detadx0 **2 + detady0 **2 + detadz0 **2)
			tempzeta = sqrt(dzetadx0**2 + dzetady0**2 + dzetadz0**2)
				
			d0  = pri_v(1,i,j,k)
			u0  = pri_v(2,i,j,k)
			v0  = pri_v(3,i,j,k)
			w0  = pri_v(4,i,j,k)
			p0  = pri_v(5,i,j,k)
			c0  = pri_v(7,i,j,k)
				
			uu = (dxidx0*u0   + dxidy0*v0   + dxidz0*w0  )/tempxi
			vv = (detadx0*u0  + detady0*v0  + detadz0*w0 )/tempeta
			ww = (dzetadx0*u0 + dzetady0*v0 + dzetadz0*w0)/tempzeta
				
			t1 = abs(uu) + c0 + 2.d0*tempxi  *amu(3,i,j,k)*4.d0/3.d0/re/d0
			t2 = abs(vv) + c0 + 2.d0*tempeta *amu(3,i,j,k)*4.d0/3.d0/re/d0
			t3 = abs(ww) + c0 + 2.d0*tempzeta*amu(3,i,j,k)*4.d0/3.d0/re/d0
				
			dt(i,j,k) = cfl/(tempxi*t1 + tempeta*t2 + tempzeta*t3) 
				
			if     (dt(i,j,k) .ge. dtmax)then
				dtmax = dt(i,j,k)
			else if(dt(i,j,k) .le. dtmin)then
				dtmin = dt(i,j,k)
			end if
		end do
		end do
		end do
	else if (iflag_time .eq. iflag_unsteady)then
		do k = ks,ke
		do j = js,je
		do i = is,ie
			dxidx0   = dxidx(1,i,j,k); dxidy0   = dxidx(2,i,j,k); dxidz0   = dxidx(3,i,j,k)
			detadx0  = dxidx(4,i,j,k); detady0  = dxidx(5,i,j,k); detadz0  = dxidx(6,i,j,k)
			dzetadx0 = dxidx(7,i,j,k); dzetady0 = dxidx(8,i,j,k); dzetadz0 = dxidx(9,i,j,k)
				
			tempxi   = sqrt(dxidx0  **2 + dxidy0  **2 + dxidz0  **2)
			tempeta  = sqrt(detadx0 **2 + detady0 **2 + detadz0 **2)
			tempzeta = sqrt(dzetadx0**2 + dzetady0**2 + dzetadz0**2)
				
			d0  = pri_v(1,i,j,k)
			u0  = pri_v(2,i,j,k)
			v0  = pri_v(3,i,j,k)
			w0  = pri_v(4,i,j,k)
			p0  = pri_v(5,i,j,k)
			c0  = pri_v(7,i,j,k)
				
			uu = (dxidx0*u0   + dxidy0*v0   + dxidz0*w0  )/tempxi
			vv = (detadx0*u0  + detady0*v0  + detadz0*w0 )/tempeta
			ww = (dzetadx0*u0 + dzetady0*v0 + dzetadz0*w0)/tempzeta
				
			t1 = abs(uu) + c0 + 2.d0*tempxi  *amu(3,i,j,k)*4.d0/3.d0/re/d0
			t2 = abs(vv) + c0 + 2.d0*tempeta *amu(3,i,j,k)*4.d0/3.d0/re/d0
			t3 = abs(ww) + c0 + 2.d0*tempzeta*amu(3,i,j,k)*4.d0/3.d0/re/d0
				
			cfl = dt(i,j,k)*(tempxi*t1 + tempeta*t2 + tempzeta*t3) 
				
			if     (cfl .ge. cflmax)then
				cflmax = cfl
			else if(cfl .le. cflmin)then
				cflmin = cfl
			end if
		end do
		end do
		end do
	end if
	!!*
	
	return
end subroutine
