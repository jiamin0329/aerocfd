!!=========================================================
!! get jacobian matrix                                     
!!                                                         
!! author: jiamin xu                                       
!! date:   2012.12.25                                      
!!=========================================================
subroutine calc_jm(a,q1,q2,q3,q4,q5,dxidx,dxidy,dxidz,gamma)
	implicit none
	
	real*8  :: gamma
	real*8  :: q1,q2,q3,q4,q5
	real*8  :: rho,mu0,mv0,mw0,e0
	real*8  :: a(5,5)

	real*8  :: dxidx,dxidy,dxidz
	real*8  :: uu,phi,a1,a2,a3
	
	a2 = gamma-1.d0
	a3 = gamma-2.d0
	
	rho = q1
	mu0 = q2
	mv0 = q3
	mw0 = q4
	e0  = q5
	
	uu  = (dxidx*mu0 + dxidy*mv0 + dxidz*mw0)/rho
	phi = 0.5d0*a2*(mu0**2+mv0**2+mw0**2)/rho**2
	a1  = gamma*e0/rho - phi
	
	!!*
	a(1,1) = 0.d0
	a(1,2) = dxidx
	a(1,3) = dxidy
	a(1,4) = dxidz
	a(1,5) = 0.d0
	!!*
	a(2,1) = phi*dxidx - uu*mu0/rho
	a(2,2) = uu - a3*dxidx*mu0/rho
	a(2,3) = dxidy*mu0/rho - a2*dxidx*mv0/rho
	a(2,4) = dxidz*mu0/rho - a2*dxidx*mw0/rho
	a(2,5) = a2*dxidx
	!!*
	a(3,1) = phi*dxidy - uu*mv0/rho
	a(3,2) = dxidx*mv0/rho - a2*dxidy*mu0/rho
	a(3,3) = uu - a3*dxidy*mv0/rho
	a(3,4) = dxidz*mv0/rho - a2*dxidy*mw0/rho
	a(3,5) = a2*dxidy
	!!*
	a(4,1) = phi*dxidz - uu*mw0/rho
	a(4,2) = dxidx*mw0/rho - a2*dxidz*mu0/rho
	a(4,3) = dxidy*mw0/rho - a2*dxidz*mv0/rho
	a(4,4) = uu - a3*dxidz*mw0/rho
	a(4,5) = a2*dxidz
	!!*
	a(5,1) = (phi-a1)*uu
	a(5,2) = a1*dxidx - a2*uu*mu0/rho
	a(5,3) = a1*dxidy - a2*uu*mv0/rho
	a(5,4) = a1*dxidz - a2*uu*mw0/rho
	a(5,5) = gamma*uu
	!!*
	return
end subroutine