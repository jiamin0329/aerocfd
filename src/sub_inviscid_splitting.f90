subroutine roe_3d (ql,qr,flux,dxidx,inv_j,gamma,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,coord)
	implicit none
	
	integer :: coord
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: gamma
	real*8  :: dxidx(9,is:ie,js:je,ks:ke)
	real*8  :: inv_j(  is:ie,js:je,ks:ke)
	real*8  :: ql   (5,is:ie,js:je,ks:ke)
	real*8  :: qr   (5,is:ie,js:je,ks:ke)
	real*8  :: flux (5,is:ie,js:je,ks:ke)
	
	integer :: i,j,k
	integer :: ii,jj,kk
	integer :: n
	
	real*8  :: dxidx1,dxidy1,dxidz1
	real*8  :: dxidx0,dxidy0,dxidz0
	real*8  :: temp
	!!left&right status variable
	real*8  :: dl,ul,vl,wl,pl,hl
	real*8  :: dr,ur,vr,wr,pr,hr
	!!left&right contravariant velocity
	!!left&right unit contravariant velocity
	real*8  :: uul,uur,ubarl,ubarr
	!!left&right status flux
	real*8  :: fl(5),fr(5)
	!!roe average term
	real*8  :: tmpp
	real*8  :: avd,avu,avv,avw,avh,ava
	
	real*8  :: avuu
	real*8  :: tmp,tmp1,tmp2,harten
	real*8  :: lamda(5)
	!!difference of left and right status flux
	real*8  :: deltad,deltau,deltav,deltaw,deltap,deltauu
	!!
	real*8  :: alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8
	!!|a~|(qr-ql)
	real*8  :: f(5)
	real*8,parameter :: delt=0.06d0
	!!*
	
	if     (coord .eq. 1)then
		ii = 1
		jj = 0
		kk = 0
	else if(coord .eq. 2)then
		ii = 0
		jj = 1
		kk = 0
	else if(coord .eq. 3)then
		ii = 0
		jj = 0
		kk = 1
	end if
	!!*		
	
	!!here, i for i+0.5, j for j+0.5, k for k+0.5
	do k = ks0-kk,ke0
	do j = js0-jj,je0
	do i = is0-ii,ie0

		dxidx1 = 0.5d0*(dxidx(coord*3-2,i,j,k)*inv_j(i,j,k) + dxidx(coord*3-2,i+ii,j+jj,k+kk)*inv_j(i+ii,j+jj,k+kk))
		dxidy1 = 0.5d0*(dxidx(coord*3-1,i,j,k)*inv_j(i,j,k) + dxidx(coord*3-1,i+ii,j+jj,k+kk)*inv_j(i+ii,j+jj,k+kk))
		dxidz1 = 0.5d0*(dxidx(coord*3  ,i,j,k)*inv_j(i,j,k) + dxidx(coord*3  ,i+ii,j+jj,k+kk)*inv_j(i+ii,j+jj,k+kk))
				
		temp   = sqrt(dxidx1**2 + dxidy1**2 + dxidz1**2)
		
		dxidx0 = dxidx1/temp
		dxidy0 = dxidy1/temp
		dxidz0 = dxidz1/temp
				
		dl = ql(1,i,j,k); dr = qr(1,i,j,k)
		ul = ql(2,i,j,k); ur = qr(2,i,j,k)
		vl = ql(3,i,j,k); vr = qr(3,i,j,k)
		wl = ql(4,i,j,k); wr = qr(4,i,j,k)
		pl = ql(5,i,j,k); pr = qr(5,i,j,k)   
				
		hl = gamma*pl/((gamma-1.d0)*dl) + 0.5d0*(ul**2 + vl**2 + wl**2)
		hr = gamma*pr/((gamma-1.d0)*dr) + 0.5d0*(ur**2 + vr**2 + wr**2)
				
		uul = ul*dxidx1 + vl*dxidy1 + wl*dxidz1
		uur = ur*dxidx1 + vr*dxidy1 + wr*dxidz1
				
		!!left and right unit contravariant velocity
		ubarl = uul/temp
		ubarr = uur/temp
				
		!!left flux fl
		fl(1) = dl*uul
		fl(2) = dl*ul*uul + dxidx1*pl
		fl(3) = dl*vl*uul + dxidy1*pl
		fl(4) = dl*wl*uul + dxidz1*pl
		fl(5) = dl*hl*uul
		!!right flux fr
		fr(1) = dr*uur
		fr(2) = dr*ur*uur + dxidx1*pr
		fr(3) = dr*vr*uur + dxidy1*pr
		fr(4) = dr*wr*uur + dxidz1*pr
		fr(5) = dr*hr*uur	
		!!

		tmpp = sqrt(dr/dl)
				
		avd  = sqrt(dl*dr)
		avu  = (ul + ur*tmpp)/(1.d0+tmpp)
		avv  = (vl + vr*tmpp)/(1.d0+tmpp)
		avw  = (wl + wr*tmpp)/(1.d0+tmpp)
		avh  = (hl + hr*tmpp)/(1.d0+tmpp)
		ava  = sqrt((gamma-1.d0)*(avh - 0.5d0*(avu**2 + avv**2 + avw**2)))
				
		avuu = avu*dxidx0 + avv*dxidy0 + avw*dxidz0
				
		lamda(1) = abs(avuu - ava)
		lamda(2) = abs(avuu)
		lamda(3) = abs(avuu)
		lamda(4) = abs(avuu)
		lamda(5) = abs(avuu + ava)
				
		do n = 1,5
			tmp    = sign(1.d0, delt-lamda(n))
					
			tmp1   = 0.5d0*(1.d0 - tmp)
			tmp2   = 0.5d0*(1.d0 + tmp)
			harten = (lamda(n)**2 + delt**2)/(2.d0*delt)
					
			lamda(n) = tmp1*lamda(n) + tmp2*harten
		end do
				
		!!
		deltad  = dr    - dl
		deltap  = pr    - pl
		deltau  = ur    - ul
		deltav  = vr    - vl
		deltaw  = wr    - wl
		deltauu = ubarr - ubarl
				
		!!
		alpha1 = lamda(2)*(deltad - deltap/ava**2)
		alpha2 = lamda(5)*(deltap + avd*ava*deltauu)/(2.d0*ava**2)
		alpha3 = lamda(1)*(deltap - avd*ava*deltauu)/(2.d0*ava**2)
		alpha4 = alpha1 + alpha2 + alpha3
		alpha5 = ava*(alpha2 - alpha3)
		alpha6 = lamda(2)*(avd*deltau - avd*deltauu*dxidx0)
		alpha7 = lamda(2)*(avd*deltav - avd*deltauu*dxidy0)
		alpha8 = lamda(2)*(avd*deltaw - avd*deltauu*dxidz0)
			
		!!f = |a~|*(qr-ql) = |a~|*deltaq
		f(1) = alpha4
		f(2) = avu*alpha4 + alpha5*dxidx0 + alpha6
		f(3) = avv*alpha4 + alpha5*dxidy0 + alpha7
		f(4) = avw*alpha4 + alpha5*dxidz0 + alpha8
		f(5) = avh*alpha4 + avuu*alpha5 + avu*alpha6 + avv*alpha7 + avw*alpha8 - ava**2*alpha1/(gamma-1.d0)
				
		!!flux = 1/2*(f(qr)+f(ql)-|a~|(qr-ql))
		flux(1,i,j,k) = 0.5d0*(fl(1) + fr(1) - temp*f(1))
		flux(2,i,j,k) = 0.5d0*(fl(2) + fr(2) - temp*f(2))
		flux(3,i,j,k) = 0.5d0*(fl(3) + fr(3) - temp*f(3))
		flux(4,i,j,k) = 0.5d0*(fl(4) + fr(4) - temp*f(4))
		flux(5,i,j,k) = 0.5d0*(fl(5) + fr(5) - temp*f(5))
		!!if (coord .eq. 1) then
		!!print *, dl,dr,i,j
		!!end if
	end do
	end do
	end do
	!!*
	
	return
end subroutine
