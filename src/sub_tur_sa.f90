!!=========================================================!!
!! spalart-allmaras turbulence model                       !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2013.01.19                                      !!
!!=========================================================!!
subroutine sa_model(nut,sa_rhs,dt,pri_v,amu,dudx,vor,dxidx,inv_j,alpha,beta,gamma, &
					length,dist,spac,re,cfl, &
					is, ie, js, je, ks, ke,  &
					is0,ie0,js0,je0,ks0,ke0, &
					is1,ie1,js1,je1,ks1,ke1)
	use flag_var
	use sa_var
	implicit none
	!!***************************************************
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
  	real*8  :: re,cfl
  	real*8  :: dist  (  is:ie,js:je,ks:ke  )
  	real*8  :: spac  (  is:ie,js:je,ks:ke  )
  	real*8  :: length(  is:ie,js:je,ks:ke  )
  	real*8  :: dxidx (3,is:ie,js:je,ks:ke,3)
  	real*8  :: dxidx2(3,is:ie,js:je,ks:ke  )
  	real*8  :: inv_j (  is:ie,js:je,ks:ke  )
  	real*8  :: alpha (3,is:ie,js:je,ks:ke  )
  	real*8  :: beta  (3,is:ie,js:je,ks:ke  )
  	real*8  :: gamma (3,is:ie,js:je,ks:ke  )
  	real*8  :: dudx  (3,is:ie,js:je,ks:ke,3)
  	real*8  :: vor   (  is:ie,js:je,ks:ke  )
 	real*8  :: pri_v (7,is:ie,js:je,ks:ke  )
  	real*8  :: dt    (  is:ie,js:je,ks:ke  )
  	real*8  :: amu   (3,is:ie,js:je,ks:ke  )
  	real*8  :: sa_rhs(  is:ie,js:je,ks:ke  )
  	real*8  :: nut   (  is:ie,js:je,ks:ke  )
  	!!***************************************************
	real*8  :: cdes = 0.65d0
	real*8  :: dist0,spac0
	
	real*8  :: ux,uy,uz
	real*8  :: vx,vy,vz
	real*8  :: wx,wy,wz
	real*8  :: d0,nut0,nu0,temp,rd,fd
	integer :: i,j,k
	real*8,allocatable,dimension(:,:,:) :: diasd
	allocate(diasd(is:ie,js:je,ks:ke))
	
	sa_rhs = 0.d0
	
	!!set turbulent length
	do k = ks,ke
	do j = js,je
	do i = is,ie
		length(i,j,k) = dist(i,j,k)
	end do
	end do
	end do
	!!*
	
	!!des option
	!!des97
	if(iflag_des .eq. iflag_des97)then
		do k = ks,ke
		do j = js,je
		do i = is,ie
			spac0 = spac(i,j,k)
			dist0 = dist(i,j,k)
			length(i,j,k) = min(dist0, cdes*spac0)
		end do
		end do
		end do
	!!ddes
	else if(iflag_des .eq. iflag_ddes)then
		do k = ks,ke
		do j = js,je
		do i = is,ie
			!!velocity gradient
			ux = dudx(1,i,j,k,1); uy = dudx(1,i,j,k,2); uz = dudx(1,i,j,k,3)
			vx = dudx(2,i,j,k,1); vy = dudx(2,i,j,k,2); vz = dudx(2,i,j,k,3)
			wx = dudx(3,i,j,k,1); wy = dudx(3,i,j,k,2); wz = dudx(3,i,j,k,3)
			!!laminar eddy and turbulent eddy
			d0    = pri_v  (1,i,j,k)
			nut0  = nut    (  i,j,k)
			nu0   = amu    (1,i,j,k)/d0
			!!
			dist0 = dist(i,j,k)
			spac0 = spac(i,j,k)
					
			temp  = uy**2 + uz**2 &
				  + vx**2 + vz**2 &
				  + wx**2 + wy**2
			temp  = sqrt(temp)
	  				
			rd = (nu0 + nut0)/(temp*re*(kar**2)*(dist0**2))
			fd = 1.d0 - tanh((8.d0*rd)**3)
			length(i,j,k) = dist0 - fd*max(0.d0,dist0-cdes*spac0)
		end do
		end do
		end do
	end if
	!!*
    
	!!convective term
	call sa_rhs1 (sa_rhs,nut,pri_v,dxidx, &
	              is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1)
	!!dissipation term
	call sa_rhs2 (sa_rhs,nut,amu,pri_v,dxidx,inv_j,alpha,beta,gamma,re, &
	              is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1)
	!!production & destruction term
	call sa_rhs3 (sa_rhs,nut,pri_v,amu,vor,length,diasd,re, &
	              is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1) 
	!!*

	!!only 1st order time marching scheme is implemented
	!!in the turbulence equation for both steady and unsteady situation
	call sa_lusgs(nut,dt,sa_rhs,pri_v,dxidx,diasd,re,cfl,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1)
	!!*
    
	deallocate(diasd)

	return
end subroutine

!!=========================================================!!
!! convective terms of sa equation                         !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2013.01.27                                      !!
!!=========================================================!!
subroutine sa_rhs1(sa_rhs,nut,pri_v,dxidx,  &
				   is, ie, js, je, ks, ke,  &
				   is0,ie0,js0,je0,ks0,ke0, &
				   is1,ie1,js1,je1,ks1,ke1)
	use flag_var
	implicit none
	!!*********************************************************
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	real*8  :: nut   (  is:ie,js:je,ks:ke)
	real*8  :: sa_rhs(  is:ie,js:je,ks:ke)
	!!*********************************************************
	integer :: i,j,k
	real*8  :: u,v,w
	real*8  :: dxidx0,dxidy0,dxidz0
	real*8  :: uu,vv,ww
	real*8  :: nutm,nutc,nutp
	real*8  :: uxi,veta,wzeta
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		u = pri_v(2,i,j,k)
		v = pri_v(3,i,j,k)
		w = pri_v(4,i,j,k)
		!!*********************
		!!xi direction
		dxidx0 = dxidx(1,i,j,k)
		dxidy0 = dxidx(2,i,j,k)
		dxidz0 = dxidx(3,i,j,k)
		
		uu = dxidx0*u + dxidy0*v + dxidz0*w
		
		nutp = nut(i+1,j,k)
		nutc = nut(i,  j,k)
		nutm = nut(i-1,j,k)
		
		uxi = -(0.5d0*(uu + abs(uu))*(nutc - nutm) &
		      + 0.5d0*(uu - abs(uu))*(nutp - nutc))
		!!eta direction
		dxidx0 = dxidx(4,i,j,k)
		dxidy0 = dxidx(5,i,j,k)
		dxidz0 = dxidx(6,i,j,k)
		
		vv = dxidx0*u + dxidy0*v + dxidz0*w
		
		nutp = nut(i,j+1,k)
		nutc = nut(i,j  ,k)
		nutm = nut(i,j-1,k)
			
		veta = -(0.5d0*(vv + abs(vv))*(nutc - nutm) &
			   + 0.5d0*(vv - abs(vv))*(nutp - nutc)) 
			   
		if (iflag_dimension .eq. iflag_3d) then
			!!zeta direction
			dxidx0 = dxidx(7,i,j,k)
			dxidy0 = dxidx(8,i,j,k)
			dxidz0 = dxidx(9,i,j,k)
			
			ww = dxidx0*u + dxidy0*v + dxidz0*w
			
			nutp = nut(i,j,k+1)
			nutc = nut(i,j,k  )
			nutm = nut(i,j,k-1)
			
			wzeta = -(0.5d0*(ww + abs(ww))*(nutc - nutm) &
			    	+ 0.5d0*(ww - abs(ww))*(nutp - nutc))
			!!*
		else
			wzeta = 0.d0
		end if

		sa_rhs(i,j,k) = uxi+veta+wzeta
	end do
	end do
	end do
	!!* 		

	return
end subroutine sa_rhs1

!!=========================================================!!
!! diffusion terms of sa equation                          !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2013.01.27                                      !!
!!=========================================================!!
subroutine sa_rhs2(sa_rhs,nut,amu,pri_v,dxidx,inv_j,alpha,beta,gamma,re, &
				   is, ie, js, je, ks, ke,  &
	               is0,ie0,js0,je0,ks0,ke0, &
	               is1,ie1,js1,je1,ks1,ke1)
	use sa_var
	use flag_var
	implicit none
	!!***************************************************************
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: re
	real*8  :: alpha (3,is:ie,js:je,ks:ke)
	real*8  :: beta  (3,is:ie,js:je,ks:ke)
	real*8  :: gamma (3,is:ie,js:je,ks:ke)
	real*8  :: inv_j (  is:ie,js:je,ks:ke)
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: amu   (3,is:ie,js:je,ks:ke)
	real*8  :: nut   (  is:ie,js:je,ks:ke)
	real*8  :: sa_rhs(  is:ie,js:je,ks:ke)
	!!***************************************************************
	integer :: i,j,k
	
	real*8  :: nut0,nu0
	real*8  :: dxidx0,  dxidy0,  dxidz0
	real*8  :: detadx0, detady0, detadz0
	real*8  :: dzetadx0,dzetady0,dzetadz0
	real*8  :: inv_j0
	real*8  :: dnutdx,dnutdy,dnutdz
	real*8  :: deltanut
	real*8,allocatable,dimension(:,:,:) :: dnutdxi,  dnutdeta,dnutdzeta
	real*8,allocatable,dimension(:,:,:) :: tempalpha,tempbeta,tempgamma
	real*8,allocatable,dimension(:,:,:) :: tempa,tempb,tempc
	real*8,allocatable,dimension(:,:,:) :: temp
	real*8,allocatable,dimension(:,:,:) :: tempxi,tempeta,tempzeta
	real*8,allocatable,dimension(:,:,:) :: tempx, tempy,  tempz
	real*8  :: diff
		
	allocate(dnutdxi   (is:ie,js:je,ks:ke))
	allocate(dnutdeta  (is:ie,js:je,ks:ke))
	allocate(dnutdzeta (is:ie,js:je,ks:ke))
	allocate(tempalpha (is:ie,js:je,ks:ke))
	allocate(tempbeta  (is:ie,js:je,ks:ke))
	allocate(tempgamma (is:ie,js:je,ks:ke))
	allocate(tempa     (is:ie,js:je,ks:ke))
	allocate(tempb     (is:ie,js:je,ks:ke))
	allocate(tempc     (is:ie,js:je,ks:ke))
	allocate(temp      (is:ie,js:je,ks:ke))
	allocate(tempxi    (is:ie,js:je,ks:ke))
	allocate(tempeta   (is:ie,js:je,ks:ke))
	allocate(tempzeta  (is:ie,js:je,ks:ke))
	allocate(tempx     (is:ie,js:je,ks:ke))
	allocate(tempy     (is:ie,js:je,ks:ke))
	allocate(tempz     (is:ie,js:je,ks:ke))

	dnutdxi   = 0.d0
	dnutdeta  = 0.d0
	dnutdzeta = 0.d0
	tempalpha = 0.d0
	tempbeta  = 0.d0
	tempgamma = 0.d0
	tempa     = 0.d0
	tempb     = 0.d0
	tempc     = 0.d0
	temp      = 0.d0
	tempxi    = 0.d0
	tempeta   = 0.d0
	tempzeta  = 0.d0
	tempx     = 0.d0
	tempy     = 0.d0
	tempz     = 0.d0
	
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (i .eq. is0) then
			dnutdxi(i,j,k) = nut(i+1,j,k) - nut(i,j,k)
		else if(i .eq. ie0) then
			dnutdxi(i,j,k) = nut(i,j,k) - nut(i-1,j,k)
		else
			dnutdxi(i,j,k) = (nut(i+1,j,k) - nut(i-1,j,k))/2.d0
		end if    
	end do
	end do
	end do	
		
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (j .eq. js0) then
			dnutdeta(i,j,k) = nut(i,j+1,k) - nut(i,j,k)
		else if(j .eq. je0) then
			dnutdeta(i,j,k) = nut(i,j,k) - nut(i,j-1,k)
		else
			dnutdeta(i,j,k) = (nut(i,j+1,k) - nut(i,j-1,k))/2.d0
		end if    
	end do
	end do
	end do

	if      (iflag_dimension .eq. iflag_2d) then
		dnutdzeta = 0.d0
	else if (iflag_dimension .eq. iflag_3d) then
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			if     (k .eq. ks0) then
				dnutdzeta(i,j,k) = nut(i,j,k+1) - nut(i,j,k)
			else if(k .eq. ke0) then
				dnutdzeta(i,j,k) = nut(i,j,k) - nut(i,j,k-1)
			else
				dnutdzeta(i,j,k) = (nut(i,j,k+1) - nut(i,j,k-1))/2.d0
			end if    
		end do
		end do
		end do	
	end if
    
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		tempalpha(i,j,k) = alpha(1,i,j,k)*dnutdxi  (i,j,k) &
				         + alpha(2,i,j,k)*dnutdeta (i,j,k) &
				         + alpha(3,i,j,k)*dnutdzeta(i,j,k)
				
		tempbeta (i,j,k) = beta (1,i,j,k)*dnutdxi  (i,j,k) &
				         + beta (2,i,j,k)*dnutdeta (i,j,k) &
				         + beta (3,i,j,k)*dnutdzeta(i,j,k)
				                 
		tempgamma(i,j,k) = gamma(1,i,j,k)*dnutdxi  (i,j,k) &
				         + gamma(2,i,j,k)*dnutdeta (i,j,k) &
				         + gamma(3,i,j,k)*dnutdzeta(i,j,k)		
	end do	                   
	end do
	end do
	
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (i .eq. is0) then
			tempa(i,j,k) = tempalpha(i+1,j,k) - tempalpha(i,j,k)
		else if(i .eq. ie0) then
			tempa(i,j,k) = tempalpha(i,j,k) - tempalpha(i-1,j,k)
		else
			tempa(i,j,k) = (tempalpha(i+1,j,k) - tempalpha(i-1,j,k))/2.d0
		end if    
	end do
	end do
	end do	
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (j .eq. js0) then
			tempb(i,j,k) = tempbeta(i,j+1,k) - tempbeta(i,j,k)
		else if(j .eq. je0) then
			tempb(i,j,k) = tempbeta(i,j,k) - tempbeta(i,j-1,k)
		else
			tempb(i,j,k) = (tempbeta(i,j+1,k) - tempbeta(i,j-1,k))/2.d0
		end if    
	end do
	end do
	end do	

	if      (iflag_dimension .eq. iflag_2d)then
		tempc = 0.d0
	else if (iflag_dimension .eq. iflag_3d)then
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			if     (k .eq. ks0) then
				tempc(i,j,k) = tempgamma(i,j,k+1) - tempgamma(i,j,k)
			else if(k .eq. ke0) then
				tempc(i,j,k) = tempgamma(i,j,k) - tempgamma(i,j,k-1)
			else
				tempc(i,j,k) = (tempgamma(i,j,k+1) - tempgamma(i,j,k-1))/2.d0
			end if    
		end do
		end do
		end do
	end if
    
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		nut0 = nut(i,j,k)
		nu0  = amu(1,i,j,k)/pri_v(1,i,j,k)
		temp(i,j,k) = nu0 + (1.d0+cb2)*nut0
	end do
	end do
	end do
	
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (i .eq. is0) then
			tempxi(i,j,k) = temp(i+1,j,k) - temp(i,j,k)
		else if(i .eq. ie0) then
			tempxi(i,j,k) = temp(i,j,k) - temp(i-1,j,k)
		else
			tempxi(i,j,k) = (temp(i+1,j,k) - temp(i-1,j,k))/2.d0
		end if    
	end do
	end do
	end do	

	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		if     (j .eq. js0) then
			tempeta(i,j,k) = temp(i,j+1,k) - temp(i,j,k)
		else if(j .eq. je0) then
			tempeta(i,j,k) = temp(i,j,k) - temp(i,j-1,k)
		else
			tempeta(i,j,k) = (temp(i,j+1,k) - temp(i,j-1,k))/2.d0
		end if    
	end do
	end do
	end do	

	if      (iflag_dimension .eq. iflag_2d)then
		tempzeta = 0.d0
	else if (iflag_dimension .eq. iflag_3d)then
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			if     (k .eq. ks0) then
				tempzeta(i,j,k) = temp(i,j,k+1) - temp(i,j,k)
			else if(k .eq. ke0) then
				tempzeta(i,j,k) = temp(i,j,k) - temp(i,j,k-1)
			else
				tempzeta(i,j,k) = (temp(i,j,k+1) - temp(i,j,k-1))/2.d0
			end if    
		end do
		end do
		end do	
	end if   
     
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		dxidx0   = dxidx(1,i,j,k)
		dxidy0   = dxidx(2,i,j,k)
		dxidz0   = dxidx(3,i,j,k)
		detadx0  = dxidx(4,i,j,k)
		detady0  = dxidx(5,i,j,k)
		detadz0  = dxidx(6,i,j,k)
		dzetadx0 = dxidx(7,i,j,k)
		dzetady0 = dxidx(8,i,j,k)
		dzetadz0 = dxidx(9,i,j,k)				
				
		tempx(i,j,k) = tempxi(i,j,k)*dxidx0 + tempeta(i,j,k)*detadx0 + tempzeta(i,j,k)*dzetadx0
		tempy(i,j,k) = tempxi(i,j,k)*dxidy0 + tempeta(i,j,k)*detady0 + tempzeta(i,j,k)*dzetady0
		tempz(i,j,k) = tempxi(i,j,k)*dxidz0 + tempeta(i,j,k)*detadz0 + tempzeta(i,j,k)*dzetadz0
	end do
	end do
	end do 
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		nut0 = nut(i,j,k)
		nu0  = amu(1,i,j,k)/pri_v(1,i,j,k)
				
		dxidx0   = dxidx(1,i,j,k)
		dxidy0   = dxidx(2,i,j,k)
		dxidz0   = dxidx(3,i,j,k)
		detadx0  = dxidx(4,i,j,k)
		detady0  = dxidx(5,i,j,k)
		detadz0  = dxidx(6,i,j,k)
		dzetadx0 = dxidx(7,i,j,k)
		dzetady0 = dxidx(8,i,j,k)
		dzetadz0 = dxidx(9,i,j,k)	
		inv_j0   = inv_j(  i,j,k)
				
		deltanut = (tempa(i,j,k) + tempb(i,j,k) + tempc(i,j,k))/inv_j0
				
		dnutdx   = dnutdxi(i,j,k)*dxidx0 + dnutdeta(i,j,k)*detadx0 + dnutdzeta(i,j,k)*dzetadx0
		dnutdy   = dnutdxi(i,j,k)*dxidy0 + dnutdeta(i,j,k)*detady0 + dnutdzeta(i,j,k)*dzetady0
		dnutdz   = dnutdxi(i,j,k)*dxidz0 + dnutdeta(i,j,k)*detadz0 + dnutdzeta(i,j,k)*dzetadz0
				
		diff = -(cb2/(re*dta))*nut0*deltanut &
			   +(1.d0/(re*dta))*(temp(i,j,k)*deltanut + tempx(i,j,k)*dnutdx + tempy(i,j,k)*dnutdy + tempz(i,j,k)*dnutdz)
			         
		!!diff = max(diff, 0.d0)
		sa_rhs(i,j,k) = sa_rhs(i,j,k) + diff
	end do
	end do
	end do
	!!*
	
	deallocate(dnutdxi,dnutdeta,dnutdzeta)
	deallocate(tempalpha,tempbeta,tempgamma)
	deallocate(tempa,tempb,tempc)
	deallocate(temp)
	deallocate(tempxi,tempeta,tempzeta)
	deallocate(tempx,tempy,tempz)

	return
end subroutine sa_rhs2


!!=========================================================
!! source terms of sa equation                             
!!                                                         
!! author: jiamin xu                                       
!! date:   2013.01.27                                      
!!=========================================================
subroutine sa_rhs3(sa_rhs,nut,pri_v,amu,vor,length,diasd,re,&
	               is, ie, js, je, ks, ke,  &
	               is0,ie0,js0,je0,ks0,ke0, &
	               is1,ie1,js1,je1,ks1,ke1)
	use sa_var
	use flag_var
	implicit none
	!!*********************************************************
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: re
	real*8  :: diasd (  is:ie,js:je,ks:ke)
	real*8  :: length(  is:ie,js:je,ks:ke)
	real*8  :: vor   (  is:ie,js:je,ks:ke)
	real*8  :: amu   (3,is:ie,js:je,ks:ke)
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: nut   (  is:ie,js:je,ks:ke)
	real*8  :: sa_rhs(  is:ie,js:je,ks:ke)
	!!*********************************************************
	integer :: i,j,k
	real*8  :: d,nut0,nu0,omega,dist0
	real*8  :: aaa,bbb,ccc
	real*8  :: sp1,sd1
	real*8  :: sp,sd
	real*8  :: spn,sdn
	real*8  :: leixn,fv1n,fv2n,ft2n,sn,gn,rn,fwn	
	real*8  :: div
	
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		d     = pri_v  (1,i,j,k)
		nut0  = nut    (  i,j,k)
		nu0   = amu    (1,i,j,k)/d
		omega = abs(vor(  i,j,k)) 
		dist0 = length (  i,j,k)
		
		leix  = max(nut0/nu0,0.001d0)
		fv1   = (leix**3)/(leix**3 + cv1**3)
		fv2   = 1.d0 - leix/(1.d0 + leix*fv1)
		
		s     = omega + nut0*fv2/(re*kar**2*dist0**2)
		s     = max(s,0.3d0*omega)
		
		div = s*re*kar**2*dist0**2
		div = max(div, 1.0d-6)
		r     = nut0*div !!nut0/(s*re*kar**2*dist0**2)
		r     = min(r,10.d0)
		g     = r + cw2*(r**6 - r)
		fw    = g*((1.d0 + cw3**6)/(g**6 + cw3**6))**(1.d0/6.d0)
		ft2   = ct3*exp(-ct4*leix**2)
		
		!!destruction term
		aaa   = cb1*((1.d0 - ft2)*fv2 + ft2)/(kar**2)
		bbb   = cw1*fw
		ccc   = nut0/dist0**2
		sd1   = (aaa - bbb)*ccc/re
		sd    = sd1*nut0
		!!production term
		sp    = cb1*(1.d0-ft2)*omega*nut0
		
		!!added into right hand side term
		sa_rhs(i,j,k) = sa_rhs(i,j,k)  + sd + sp
		!!*
		!!print *, sd+sp
        
		!!get fv2n
		leixn =  1.d0/nu0
		fv1n  = (3.d0*leix**2*leixn*(leix**3+cv1**3) - leix**3*(3.d0*leix**2*leixn))/(leix**3+cv1**3)**2
		fv2n  = (leixn*(1.d0+leix*fv1) - leix*(leixn*fv1+leix*fv1n))/(1.d0+leix*fv1)**2
		!!*
		
		!!get ft2n
		ft2n  = -2.d0*ct4*leix*leixn*ft2
		!!*
		
		!!get fwn
		sn    = 1.d0/(re*kar**2*dist0**2)*(fv2 + fv2n*nut0)
		div   = s**2*(re*kar**2*dist0**2)
		div   = max(div,1.d-6)
		rn    = (s-sn*nut0)/div !!(s-sn*nut0)/s**2/(re*kar**2*dist0**2)
		gn    = (1.d0 + cw2*6.d0*r**5 - cw2)*rn
		fwn   = -(1.d0/6.d0)*((g**(-6.d0)+cw3**(-6.d0))/(1.d0+cw3**(-6.d0)))**(-7.d0/6.d0)
		fwn   = fwn/(1.d0+cw3**(-6.d0))*(-6.d0*g**(-7.d0))*gn
		
		spn = cb1*(1.d0-ft2)*omega + cb1*(-ft2n)*omega*nut0
		sdn = (1.d0/re)*(cb1*(fv2  - ft2 *fv2            + ft2 )/kar**2 - cw1*fw )*2.d0*nut0/dist0**2 
		sdn = (1.d0/re)*(cb1*(fv2n - ft2n*fv2 - ft2*fv2n + ft2n)/kar**2 - cw1*fwn)*(nut0**2 /dist0**2) + sdn
		
		diasd(i,j,k) = min(spn+sdn,0.d0)       
	end do
	end do
	end do
	!!*
    
	return
end subroutine sa_rhs3

!!=========================================================
!! implicit time advancement                               
!! of one equation sa turbulence model                     
!!                                                         
!! author: jiamin xu                                       
!! date:   2013.01.19                                      
!!=========================================================
subroutine sa_lusgs(nut,dt,sa_rhs,pri_v,dxidx,diasd,re,cfl,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1)
	use flag_var
	implicit none
	!!****************************************
	integer :: is,ie,js,je,ks,ke
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: re,cfl
	real*8  :: diasd (  is:ie,js:je,ks:ke)
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: sa_rhs(  is:ie,js:je,ks:ke)
	real*8  :: dt    (  is:ie,js:je,ks:ke)
	real*8  :: nut   (  is:ie,js:je,ks:ke)
	!!****************************************
	integer :: i,j,k
	real*8  :: u,v,w
	real*8  :: dxidx0,dxidy0,dxidz0
	real*8  :: phi
	real*8  :: tempa,tempb,tempc
	real*8  :: dxidx2
	real*8  :: tm
	
	real*8,allocatable,dimension(:,:,:) :: contrau,contrav,contraw
	real*8,allocatable,dimension(:,:,:) :: dnut
	real*8,allocatable,dimension(:,:,:) :: dnut1
	real*8,allocatable,dimension(:,:,:) :: diag
	
	allocate (dnut   (is:ie,js:je,ks:ke))
	allocate (dnut1  (is:ie,js:je,ks:ke))
	allocate (diag   (is:ie,js:je,ks:ke))
	allocate (contrau(is:ie,js:je,ks:ke))
	allocate (contrav(is:ie,js:je,ks:ke))
	allocate (contraw(is:ie,js:je,ks:ke))
	
	dnut  = 0.d0
	dnut1 = 0.d0
	diag  = 0.d0
	
	if(cfl .ge. 1.d3 .and. iflag_time .eq. iflag_steady)then
		phi = 0.d0
	else 
		phi = 1.d0
	end if    
	
	!!get diag
	if (iflag_dimension .eq. iflag_2d) then
		k = 1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			u  = pri_v(2,i,j,k)
			v  = pri_v(3,i,j,k)
			w  = pri_v(4,i,j,k)

			dxidx0 = dxidx(1,i,j,k)
			dxidy0 = dxidx(2,i,j,k)
			dxidz0 = dxidx(3,i,j,k)
			dxidx2 = dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0

			contrau(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w
				
			dxidx0 = dxidx(4,i,j,k)
			dxidy0 = dxidx(5,i,j,k)
			dxidz0 = dxidx(6,i,j,k)
		
			dxidx2 = dxidx2 + dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0
			contrav(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w
				
			dxidx0 = dxidx(7,i,j,k)
			dxidy0 = dxidx(8,i,j,k)
			dxidz0 = dxidx(9,i,j,k)
			
			dxidx2 = dxidx2 + dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0
			contraw(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w				
			!!*
				
			!!viscous spectral radious
			tm = 4.d0*dxidx2*nut(i,j,k)/re
			diag(i,j,k) = 1.d0/dt(i,j,k)*phi + abs(contrau(i,j,k)) + abs(contrav(i,j,k)) + abs(contraw(i,j,k)) + tm - diasd(i,j,k)
		end do
		end do
	else if (iflag_dimension .eq. iflag_3d) then
		do k = ks1-1,ke1+1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			u  = pri_v(2,i,j,k)
			v  = pri_v(3,i,j,k)
			w  = pri_v(4,i,j,k)
				
			dxidx0 = dxidx(1,i,j,k)
			dxidy0 = dxidx(2,i,j,k)
			dxidz0 = dxidx(3,i,j,k)
			dxidx2 = dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0

			contrau(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w
				
			dxidx0 = dxidx(4,i,j,k)
			dxidy0 = dxidx(5,i,j,k)
			dxidz0 = dxidx(6,i,j,k)
		
			dxidx2 = dxidx2 + dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0
			contrav(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w
				
			dxidx0 = dxidx(7,i,j,k)
			dxidy0 = dxidx(8,i,j,k)
			dxidz0 = dxidx(9,i,j,k)
			
			dxidx2 = dxidx2 + dxidx0*dxidx0 + dxidy0*dxidy0 + dxidz0*dxidz0
			contraw(i,j,k) = dxidx0*u + dxidy0*v + dxidz0*w				
			!!*
				
			!!viscous spectral radious
			tm = 4.d0*dxidx2*nut(i,j,k)/re
			diag(i,j,k) = 1.d0/dt(i,j,k)*phi + abs(contrau(i,j,k)) + abs(contrav(i,j,k)) + abs(contraw(i,j,k)) + tm - diasd(i,j,k)
		end do
		end do
		end do
	end if
	!!*
	
	if (iflag_dimension .eq. iflag_3d) then
		!!forward sweep
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			tempa = 0.5d0*(contrau(i,j,k) + abs(contrau(i,j,k)))*dnut1(i-1,j  ,k  )
			tempb = 0.5d0*(contrav(i,j,k) + abs(contrav(i,j,k)))*dnut1(i,  j-1,k  )
			tempc = 0.5d0*(contraw(i,j,k) + abs(contraw(i,j,k)))*dnut1(i,  j  ,k-1)
				
			dnut1(i,j,k) = (sa_rhs(i,j,k) + (tempa + tempb + tempc))/diag(i,j,k)
		end do
		end do
		end do
		!!*
	
		!!backward sweep
		do k = ke1,ks1,-1
		do j = je1,js1,-1
		do i = ie1,is1,-1
			tempa = 0.5d0*(contrau(i,j,k) - abs(contrau(i,j,k)))*dnut(i+1,j  ,k  )
			tempb = 0.5d0*(contrav(i,j,k) - abs(contrav(i,j,k)))*dnut(i,  j+1,k  )
			tempc = 0.5d0*(contraw(i,j,k) - abs(contraw(i,j,k)))*dnut(i,  j  ,k+1)
				
			dnut(i,j,k) = dnut1(i,j,k) - (tempa + tempb + tempc)/diag(i,j,k)
		end do
		end do
		end do
		!!*
	else 
		!!forward sweep
		k = 1
		do j = js1,je1
		do i = is1,ie1
			tempa = 0.5d0*(contrau(i,j,k) + abs(contrau(i,j,k)))*dnut1(i-1,j  ,k)
			tempb = 0.5d0*(contrav(i,j,k) + abs(contrav(i,j,k)))*dnut1(i,  j-1,k)
			
			dnut1(i,j,k) = (sa_rhs(i,j,k) + (tempa + tempb))/diag(i,j,k)
		end do
		end do
		!!*
		
		!!backward sweep
		do j = je1,js1,-1
		do i = ie1,is1,-1
			tempa = 0.5d0*(contrau(i,j,k) - abs(contrau(i,j,k)))*dnut(i+1,j  ,k)
			tempb = 0.5d0*(contrav(i,j,k) - abs(contrav(i,j,k)))*dnut(i,  j+1,k)
					
			dnut(i,j,k) = dnut1(i,j,k) - (tempa + tempb)/diag(i,j,k)
		end do
		end do
		!!*
	end if
		
	!!update nu
	!!nu(n+1) = nu(n) + dnu(n)
	do k = ks1, ke1
	do j = js1, je1
	do i = is1, ie1
		nut(i,j,k) = nut(i,j,k) + dnut(i,j,k)
	end do
	end do
	end do
	!!*			

	deallocate(contrau,contrav,contraw)
	deallocate(dnut,dnut1,diag)

	return
end subroutine