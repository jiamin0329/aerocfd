!!*********************************************************
!! viscous flux to get rhsv
!! for both 2d/3d case                                
!! 2nd central discretization                              
!!                                                         
!! author: jiamin xu                                       
!! date:   2012.11.22                                      
!!*********************************************************
subroutine viscous_flux(rhsv,pri_v,amu,vor,dudx,tao,dxidx,inv_j,tinf,prl,prt,gamma,cv,re, &
						is,ie,js,je,ks,ke, &
						is0,ie0,js0,je0,ks0,ke0, &
						is1,ie1,js1,je1,ks1,ke1,mm)
	use flag_var
	implicit none
	!!***************************************
	integer :: mm
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: tinf,prl,prt,gamma,cv,re
	real*8  :: inv_j(  is:ie,js:je,ks:ke  )
	real*8  :: dxidx(9,is:ie,js:je,ks:ke)
	real*8  :: tao  (6,is:ie,js:je,ks:ke  )
	real*8  :: dudx (9,is:ie,js:je,ks:ke)
	real*8  :: vor  (  is:ie,js:je,ks:ke  )
	real*8  :: amu  (3,is:ie,js:je,ks:ke  )
	real*8  :: pri_v(7,is:ie,js:je,ks:ke  )
	real*8  :: rhsv (5,is:ie,js:je,ks:ke  )
	!!***************************************
	integer :: i,j,k,n
	real*8  :: u0,v0,w0,t0
	real*8  :: ki
	real*8  :: tx,ty,tz
	character(len = 180) :: debugFile
	
	real*8,dimension(:,:,:),  allocatable :: u,v,w,t
	real*8,dimension(:,:,:),  allocatable :: txi,teta,tzeta
	real*8,dimension(:,:,:),  allocatable :: uxi,ueta,uzeta
	real*8,dimension(:,:,:),  allocatable :: vxi,veta,vzeta
	real*8,dimension(:,:,:),  allocatable :: wxi,weta,wzeta
	real*8,dimension(:,:,:),  allocatable :: bx,by,bz
	real*8,dimension(:,:,:,:),allocatable :: fv,gv,hv,fvxi,gveta,hvzeta
	
	allocate(u     (  is:ie,js:je,ks:ke))
	allocate(v     (  is:ie,js:je,ks:ke))
	allocate(w     (  is:ie,js:je,ks:ke))
	allocate(t     (  is:ie,js:je,ks:ke))
	allocate(txi   (  is:ie,js:je,ks:ke))
	allocate(teta  (  is:ie,js:je,ks:ke))
	allocate(tzeta (  is:ie,js:je,ks:ke))
	allocate(uxi   (  is:ie,js:je,ks:ke))
	allocate(ueta  (  is:ie,js:je,ks:ke))
	allocate(uzeta (  is:ie,js:je,ks:ke))
	allocate(vxi   (  is:ie,js:je,ks:ke))
	allocate(veta  (  is:ie,js:je,ks:ke))
	allocate(vzeta (  is:ie,js:je,ks:ke))
	allocate(wxi   (  is:ie,js:je,ks:ke))
	allocate(weta  (  is:ie,js:je,ks:ke))
	allocate(wzeta (  is:ie,js:je,ks:ke))
    
	allocate(bx    (  is:ie,js:je,ks:ke))
	allocate(by    (  is:ie,js:je,ks:ke))
	allocate(bz    (  is:ie,js:je,ks:ke))
	
	allocate(fv    (5,is:ie,js:je,ks:ke))
	allocate(gv    (5,is:ie,js:je,ks:ke))
	allocate(hv    (5,is:ie,js:je,ks:ke))
	allocate(fvxi  (5,is:ie,js:je,ks:ke))
	allocate(gveta (5,is:ie,js:je,ks:ke))
	allocate(hvzeta(5,is:ie,js:je,ks:ke))

	u    = 0.d0; v     = 0.d0; w      = 0.d0; 
	t    = 0.d0;
	txi  = 0.d0; teta  = 0.d0; tzeta  = 0.d0;
	uxi  = 0.d0; ueta  = 0.d0; uzeta  = 0.d0;
	vxi  = 0.d0; veta  = 0.d0; vzeta  = 0.d0;
	wxi  = 0.d0; weta  = 0.d0; wzeta  = 0.d0;
	bx   = 0.d0; by    = 0.d0; bz     = 0.d0;
	fv   = 0.d0; gv    = 0.d0; hv     = 0.d0;
	fvxi = 0.d0; gveta = 0.d0; hvzeta = 0.d0;

	if (iflag_dimension .eq. iflag_3d) then
		do k = ks1-1,ke1+1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			u(i,j,k) = pri_v(2,i,j,k)
			v(i,j,k) = pri_v(3,i,j,k)
			w(i,j,k) = pri_v(4,i,j,k)
			t(i,j,k) = pri_v(6,i,j,k)
		end do
		end do
		end do
	else 
		k = 1
		do j = js1-1,je1+1
		do i = is1-1,ie1+1
			u(i,j,k) = pri_v(2,i,j,k)
			v(i,j,k) = pri_v(3,i,j,k)
			w(i,j,k) = pri_v(4,i,j,k)
			t(i,j,k) = pri_v(6,i,j,k)
		end do
		end do
	end if
	!!*
	
	!! xi direction
	!!call secondordercentral_1d(txi, t,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1)
	!!call secondordercentral_1d(uxi, u,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1)
	!!call secondordercentral_1d(vxi, v,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1)

	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		if      (i .eq. is0) then
			txi(i,j,k) = t(i+1,j,k) - t(i,j,k)
			uxi(i,j,k) = u(i+1,j,k) - u(i,j,k)
			vxi(i,j,k) = v(i+1,j,k) - v(i,j,k)
		else if (i .eq. ie0) then
			txi(i,j,k) = t(i,j,k) - t(i-1,j,k)
			uxi(i,j,k) = u(i,j,k) - u(i-1,j,k)
			vxi(i,j,k) = v(i,j,k) - v(i-1,j,k)
		else 
			txi(i,j,k) = (t(i+1,j,k) - t(i-1,j,k))/2.d0
			uxi(i,j,k) = (u(i+1,j,k) - u(i-1,j,k))/2.d0
			vxi(i,j,k) = (v(i+1,j,k) - v(i-1,j,k))/2.d0
		end if
	end do
 	end do   
 	end do

	!! eta direction
	!!call secondordercentral_1d(teta,t,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2)
	!!call secondordercentral_1d(ueta,u,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2)
	!!call secondordercentral_1d(veta,v,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2)
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		if (j .eq. js) then
			teta(i,j,k) = t(i,j+1,k) - t(i,j,k)
			ueta(i,j,k) = u(i,j+1,k) - u(i,j,k)
			veta(i,j,k) = v(i,j+1,k) - v(i,j,k)
		else if(j .eq. je) then
			teta(i,j,k) = t(i,j,k) - t(i,j-1,k)
			ueta(i,j,k) = u(i,j,k) - u(i,j-1,k)
			veta(i,j,k) = v(i,j,k) - v(i,j-1,k)
		else 
			teta(i,j,k) = (t(i,j+1,k) - t(i,j-1,k))/2.d0
			ueta(i,j,k) = (u(i,j+1,k) - u(i,j-1,k))/2.d0
			veta(i,j,k) = (v(i,j+1,k) - v(i,j-1,k))/2.d0
		end if    
 	end do
	end do
	end do


	!! zeta direction
	if (iflag_dimension .eq. iflag_3d) then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if     (k .eq. ks) then
				teta(i,j,k) = t(i,j,k+1) - t(i,j,k)
				ueta(i,j,k) = u(i,j,k+1) - u(i,j,k)
				veta(i,j,k) = v(i,j,k+1) - v(i,j,k)
			else if(k .eq. ke) then
				teta(i,j,k) = t(i,j,k) - t(i,j,k-1)
				ueta(i,j,k) = u(i,j,k) - u(i,j,k-1)
				veta(i,j,k) = v(i,j,k) - v(i,j,k-1)
			else 
				teta(i,j,k) = (t(i,j,k+1) - t(i,j,k-1))/2.d0
				ueta(i,j,k) = (u(i,j,k+1) - u(i,j,k-1))/2.d0
				veta(i,j,k) = (v(i,j,k+1) - v(i,j,k-1))/2.d0
			end if    
		 end do
		end do
		end do
	end if

	if (isDebug .eq. 1) then
		write(debugFile,"('result/debug_vis_ueta_'I10.10'.dat')"),mm
		write(*,*) "Writing ", debugFile
		call DEBUG_OUTPUT_2D_1(debugFile,mm,ueta,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
	end if
		
	wxi   = 0.d0
	weta  = 0.d0
	tzeta = 0.d0
	uzeta = 0.d0
	vzeta = 0.d0
	wzeta = 0.d0
	!!*
	
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		u0 = u(i,j,k)
		v0 = v(i,j,k)
		w0 = w(i,j,k)
		t0 = t(i,j,k)
				
		amu (1,i,j,k) = (1.d0+124.d0/tinf)/(t0 + 124.d0/tinf)*t0**1.5d0
		amu (3,i,j,k) = amu(1,i,j,k) + amu(2,i,j,k)
				
		ki = gamma*cv*(amu(1,i,j,k)/prl + amu(2,i,j,k)/prt)
		!!=====================================================================================
		tx = txi  (i,j,k)*dxidx(1,i,j,k) &
		   + teta (i,j,k)*dxidx(4,i,j,k) &
		   + tzeta(i,j,k)*dxidx(7,i,j,k)
		ty = txi  (i,j,k)*dxidx(2,i,j,k) &
		   + teta (i,j,k)*dxidx(5,i,j,k) &
		   + tzeta(i,j,k)*dxidx(8,i,j,k)
		tz = txi  (i,j,k)*dxidx(3,i,j,k) &
		   + teta (i,j,k)*dxidx(6,i,j,k) &
		   + tzeta(i,j,k)*dxidx(9,i,j,k)
		!!=====================================================================================
		dudx(1,i,j,k) = uxi  (i,j,k)*dxidx(1,i,j,k) &
				      + ueta (i,j,k)*dxidx(4,i,j,k) &
				      + uzeta(i,j,k)*dxidx(7,i,j,k)
		dudx(2,i,j,k) = uxi  (i,j,k)*dxidx(2,i,j,k) &
				      + ueta (i,j,k)*dxidx(5,i,j,k) &
				      + uzeta(i,j,k)*dxidx(8,i,j,k)  
		dudx(3,i,j,k) = uxi  (i,j,k)*dxidx(3,i,j,k) &
				      + ueta (i,j,k)*dxidx(6,i,j,k) &
				      + uzeta(i,j,k)*dxidx(9,i,j,k)  
				                 			                 
		dudx(4,i,j,k) = vxi  (i,j,k)*dxidx(1,i,j,k) &
				      + veta (i,j,k)*dxidx(4,i,j,k) &
				      + vzeta(i,j,k)*dxidx(7,i,j,k)
		dudx(5,i,j,k) = vxi  (i,j,k)*dxidx(2,i,j,k) &
				      + veta (i,j,k)*dxidx(5,i,j,k) &
				      + vzeta(i,j,k)*dxidx(8,i,j,k)  
		dudx(6,i,j,k) = vxi  (i,j,k)*dxidx(3,i,j,k) &
				      + veta (i,j,k)*dxidx(6,i,j,k) &
				      + vzeta(i,j,k)*dxidx(9,i,j,k)
				                 
		dudx(7,i,j,k) = wxi  (i,j,k)*dxidx(1,i,j,k) &
				      + weta (i,j,k)*dxidx(4,i,j,k) &
				      + wzeta(i,j,k)*dxidx(7,i,j,k)
		dudx(8,i,j,k) = wxi  (i,j,k)*dxidx(2,i,j,k) &
				      + weta (i,j,k)*dxidx(5,i,j,k) &
				      + wzeta(i,j,k)*dxidx(8,i,j,k)  
		dudx(9,i,j,k) = wxi  (i,j,k)*dxidx(3,i,j,k) &
				      + weta (i,j,k)*dxidx(6,i,j,k) &
				      + wzeta(i,j,k)*dxidx(9,i,j,k)					                 
		!!=====================================================================================
		tao(1,i,j,k) = amu(3,i,j,k)*(4.0d0/3.0d0*dudx(1,i,j,k) - 2.0d0/3.0d0*(dudx(5,i,j,k) + dudx(9,i,j,k)))
		tao(2,i,j,k) = amu(3,i,j,k)*(dudx(2,i,j,k) + dudx(4,i,j,k))
		tao(3,i,j,k) = amu(3,i,j,k)*(dudx(7,i,j,k) + dudx(3,i,j,k))
				
		tao(4,i,j,k) = amu(3,i,j,k)*(4.0d0/3.0d0*dudx(5,i,j,k) - 2.0d0/3.0d0*(dudx(1,i,j,k) + dudx(9,i,j,k)))
		tao(5,i,j,k) = amu(3,i,j,k)*(dudx(8,i,j,k) + dudx(6,i,j,k))
				
		tao(6,i,j,k) = amu(3,i,j,k)*(4.0d0/3.0d0*dudx(9,i,j,k) - 2.0d0/3.0d0*(dudx(1,i,j,k) + dudx(5,i,j,k)))
		!!=====================================================================================
		bx(i,j,k) = u0*tao(1,i,j,k) + v0*tao(2,i,j,k) + w0*tao(3,i,j,k) + ki*tx
		by(i,j,k) = u0*tao(2,i,j,k) + v0*tao(4,i,j,k) + w0*tao(5,i,j,k) + ki*ty
		bz(i,j,k) = u0*tao(3,i,j,k) + v0*tao(5,i,j,k) + w0*tao(6,i,j,k) + ki*tz
	end do
	end do
	end do
	!!*
	
	!!get vorticity
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		vor(i,j,k) = sqrt((dudx(6,i,j,k) - dudx(8,i,j,k))**2 &
				         +(dudx(3,i,j,k) - dudx(7,i,j,k))**2 &
				         +(dudx(2,i,j,k) - dudx(4,i,j,k))**2)
	end do
	end do
	end do
	!!*

	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		fv(1,i,j,k) = 0.d0
		gv(1,i,j,k) = 0.d0
		hv(1,i,j,k) = 0.d0
				
		fv(2,i,j,k) = (dxidx(1,i,j,k)*tao(1,i,j,k) + dxidx(2,i,j,k)*tao(2,i,j,k) + dxidx(3,i,j,k)*tao(3,i,j,k))*inv_j(i,j,k)
		gv(2,i,j,k) = (dxidx(4,i,j,k)*tao(1,i,j,k) + dxidx(5,i,j,k)*tao(2,i,j,k) + dxidx(6,i,j,k)*tao(3,i,j,k))*inv_j(i,j,k)
		hv(2,i,j,k) = (dxidx(7,i,j,k)*tao(1,i,j,k) + dxidx(8,i,j,k)*tao(2,i,j,k) + dxidx(9,i,j,k)*tao(3,i,j,k))*inv_j(i,j,k)

		fv(3,i,j,k) = (dxidx(1,i,j,k)*tao(2,i,j,k) + dxidx(2,i,j,k)*tao(4,i,j,k) + dxidx(3,i,j,k)*tao(5,i,j,k))*inv_j(i,j,k)
		gv(3,i,j,k) = (dxidx(4,i,j,k)*tao(2,i,j,k) + dxidx(5,i,j,k)*tao(4,i,j,k) + dxidx(6,i,j,k)*tao(5,i,j,k))*inv_j(i,j,k)
		hv(3,i,j,k) = (dxidx(7,i,j,k)*tao(2,i,j,k) + dxidx(8,i,j,k)*tao(4,i,j,k) + dxidx(9,i,j,k)*tao(5,i,j,k))*inv_j(i,j,k)

		fv(4,i,j,k) = (dxidx(1,i,j,k)*tao(3,i,j,k) + dxidx(2,i,j,k)*tao(5,i,j,k) + dxidx(3,i,j,k)*tao(6,i,j,k))*inv_j(i,j,k)
		gv(4,i,j,k) = (dxidx(4,i,j,k)*tao(3,i,j,k) + dxidx(5,i,j,k)*tao(5,i,j,k) + dxidx(6,i,j,k)*tao(6,i,j,k))*inv_j(i,j,k)
		hv(4,i,j,k) = (dxidx(7,i,j,k)*tao(3,i,j,k) + dxidx(8,i,j,k)*tao(5,i,j,k) + dxidx(9,i,j,k)*tao(6,i,j,k))*inv_j(i,j,k)

		fv(5,i,j,k) = (dxidx(1,i,j,k)*bx(i,j,k)    + dxidx(2,i,j,k)*by(i,j,k)    + dxidx(3,i,j,k)*bz(i,j,k)   )*inv_j(i,j,k)
		gv(5,i,j,k) = (dxidx(4,i,j,k)*bx(i,j,k)    + dxidx(5,i,j,k)*by(i,j,k)    + dxidx(6,i,j,k)*bz(i,j,k)   )*inv_j(i,j,k)
		hv(5,i,j,k) = (dxidx(7,i,j,k)*bx(i,j,k)    + dxidx(8,i,j,k)*by(i,j,k)    + dxidx(9,i,j,k)*bz(i,j,k)   )*inv_j(i,j,k)
	end do
	end do
	end do
	!!*
	
	!!call secondordercentral_nd(fvxi, fv,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,5,1)
	!!call secondordercentral_nd(gveta,gv,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,5,2)
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
	do n = 1, 5
		if (i .eq. is) then
			fvxi(n,i,j,k) = fv(n,i+1,j,k) - fv(n,i,j,k)
		else if (i .eq. ie) then
			fvxi(n,i,j,k) = fv(n,i,j,k) - fv(n,i-1,j,k)
		else
			fvxi(n,i,j,k) = (fv(n,i+1,j,k) - fv(n,i-1,j,k))/2.d0
		end if
	end do
	end do
	end do
	end do

	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		if (j .eq. js) then
			gveta(:,i,j,k) = gv(:,i,j+1,k) - gv(:,i,j,k)
		else if (j .eq. je) then 
			gveta(:,i,j,k) = gv(:,i,j,k) - gv(:,i,j-1,k)
		else
			gveta(:,i,j,k) = (gv(:,i,j+1,k) - gv(:,i,j-1,k))/2.d0
		end if
	end do
	end do
	end do

	if (iflag_dimension .eq. iflag_3d) then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (k .eq. ks) then
				hvzeta(:,i,j,k) = hv(:,i,j,k+1) - hv(:,i,j,k)
			else if (k .eq. ke) then 
				hvzeta(:,i,j,k) = hv(:,i,j,k) - hv(:,i,j,k-1)
			else
				hvzeta(:,i,j,k) = (hv(:,i,j,k+1) - hv(:,i,j,k-1))/2.d0
			end if
		end do
		end do
		end do
	else 
		hvzeta = 0.d0
	end if

	if (isDebug .eq. 1) then
		print *, "***", is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1
		print *, "***", fv(2,246,2,1), fv(2,248,2,1), fvxi(2,247,2,1)
		print *, "***", fv(2,  0,2,1), fv(2,  2,2,1), fvxi(2,  1,2,1)

		write(debugFile,"('result/debug_vis_fv_'I10.10'.dat')"),mm
		write(*,*) "Writing ", debugFile
		call DEBUG_OUTPUT_2D_5(debugFile,mm,fv,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		write(debugFile,"('result/debug_vis_gv_'I10.10'.dat')"),mm
		write(*,*) "Writing ", debugFile
		call DEBUG_OUTPUT_2D_5(debugFile,mm,gv,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		write(debugFile,"('result/debug_vis_fvxi_'I10.10'.dat')"),mm
		write(*,*) "Writing ", debugFile
		call DEBUG_OUTPUT_2D_5(debugFile,mm,fvxi,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		write(debugFile,"('result/debug_vis_gveta_'I10.10'.dat')"),mm
		write(*,*) "Writing ", debugFile
		call DEBUG_OUTPUT_2D_5(debugFile,mm,gveta,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
	end if
	!!*
	
	!!get rhsv
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		rhsv(1,i,j,k) = (fvxi(1,i,j,k) + gveta(1,i,j,k) + hvzeta(1,i,j,k))/re
		rhsv(2,i,j,k) = (fvxi(2,i,j,k) + gveta(2,i,j,k) + hvzeta(2,i,j,k))/re
		rhsv(3,i,j,k) = (fvxi(3,i,j,k) + gveta(3,i,j,k) + hvzeta(3,i,j,k))/re
		rhsv(4,i,j,k) = (fvxi(4,i,j,k) + gveta(4,i,j,k) + hvzeta(4,i,j,k))/re
		rhsv(5,i,j,k) = (fvxi(5,i,j,k) + gveta(5,i,j,k) + hvzeta(5,i,j,k))/re
	end do
	end do
	end do
	!!*
	
	deallocate(u,v,w,t)
	deallocate(txi,teta,tzeta)
	deallocate(uxi,ueta,uzeta)
	deallocate(vxi,veta,vzeta)
	deallocate(wxi,weta,wzeta)
	deallocate(bx,by,bz)
	deallocate(fv,gv,hv)
	deallocate(fvxi,gveta,hvzeta)
	
	return
end subroutine
