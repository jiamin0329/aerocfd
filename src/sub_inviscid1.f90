!!=========================================================!!
!! inviscid flux to get rhsi                               !!
!! upwind shceme with fds/fvs                              !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.12.12                                      !!
!!=========================================================!!
subroutine inviscid_flux_1(rhsi,q,dxidx,inv_j, &
						   is,ie,js,je,ks,ke, &
						   is0,ie0,js0,je0,ks0,ke0, &
						   is1,ie1,js1,je1,ks1,ke1,gamma)
    use flag_var
    use blk_var
    implicit none
	
    real*8  :: gamma
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: dxidx  (9,is:ie,js:je,ks:ke)
	real*8  :: inv_j  (  is:ie,js:je,ks:ke)
	real*8  :: q      (5,is:ie,js:je,ks:ke)
	real*8  :: rhsi   (5,is:ie,js:je,ks:ke)
	
	integer :: i,j,k
	real*8, dimension(:,:,:,:),allocatable :: ql,qr,ul,ur,flux

	!! for debug
	real*8, dimension(:,:,:  ),allocatable :: temp
	real*8, dimension(:,:,:,:),allocatable :: temp5

	allocate(ql  (5,is:ie,js:je,ks:ke),qr(5,is:ie,js:je,ks:ke))
	allocate(ul  (5,is:ie,js:je,ks:ke),ur(5,is:ie,js:je,ks:ke))
	allocate(flux(5,is:ie,js:je,ks:ke))

	allocate(temp (  is:ie,js:je,ks:ke))
	allocate(temp5(5,is:ie,js:je,ks:ke))
	
	!!initialization of the rhs
	rhsi = 0.d0
	
	ql = 0.d0
	qr = 0.d0
	ul = 0.d0
	ur = 0.d0
	flux = 0.d0
	temp = 0.d0
	temp5 = 0.d0
	!!*

	!!xi direction
	!!reconstruction	
	call thirdorder_weno(ur,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1, 1) !! ritht main block
	call thirdorder_weno(ul,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1,-1) !! left
	!!* 

	!!left and right status variables
	do k = ks1,  ke1
	do j = js1,  je1
	do i = is1-1,ie1
		ql(1,i,j,k) = ul(1,i,j,k)
		ql(2,i,j,k) = ul(2,i,j,k)/ul(1,i,j,k)
		ql(3,i,j,k) = ul(3,i,j,k)/ul(1,i,j,k)
		ql(4,i,j,k) = ul(4,i,j,k)/ul(1,i,j,k)
		ql(5,i,j,k) =(ul(5,i,j,k) - 0.5d0*(ul(2,i,j,k)**2 + ul(3,i,j,k)**2 + ul(4,i,j,k)**2)/ul(1,i,j,k))*(gamma-1.d0)
				
		qr(1,i,j,k) = ur(1,i,j,k)
		qr(2,i,j,k) = ur(2,i,j,k)/ur(1,i,j,k)
		qr(3,i,j,k) = ur(3,i,j,k)/ur(1,i,j,k)
		qr(4,i,j,k) = ur(4,i,j,k)/ur(1,i,j,k)
		qr(5,i,j,k) =(ur(5,i,j,k) - 0.5d0*(ur(2,i,j,k)**2 + ur(3,i,j,k)**2 + ur(4,i,j,k)**2)/ur(1,i,j,k))*(gamma-1.d0)
	end do
	end do
	end do
	!!*

	!!splitting
	call roe_3d (ql,qr,flux,dxidx,inv_j,gamma,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,1)
	!!*
  	
	!!inviscid right hand side term calculation
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		rhsi(:,i,j,k) = -(flux(:,i,j,k) - flux(:,i-1,j,k))
	end do
	end do
	end do
	!!*

	ql = 0.d0
	qr = 0.d0
	ul = 0.d0
	ur = 0.d0
	flux = 0.d0
	temp = 0.d0
	temp5 = 0.d0
	
	!!eta direction
	!!reconstruction
	call thirdorder_weno(ur,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2, 1) !! ritht
	call thirdorder_weno(ul,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2,-1) !! left	
	!!* 
		
	do k = ks1  ,ke1
	do j = js1-1,je1
	do i = is1  ,ie1
		ql(1,i,j,k) = ul(1,i,j,k)
		ql(2,i,j,k) = ul(2,i,j,k)/ul(1,i,j,k)
		ql(3,i,j,k) = ul(3,i,j,k)/ul(1,i,j,k)
		ql(4,i,j,k) = ul(4,i,j,k)/ul(1,i,j,k)
		ql(5,i,j,k) =(ul(5,i,j,k) - 0.5d0*(ul(2,i,j,k)**2 + ul(3,i,j,k)**2 + ul(4,i,j,k)**2)/ul(1,i,j,k))*(gamma-1.d0)
				
		qr(1,i,j,k) = ur(1,i,j,k)
		qr(2,i,j,k) = ur(2,i,j,k)/ur(1,i,j,k)
		qr(3,i,j,k) = ur(3,i,j,k)/ur(1,i,j,k)
		qr(4,i,j,k) = ur(4,i,j,k)/ur(1,i,j,k)
		qr(5,i,j,k) =(ur(5,i,j,k) - 0.5d0*(ur(2,i,j,k)**2 + ur(3,i,j,k)**2 + ur(4,i,j,k)**2)/ur(1,i,j,k))*(gamma-1.d0)
	end do
	end do
	end do
	!!*
	
	!!splitting
	call roe_3d (ql,qr,flux,dxidx,inv_j,gamma,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,2)
	!!*
		
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		rhsi(:,i,j,k) = rhsi(:,i,j,k) - (flux(:,i,j,k) - flux(:,i,j-1,k))
	end do
	end do
	end do
	!!*

	if (iflag_dimension .eq. iflag_3d) then
		ql = 0.d0
		qr = 0.d0
		ul = 0.d0
		ur = 0.d0
		flux = 0.d0
		temp = 0.d0
		temp5 = 0.d0
		
		!!zeta direction
		!!reconstruction
		call thirdorder_weno(ur,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,3, 1) !! ritht
		call thirdorder_weno(ul,q,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,3,-1) !! left	
		!!* 
			
		do k = ks1-1,ke1
		do j = js1  ,je1
		do i = is1  ,ie1
			ql(1,i,j,k) = ul(1,i,j,k)
			ql(2,i,j,k) = ul(2,i,j,k)/ul(1,i,j,k)
			ql(3,i,j,k) = ul(3,i,j,k)/ul(1,i,j,k)
			ql(4,i,j,k) = ul(4,i,j,k)/ul(1,i,j,k)
			ql(5,i,j,k) =(ul(5,i,j,k) - 0.5d0*(ul(2,i,j,k)**2 + ul(3,i,j,k)**2 + ul(4,i,j,k)**2)/ul(1,i,j,k))*(gamma-1.d0)
					
			qr(1,i,j,k) = ur(1,i,j,k)
			qr(2,i,j,k) = ur(2,i,j,k)/ur(1,i,j,k)
			qr(3,i,j,k) = ur(3,i,j,k)/ur(1,i,j,k)
			qr(4,i,j,k) = ur(4,i,j,k)/ur(1,i,j,k)
			qr(5,i,j,k) =(ur(5,i,j,k) - 0.5d0*(ur(2,i,j,k)**2 + ur(3,i,j,k)**2 + ur(4,i,j,k)**2)/ur(1,i,j,k))*(gamma-1.d0)
		end do
		end do
		end do
		!!*
		
		!!splitting
		call roe_3d (ql,qr,flux,dxidx,inv_j,gamma,is,ie,js,je,ks,ke,is1,ie1,js1,je1,ks1,ke1,3)
		!!*
			
		do k = ks1,ke1
		do j = js1,je1
		do i = is1,ie1
			rhsi(:,i,j,k) = rhsi(:,i,j,k) - (flux(:,i,j,k) - flux(:,i,j,k-1))
		end do
		end do
		end do
		!!*
	end if

	deallocate(ql,qr,ul,ur,flux)
	deallocate(temp,temp5)

	return
end subroutine
