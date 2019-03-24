!!*************************************************************************************************
!! 	grid transformation                                                                            
!!
!!  purpose: to compute jacobian metrix
!!           VCL and SCL should be satisfied by CMM
!!	                                                                                               
!!  input:                                                                                         
!!  coordination: x(i,j,k),y(i,j,k),z(i,j,k) 
!!  block dimension : is,ie,js,je,ks,ke
!!                                                                                                 
!!	output:                                                                                        
!!  jacobian matrix:  dxidx(1,:,:,:,1): dxidx||dxidx(2,:,:,:,1): detadx||dxidx(3,:,:,:,1): dzetadx 
!!                    dxidx(1,:,:,:,2): dxidy||dxidx(2,:,:,:,2): detady||dxidx(3,:,:,:,2): dzetady 
!!                    dxidx(1,:,:,:,3): dxidz||dxidx(2,:,:,:,3): detadz||dxidx(3,:,:,:,3): dzetadz 
!!  jacobian:         inv_j(i,j,k)                                                                 
!!  dxidx2:           dxidx2(1,:,:,:) = dxidx  *dxidx   + dxidy  *dxidy   + dxidz  *dxidz          
!!                    dxidx2(2,:,:,:) = detadx *detadx  + detady *detady  + detadz *detadz         
!!                    dxidx2(3,:,:,:) = dzetadx*dzetadx + dzetady*dzetady + dzetadz*dzetadz        
!!  alpha/beta/gamma: div xi(i)*div xi(j)   
!!  dxidxdxi:         1,i,j,k,1 ==> d(dxidx/j)/dxi  2,i,j,k,1 ==> d(detadx/j)/deta   3,i,j,k,1 ==> d(dzetadx/j)/dzeta
!!                    1,i,j,k,2 ==> d(dxidy/j)/dxi  2,i,j,k,2 ==> d(detady/j)/deta   3,i,j,k,2 ==> d(dzetady/j)/dzeta
!!                    1,i,j,k,3 ==> d(dxidz/j)/dxi  2,i,j,k,3 ==> d(detadz/j)/deta   3,i,j,k,3 ==> d(dzetadz/j)/dzeta
!!                                                                                                 
!!	Author: Jiamin Xu                                                                              
!!	Date:   2013.06.11                                                                             
!!*************************************************************************************************
subroutine gridtransformation(dxidx,inv_j,alpha,beta,gamma,x,y,z, &
	                          is,ie,js,je,ks,ke,                  &
	                          is0,ie0,js0,je0,ks0,ke0,m0)
	use flag_var
	use mpi_var
	implicit none
	!!**********************************************************
	integer :: m0
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: x     (  is:ie,js:je,ks:ke)
	real*8  :: y     (  is:ie,js:je,ks:ke)
	real*8  :: z     (  is:ie,js:je,ks:ke)
	real*8  :: alpha (3,is:ie,js:je,ks:ke)
	real*8  :: beta  (3,is:ie,js:je,ks:ke)
	real*8  :: gamma (3,is:ie,js:je,ks:ke)
	real*8  :: inv_j (  is:ie,js:je,ks:ke)
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	!!**********************************************************
	integer :: i,j,k
	
	real*8  :: dxdxi0, dxdeta0, dxdzeta0
	real*8  :: dydxi0, dydeta0, dydzeta0
	real*8  :: dzdxi0, dzdeta0, dzdzeta0
		
	real*8  :: dxidx0,  dxidy0,  dxidz0
	real*8  :: detadx0, detady0, detadz0
	real*8  :: dzetadx0,dzetady0,dzetadz0
	real*8  :: inv_j0
	
	real*8,allocatable,dimension(:,:,:,:) :: xyz
	real*8,allocatable,dimension(:,:,:,:) :: dxdxi,dxdeta,dxdzeta

	real*8,allocatable,dimension(:,:,:,:) :: tempxi,  tempeta,   tempzeta
	real*8,allocatable,dimension(:,:,:,:) :: tempxixi,tempetaeta,tempzetazeta

	character(len = 180) :: debugFile
	
	allocate(dxdxi       (3,is:ie,js:je,ks:ke))
	allocate(dxdeta      (3,is:ie,js:je,ks:ke))
	allocate(dxdzeta     (3,is:ie,js:je,ks:ke))
	allocate(xyz         (3,is:ie,js:je,ks:ke))
	allocate(tempxi      (6,is:ie,js:je,ks:ke))
	allocate(tempeta     (6,is:ie,js:je,ks:ke))
	allocate(tempzeta    (6,is:ie,js:je,ks:ke))
	allocate(tempxixi    (6,is:ie,js:je,ks:ke))
	allocate(tempetaeta  (6,is:ie,js:je,ks:ke))
	allocate(tempzetazeta(6,is:ie,js:je,ks:ke))
	!!**********************************************************!!
	!!                   load x,y,z to xyz                      !!
	!!**********************************************************!!
	!! we do NOT have coordinates at corner buffer blocks,
	!! 2D - 4 corner blocks
	!! 3D - 8 corner blocks
	do k = ks,ke
	do j = js,je
	do i = is,ie
		xyz(1,i,j,k) = x(i,j,k)
		xyz(2,i,j,k) = y(i,j,k)
		xyz(3,i,j,k) = z(i,j,k)
	end do
	end do
	end do
	!!*

	if (isDebug .eq. 1) then
		write(debugFile,"('result/DEBUG_xyz_'I10.10'_'I10.10'.dat')"),myid,m0
		write(*,*) "Writing ", debugFile
		
		if      (iflag_dimension .eq. iflag_3d) then
			call DEBUG_OUTPUT_3D_3(debugFile,m0,xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
		else if (iflag_dimension .eq. iflag_2d) then
			call DEBUG_OUTPUT_2D_3(debugFile,m0,xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		end if
	end if

	!!**********************************************************!!
	!!       get dxdxi,dxdeta,dxdzeta using cmpt scheme         !!
	!!**********************************************************!! 
	
	!!xi direction
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		if      (i .eq. is0) then
			dxdxi(:,i,j,k) = xyz(:,i+1,j,k) - xyz(:,i,j,k)
		else if (i .eq. ie0) then
			dxdxi(:,i,j,k) = xyz(:,i,j,k) - xyz(:,i-1,j,k)
		else 
			dxdxi(:,i,j,k) = (xyz(:,i+1,j,k) - xyz(:,i-1,j,k))/2.d0
		end if
	end do
	end do   
	end do

	!! eta direction
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		if      (j .eq. js0) then
			dxdeta(:,i,j,k) = xyz(:,i,j+1,k) - xyz(:,i,j,k)
		else if (j .eq. je0) then
			dxdeta(:,i,j,k) = xyz(:,i,j,k) - xyz(:,i,j-1,k)
		else 
			dxdeta(:,i,j,k) = (xyz(:,i,j+1,k) - xyz(:,i,j-1,k))/2.d0
		end if
	end do
	end do   
	end do

	!!call secondordercentral_nd(dxdxi, xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,3,1)
	!!call secondordercentral_nd(dxdeta,xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,3,2)
	if      (iflag_dimension .eq. iflag_2d)then
		dxdzeta(1,:,:,:) = 0.d0
		dxdzeta(2,:,:,:) = 0.d0
		dxdzeta(3,:,:,:) = 1.d0
	else if (iflag_dimension .eq. iflag_3d)then
		!!call secondordercentral_nd(dxdxi,  xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,3,1)
		!!call secondordercentral_nd(dxdeta, xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,3,2)
		!!call secondordercentral_nd(dxdzeta,xyz,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,3,3)

		!! zeta direction
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (k .eq. ks0) then
				dxdzeta(:,i,j,k) = xyz(:,i,j,k+1) - xyz(:,i,j,k)
			else if (k .eq. ke0) then
				dxdzeta(:,i,j,k) = xyz(:,i,j,k) - xyz(:,i,j,k-1)
			else 
				dxdzeta(:,i,j,k) = (xyz(:,i,j,k+1) - xyz(:,i,j,k-1))/2.d0
			end if
		end do
		end do   
		end do
	end if
	!!*
	
	if (isDebug .eq. 1) then
		write(debugFile,"('result/debug_dxdxi_'I10.10'_'I10.10'.dat')"),myid,m0
		write(*,*) "Writing ", debugFile
		open (99,file = debugFile)
!   
		write(99,*) "variables = x,y,z,dxdxi,dydxi,dzdxi,dxdeta,dydeta,dzdeta,dxdzeta,dydzeta,dzdzeta"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
		
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xyz    (1,i,j,k),xyz    (2,i,j,k),xyz    (3,i,j,k), &
						dxdxi  (1,i,j,k),dxdxi  (2,i,j,k),dxdxi  (3,i,j,k), &
						dxdeta (1,i,j,k),dxdeta (2,i,j,k),dxdeta (3,i,j,k), &
						dxdzeta(1,i,j,k),dxdzeta(2,i,j,k),dxdzeta(3,i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	!!**********************************************************!!

	!!**********************************************************!!
	!!                       get 1/jac                          !!
	!!**********************************************************!! 	
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		dxdxi0 = dxdxi(1,i,j,k); dxdeta0 = dxdeta(1,i,j,k); dxdzeta0 = dxdzeta(1,i,j,k)
		dydxi0 = dxdxi(2,i,j,k); dydeta0 = dxdeta(2,i,j,k); dydzeta0 = dxdzeta(2,i,j,k)
		dzdxi0 = dxdxi(3,i,j,k); dzdeta0 = dxdeta(3,i,j,k); dzdzeta0 = dxdzeta(3,i,j,k)
				
		inv_j (i,j,k) = dxdxi0  *dydeta0 *dzdzeta0 &
				      + dxdeta0 *dydzeta0*dzdxi0   &
				      + dxdzeta0*dydxi0  *dzdeta0  &
				      - dxdzeta0*dydeta0 *dzdxi0   &
				      - dxdxi0  *dydzeta0*dzdeta0  &
				      - dxdeta0 *dydxi0  *dzdzeta0
				
		!!check negative volumes              
		if(inv_j(i,j,k) .le. 0)then
			write(*,* ) "negative volume detected!!!", m0
			write(*,10) xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k)
10			format(" x =",f9.6," y =",f9.6, " z =",f9.6)
		end if
	end do
	end do
	end do
	!!*

	if (isDebug .eq. 1) then
		write(debugFile,"('result/debug_jacobian_'I10.10'_'I10.10'.dat')"),myid,m0
		open (99,file = debugFile)
		write(99,*) "variables = x,y,z,inv_j"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
		
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k),inv_j(i,j,k)
		end do
		end do
		end do
		close(99)
	end if

	!!***************************************************************************!!
	!! get temp                                                                  !!
	!!                                                                           !!
	!! tempxi  (1) = dydzeta*z  tempxi  (2) = dzdzeta*x  tempxi  (3) = dxdzeta*y !!
	!! tempxi  (4) = dydeta *z  tempxi  (5) = dzdeta *x  tempxi  (6) = dxdeta *y !!
	!!                                                                           !!
	!! tempeta (1) = dydzeta*z  tempeta (2) = dzdzeta*x  tempeta (3) = dxdzeta*y !!
	!! tempeta (4) = dydxi  *z  tempeta (5) = dzdxi  *x  tempeta (6) = dxdxi  *y !!
	!!                                                                           !!
	!! tempzeta(1) = dydeta *z  tempzeta(2) = dzdeta *x  tempzeta(3) = dxdeta *y !!
	!! tempzeta(4) = dydxi  *z  tempzeta(5) = dzdxi  *x  tempzeta(6) = dxdxi  *y !!
	!!                                                                           !!
	!!***************************************************************************!!
	if      (iflag_dimension .eq. iflag_2d)then
		k = 1
		do j = js0,je0
		do i = is0,ie0
			dxidx(1,i,j,k) =  dxdeta(2,i,j,k)/inv_j(i,j,k)
			dxidx(2,i,j,k) = -dxdeta(1,i,j,k)/inv_j(i,j,k)
			dxidx(3,i,j,k) =  0.d0
			dxidx(4,i,j,k) = -dxdxi (2,i,j,k)/inv_j(i,j,k)
			dxidx(5,i,j,k) =  dxdxi (1,i,j,k)/inv_j(i,j,k)
			dxidx(6,i,j,k) =  0.d0
			dxidx(7,i,j,k) =  0.d0
			dxidx(8,i,j,k) =  0.d0
			dxidx(9,i,j,k) =  1.d0
		end do
		end do
	else if (iflag_dimension .eq. iflag_3d)then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			tempxi  (1,i,j,k) = dxdzeta(2,i,j,k)*xyz(3,i,j,k)
			tempxi  (2,i,j,k) = dxdzeta(3,i,j,k)*xyz(1,i,j,k)
			tempxi  (3,i,j,k) = dxdzeta(1,i,j,k)*xyz(2,i,j,k)
			tempxi  (4,i,j,k) = dxdeta (2,i,j,k)*xyz(3,i,j,k)
			tempxi  (5,i,j,k) = dxdeta (3,i,j,k)*xyz(1,i,j,k)
			tempxi  (6,i,j,k) = dxdeta (1,i,j,k)*xyz(2,i,j,k)
						
			tempeta (1,i,j,k) = dxdzeta(2,i,j,k)*xyz(3,i,j,k)
			tempeta (2,i,j,k) = dxdzeta(3,i,j,k)*xyz(1,i,j,k)
			tempeta (3,i,j,k) = dxdzeta(1,i,j,k)*xyz(2,i,j,k)
			tempeta (4,i,j,k) = dxdxi  (2,i,j,k)*xyz(3,i,j,k)
			tempeta (5,i,j,k) = dxdxi  (3,i,j,k)*xyz(1,i,j,k)
			tempeta (6,i,j,k) = dxdxi  (1,i,j,k)*xyz(2,i,j,k)				
						
			tempzeta(1,i,j,k) = dxdeta (2,i,j,k)*xyz(3,i,j,k)
			tempzeta(2,i,j,k) = dxdeta (3,i,j,k)*xyz(1,i,j,k)
			tempzeta(3,i,j,k) = dxdeta (1,i,j,k)*xyz(2,i,j,k)
			tempzeta(4,i,j,k) = dxdxi  (2,i,j,k)*xyz(3,i,j,k)
			tempzeta(5,i,j,k) = dxdxi  (3,i,j,k)*xyz(1,i,j,k)
			tempzeta(6,i,j,k) = dxdxi  (1,i,j,k)*xyz(2,i,j,k)
		end do
		end do
		end do
			
		!! xi direction
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (i .eq. is0) then
				tempxixi(:,i,j,k) = tempxi(:,i+1,j,k) - tempxi(:,i,j,k)
			else if (i .eq. ie0) then
				tempxixi(:,i,j,k) = tempxi(:,i,j,k) - tempxi(:,i-1,j,k)
			else 
				tempxixi(:,i,j,k) = (tempxi(:,i+1,j,k) - tempxi(:,i-1,j,k))/2.d0
			end if
		end do
		end do   
		end do

		!! eta direction
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (j .eq. js0) then
				tempetaeta(:,i,j,k) = tempeta(:,i,j+1,k) - tempeta(:,i,j,k)
			else if (j .eq. je0) then
				tempetaeta(:,i,j,k) = tempeta(:,i,j,k) - tempeta(:,i,j-1,k)
			else 
				tempetaeta(:,i,j,k) = (tempeta(:,i,j+1,k) - tempeta(:,i,j-1,k))/2.d0
			end if
		end do
		end do   
		end do

		!! zeta direction
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (k .eq. ks0) then
				tempzetazeta(:,i,j,k) = tempzeta(:,i,j,k+1) - tempzeta(:,i,j,k)
			else if (k .eq. ke0) then
				tempzetazeta(:,i,j,k) = tempzeta(:,i,j,k) - tempzeta(:,i,j,k-1)
			else 
				tempzetazeta(:,i,j,k) = (tempzeta(:,i,j,k+1) - tempzeta(:,i,j,k-1))/2.d0
			end if
		end do
		end do   
		end do

		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			!!dxdxi0 = dxdxi(1,i,j,k); dxdeta0 = dxdeta(1,i,j,k); dxdzeta0 = dxdzeta(1,i,j,k)
			!!dydxi0 = dxdxi(2,i,j,k); dydeta0 = dxdeta(2,i,j,k); dydzeta0 = dxdzeta(2,i,j,k)
			!!dzdxi0 = dxdxi(3,i,j,k); dzdeta0 = dxdeta(3,i,j,k); dzdzeta0 = dxdzeta(3,i,j,k)

			!!dxidx(1,i,j,k) = (dydeta0*dzdzeta0 - dydzeta0*dzdeta0)/inv_j(i,j,k)!!(tempzetazeta(1,i,j,k) - tempetaeta  (1,i,j,k))/inv_j(i,j,k)
			!!dxidx(2,i,j,k) = (dxdzeta0*dzdeta0 - dxdeta0*dzdzeta0)/inv_j(i,j,k)!!(tempzetazeta(2,i,j,k) - tempetaeta  (2,i,j,k))/inv_j(i,j,k)
			!!dxidx(3,i,j,k) = (dxdeta0*dydzeta0 - dxdzeta0*dydeta0)/inv_j(i,j,k)!!(tempzetazeta(3,i,j,k) - tempetaeta  (3,i,j,k))/inv_j(i,j,k)
			!!				
			!!dxidx(4,i,j,k) = (dydzeta0*dzdxi0 - dydxi0*dzdzeta0)/inv_j(i,j,k)!!(tempxixi    (1,i,j,k) - tempzetazeta(4,i,j,k))/inv_j(i,j,k)
			!!dxidx(5,i,j,k) = (dxdxi0*dzdzeta0 - dxdzeta0*dzdxi0)/inv_j(i,j,k)!!(tempxixi    (2,i,j,k) - tempzetazeta(5,i,j,k))/inv_j(i,j,k)
			!!dxidx(6,i,j,k) = (dxdzeta0*dydxi0 - dxdxi0*dydzeta0)/inv_j(i,j,k)!!(tempxixi    (3,i,j,k) - tempzetazeta(6,i,j,k))/inv_j(i,j,k)
			!!				
			!!dxidx(7,i,j,k) = (dydxi0*dzdeta0 - dydeta0*dzdxi0)/inv_j(i,j,k)!!(tempetaeta  (4,i,j,k) - tempxixi    (4,i,j,k))/inv_j(i,j,k)
			!!dxidx(8,i,j,k) = (dxdeta0*dzdxi0 - dxdxi0*dzdeta0)/inv_j(i,j,k)!!(tempetaeta  (5,i,j,k) - tempxixi    (5,i,j,k))/inv_j(i,j,k)
			!!dxidx(9,i,j,k) = (dxdxi0*dydeta0 - dxdeta0*dydxi0)/inv_j(i,j,k)!!(tempetaeta  (6,i,j,k) - tempxixi    (6,i,j,k))/inv_j(i,j,k)

			dxidx(1,i,j,k) = (tempzetazeta(1,i,j,k) - tempetaeta  (1,i,j,k))/inv_j(i,j,k)
			dxidx(2,i,j,k) = (tempzetazeta(2,i,j,k) - tempetaeta  (2,i,j,k))/inv_j(i,j,k)
			dxidx(3,i,j,k) = (tempzetazeta(3,i,j,k) - tempetaeta  (3,i,j,k))/inv_j(i,j,k)
							
			dxidx(4,i,j,k) = (tempxixi    (1,i,j,k) - tempzetazeta(4,i,j,k))/inv_j(i,j,k)
			dxidx(5,i,j,k) = (tempxixi    (2,i,j,k) - tempzetazeta(5,i,j,k))/inv_j(i,j,k)
			dxidx(6,i,j,k) = (tempxixi    (3,i,j,k) - tempzetazeta(6,i,j,k))/inv_j(i,j,k)
							
			dxidx(7,i,j,k) = (tempetaeta  (4,i,j,k) - tempxixi    (4,i,j,k))/inv_j(i,j,k)
			dxidx(8,i,j,k) = (tempetaeta  (5,i,j,k) - tempxixi    (5,i,j,k))/inv_j(i,j,k)
			dxidx(9,i,j,k) = (tempetaeta  (6,i,j,k) - tempxixi    (6,i,j,k))/inv_j(i,j,k)
		end do
		end do
		end do
	end if

	if (isDebug .eq. 1) then
		write(debugFile,"('result/debug_dxidx_'I10.10'_'I10.10'.dat')"),myid,m0
		open (99,file = debugFile)
		write(99,*) "variables = x,y,z,dxidx, dxidy, dxidz, detadx, detady, detadz, dzetadx, dzetady, dzetadz"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
		
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xyz  (1,i,j,k),xyz  (2,i,j,k),xyz  (3,i,j,k), &
						dxidx(1,i,j,k),dxidx(2,i,j,k),dxidx(3,i,j,k), &
						dxidx(4,i,j,k),dxidx(5,i,j,k),dxidx(6,i,j,k), &
						dxidx(7,i,j,k),dxidx(8,i,j,k),dxidx(9,i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	!!*
	
 	!!**********************************************************!!
 	!!               get dxidx2,alpha,beta,gamma                !!
 	!!**********************************************************!! 
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
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
				
		alpha(1,i,j,k) = (dxidx0  *dxidx0   + dxidy0  *dxidy0   + dxidz0  *dxidz0  )*inv_j0
		alpha(2,i,j,k) = (dxidx0  *detadx0  + dxidy0  *detady0  + dxidz0  *detadz0 )*inv_j0
		alpha(3,i,j,k) = (dxidx0  *dzetadx0 + dxidy0  *dzetady0 + dxidz0  *dzetadz0)*inv_j0
				
		beta (1,i,j,k) = (detadx0 *dxidx0   + detady0 *dxidy0   + detadz0 *dxidz0  )*inv_j0
		beta (2,i,j,k) = (detadx0 *detadx0  + detady0 *detady0  + detadz0 *detadz0 )*inv_j0
		beta (3,i,j,k) = (detadx0 *dzetadx0 + detady0 *dzetady0 + detadz0 *dzetadz0)*inv_j0
				
		gamma(1,i,j,k) = (dzetadx0*dxidx0   + dzetady0*dxidy0   + dzetadz0*dxidz0  )*inv_j0
		gamma(2,i,j,k) = (dzetadx0*detadx0  + dzetady0*detady0  + dzetadz0*detadz0 )*inv_j0
		gamma(3,i,j,k) = (dzetadx0*dzetadx0 + dzetady0*dzetady0 + dzetadz0*dzetadz0)*inv_j0
	end do
	end do
	end do

	if (isDebug .eq. 1) then
		write(debugFile,"('result/debug_alpha_beta_gamma_'I10.10'_'I10.10'.dat')"),myid,m0
		open (99,file = debugFile)
		write(99,*) "variables = x,y,z,alpha1,alpha2,alpha3,beta1,beta2,beta3,gamma1,gamma2,gamma3"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
		
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xyz  (1,i,j,k),xyz  (2,i,j,k),xyz  (3,i,j,k), & 
						alpha(1,i,j,k),alpha(2,i,j,k),alpha(3,i,j,k), &
						beta (1,i,j,k),beta (2,i,j,k),beta (3,i,j,k), &
						gamma(1,i,j,k),gamma(2,i,j,k),gamma(3,i,j,k)
		end do
		end do
		end do
		close(99)
	end if

	deallocate(tempxixi,tempetaeta,tempzetazeta)
	deallocate(tempxi,tempeta,tempzeta)
	deallocate(dxdxi,dxdeta,dxdzeta)
	deallocate(xyz)
	!!*
	
	return
end subroutine
