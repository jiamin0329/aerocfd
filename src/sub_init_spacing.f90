!!=======================================================================
!! get grid spacing and volume                             							 
!!                                                         							 
!! purpose��get grid spacing and volume in each block      							 
!!                                                         							 
!! input:                                                  							 
!! coordination: x,y,z                                     							 
!! block dimension: is,ie,js,je,ks,ke                      							 
!!                                                         							 
!! output:                                                 							 
!! primitive variables: spci,spcj,spck,spacing,vol         							 
!!                                                         							 
!! author: jiamin xu                                       							 
!! date:   2012.11.22                                      							 
!!=======================================================================
subroutine init_spacing(spci,spcj,spck,spacing,vol,x,y,z,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	implicit none
	!!***********************************
	!!input parameters
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: x      (is:ie,js:je,ks:ke)
	real*8  :: y      (is:ie,js:je,ks:ke)
	real*8  :: z      (is:ie,js:je,ks:ke)
	!!output parameters
	real*8  :: spci   (is:ie,js:je,ks:ke)	
	real*8  :: spcj   (is:ie,js:je,ks:ke)	
	real*8  :: spck   (is:ie,js:je,ks:ke)	
	real*8  :: spacing(is:ie,js:je,ks:ke)	
	real*8  :: vol    (is:ie,js:je,ks:ke)	
	!!***********************************
	integer :: i,j,k
	real*8  :: dx,dy,dz
	!!***********************************

	spci = 0.d0; spcj = 0.d0 ; spck = 0.d0;
	spacing = 0.d0; vol = 0.d0;

	!!grid space calculation
	!!xi direction grid spacing
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0-1
		dx = x(i+1,j,k) - x(i,j,k)
		dy = y(i+1,j,k) - y(i,j,k)
		dz = z(i+1,j,k) - z(i,j,k)
		spci(i,j,k) = sqrt(dx**2 + dy**2 + dz**2)
	end do
	end do
	end do
	
	do k = ks0,ke0
	do j = js0,je0
		spci(ie0,j,k) = spci(ie0-1,j,k)
	end do
	end do
	!!*end xi direction grid spacing
	
	!!eta direction grid spacing
	do k = ks0,ke0
	do j = js0,je0-1
	do i = is0,ie0
		dx = x(i,j+1,k) - x(i,j,k)
		dy = y(i,j+1,k) - y(i,j,k)
		dz = z(i,j+1,k) - z(i,j,k)
		spcj(i,j,k) = sqrt(dx**2 + dy**2 + dz**2)
	end do
	end do
	end do
	
	do k = ks0,ke0
	do i = is0,ie0
		spcj(i,je0,k) = spcj(i,je0-1,k)
	end do
	end do
	!!*end eta direction grid spacing
	
	!!zeta direction grid spacing
	if ( ke0 .gt. ks0 ) then
		do k = ks0,ke0-1
		do j = js0,je0
		do i = is0,ie0
			dx = x(i,j,k+1) - x(i,j,k)
			dy = y(i,j,k+1) - y(i,j,k)
			dz = z(i,j,k+1) - z(i,j,k)
			spck(i,j,k) = sqrt(dx**2 + dy**2 + dz**2)
		end do
		end do
		end do
	
		do j = js0,je0
		do i = is0,ie0
			spck(i,j,ke0) = spck(i,j,ke0-1)
		end do
		end do
	end if
	!!*end zeta direction grid spacing
	!!***********************************
	!!maximum grid spacing and volume
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		spacing(i,j,k) = max(spci(i,j,k),spcj(i,j,k),spck(i,j,k))
		if ( ke0 .eq. ks0 ) then
			vol(i,j,k) = spci(i,j,k)*spcj(i,j,k)
		else
			vol(i,j,k) = spci(i,j,k)*spcj(i,j,k)*spck(i,j,k)
		end if
	end do
	end do
	end do
	!!***********************************
	return
end subroutine