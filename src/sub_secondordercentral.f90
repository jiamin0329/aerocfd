subroutine secondordercentral_1d (df,f,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,coord)
	implicit none
	
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
 	integer :: coord
 	real*8  :: f(is:ie,js:je,ks:ke), df(is:ie,js:je,ks:ke)
 	    
 	integer :: i,j,k
	
 	!!for xi direction
 	if (coord .eq. 1) then
 		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if (i .eq. is) then
				df(i,j,k) = f(i+1,j,k) - f(i,j,k)
			else if (i .eq. ie) then
				df(i,j,k) = f(i,j,k) - f(i-1,j,k)
			else 
 				df(i,j,k) = (f(i+1,j,k) - f(i-1,j,k))/2.d0
			end if
		end do
 		end do   
 		end do
 	end if
 	!!*end xi direction

 	!!for eta direction	
 	if (coord .eq. 2) then
 		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if (j .eq. js) then
				df(i,j,k) = f(i,j+1,k) - f(i,j,k)
			else if(j .eq. je) then
				df(i,j,k) = f(i,j,k) - f(i,j-1,k)
			else 
				df(i,j,k) = (f(i,j+1,k) - f(i,j-1,k))/2.d0
			end if    
 		end do
		end do
		end do
 	end if
 	!!*end eta direction

 	!!for zeta direction
 	if (coord .eq. 3) then
		do k = ks0,ke0
 		do j = js0,je0
		do i = is0,ie0
			if (k .eq. ks) then
 				df(i,j,k) = f(i,j,k+1) - f(i,j,k)
			else if (k .eq. ke) then
				df(i,j,k) = f(i,j,k) - f(i,j,k-1) 
			else
 				df(i,j,k) = (f(i,j,k+1) - f(i,j,k-1))/2.d0   
			end if
		end do
		end do
 		end do
 	end if 
 	!!*end zeta direction
 	
 	return
end subroutine

!!=========================================================!!
!! second order central scheme for ndim = n                !!
!!                                                         !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.11.30                                      !!
!!=========================================================!!
subroutine secondordercentral_nd (df,f,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,ndim,coord)
	implicit none
 	
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
    integer :: ndim
	integer :: coord
	real*8  :: f (ndim,is:ie,js:je,ks:ke)
	real*8  :: df(ndim,is:ie,js:je,ks:ke)
 	
	integer :: i,j,k,n
 	
	!!for xi direction
	if (coord .eq. 1) then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
		do n = 1, ndim
			if (i .eq. is) then
				df(n,i,j,k) = f(n,i+1,j,k) - f(n,i,j,k)
			else if (i .eq. ie) then
				df(n,i,j,k) = f(n,i,j,k) - f(n,i-1,j,k)
			else
				df(n,i,j,k) = (f(n,i+1,j,k) - f(n,i-1,j,k))/2.d0
			end if
		end do
		end do
		end do
		end do
	end if
	!!end xi direction
 	
	!!for eta direction    	
	if (coord .eq. 2) then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if (j .eq. js) then
				df(:,i,j,k) = f(:,i,j+1,k) - f(:,i,j,k)
			else if (j .eq. je) then 
				df(:,i,j,k) = f(:,i,j,k) - f(:,i,j-1,k)
			else
				df(:,i,j,k) = (f(:,i,j+1,k) - f(:,i,j-1,k))/2.d0
			end if
		end do
		end do
		end do
	end if
	!!end eta direction
	
	!!for zeta direction    	
	if (coord .eq. 3) then
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			if      (k .eq. ks) then
				df(:,i,j,k) = f(:,i,j,k+1) - f(:,i,j,k)
			else if (k .eq. ke) then
				df(:,i,j,k) = f(:,i,j,k) - f(:,i,j,k-1)
			else
				df(:,i,j,k) = (f(:,i,j,k+1) - f(:,i,j,k-1))/2.d0	
			end if
		end do
		end do
		end do
	end if
	!!end zeta direction 	

 	return
end subroutine