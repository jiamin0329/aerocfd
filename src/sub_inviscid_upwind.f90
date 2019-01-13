subroutine thirdorder_weno (df,f,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,coord,direction)
	implicit none
	
	integer :: coord,direction
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: df(5,is:ie,js:je,ks:ke),f(5,is:ie,js:je,ks:ke)
	
	integer :: i,j,k
	real*8  :: temp(5,is:ie,js:je,ks:ke)
	real*8  :: eps
	real*8  :: u1(5),u2(5),u3(5)
	real*8  :: s1(5),s2(5)
	real*8  :: a1(5),a2(5)
	real*8  :: w1(5),w2(5)
	
	eps = 1.d-6
	
	!!xi direction
	if (coord .eq. 1)then
		if (direction .eq. 1) then  !!right reconstruction
			do k = ks0,ke0
			do j = js0,je0
			do i = is0-1,ie0
				if (i .eq. ie-1) then
					df(:,i,j,k) = f(:,i+1,j,k)
				else 
					u1(:) = f(:,i  ,j,k)
					u2(:) = f(:,i+1,j,k)
					u3(:) = f(:,i+2,j,k)
							
					s1(:) = (u2(:)-u3(:))**2
					s2(:) = (u1(:)-u2(:))**2
							
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:)) 
							
					df(:,i,j,k) = w1(:)*(-u3(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u1(:)/2.d0)
				end if
			end do
			end do
			end do
		else if(direction .eq. -1) then  !!left reconstruction
			do k = ks0,ke0
			do j = js0,je0
			do i = is0-1,ie0
				if (i .eq. is) then
					df(:,i,j,k) = f(:,i,j,k)
				else 
					u1(:) = f(:,i-1,j,k)
					u2(:) = f(:,i  ,j,k)
					u3(:) = f(:,i+1,j,k)
							
					s1(:) = (u2(:)-u1(:))**2
					s2(:) = (u3(:)-u2(:))**2
							
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:))
							
					df(:,i,j,k) = w1(:)*(-u1(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u3(:)/2.d0)
				end if
			end do
			end do
			end do
		end if
	else if(coord .eq. 2)then
		if (direction .eq. 1) then  !!right reconstruction
			do k = ks0,ke0
			do j = js0-1,je0
			do i = is0,ie0
				if (j .eq. je-1) then
					df(:,i,j,k) = f(:,i,j+1,k)
				else 
					u1(:) = f(:,i,j  ,k)
					u2(:) = f(:,i,j+1,k)
					u3(:) = f(:,i,j+2,k)
							
					s1(:) = (u2(:)-u3(:))**2
					s2(:) = (u1(:)-u2(:))**2
							
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:)) 
							
					df(:,i,j,k) = w1(:)*(-u3(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u1(:)/2.d0)
				end if
			end do
			end do
			end do
		else if(direction .eq. -1) then  !!left reconstruction
			do k = ks0,ke0
			do j = js0-1,je0
			do i = is0,ie0
				if (j .eq. js) then
					df(:,i,j,k) = f(:,i,j,k)
				else 
					u1(:) = f(:,i,j-1,k)
					u2(:) = f(:,i,j  ,k)
					u3(:) = f(:,i,j+1,k)
							
					s1(:) = (u2(:)-u1(:))**2
					s2(:) = (u3(:)-u2(:))**2
							
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:))
							
					df(:,i,j,k) = w1(:)*(-u1(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u3(:)/2.d0)
				end if
			end do
			end do
			end do
		end if
	else if(coord .eq. 3)then
		if (direction .eq. 1) then  !!right reconstruction
			do k = ks0-1,ke0
			do j = js0,je0
			do i = is0,ie0
				if (k .eq. ke-1) then
					df(:,i,j,k) = f(:,i,j,k+1)
				else 
					u1(:) = f(:,i,j,k  )
					u2(:) = f(:,i,j,k+1)
					u3(:) = f(:,i,j,k+2)
							
					s1(:) = (u2(:)-u3(:))**2
					s2(:) = (u1(:)-u2(:))**2
						
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:)) 
							
					df(:,i,j,k) = w1(:)*(-u3(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u1(:)/2.d0)
				end if
			end do
			end do
			end do
		else if(direction .eq. -1) then  !!left reconstruction
			do k = ks0-1,ke0
			do j = js0,je0
			do i = is0,ie0
				if (k .eq. ks) then
					df(:,i,j,k) = f(:,i,j,k)
				else 
					u1(:) = f(:,i,j,k-1)
					u2(:) = f(:,i,j,k  )
					u3(:) = f(:,i,j,k+1)
							
					s1(:) = (u2(:)-u1(:))**2
					s2(:) = (u3(:)-u2(:))**2
							
					a1(:) = 1.d0/(3.d0*(eps + s1(:))**2)
					a2(:) = 2.d0/(3.d0*(eps + s2(:))**2)
							
					w1(:) = a1(:)/(a1(:)+a2(:))
					w2(:) = a2(:)/(a1(:)+a2(:))
							
					df(:,i,j,k) = w1(:)*(-u1(:)/2.d0 + 3.d0*u2(:)/2.d0) + w2(:)*(u2(:)/2.d0 + u3(:)/2.d0)
				end if
			end do
			end do
			end do
		end if
	end if
	!!*
		
	return
end subroutine