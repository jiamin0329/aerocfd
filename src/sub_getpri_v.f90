!!=========================================================!!
!! get primitive variables                                 !!
!!                                                         !!
!! purpose: to get primitive variables from                !!
!!          the conservative variables                     !!
!!                                                         !!
!!                                                         !!
!! input:                                                  !!
!! conservative variables: q                               !!
!! gas constant: gamma,cv                                  !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!!                                                         !!
!! output:                                                 !!
!! primitive variables: pri_v                              !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.11.22                                      !!
!!=========================================================!!
subroutine get_primitivevariables(pri_v,q,gamma,cv,is,ie,js,je,ks,ke)
	implicit none
	
	integer :: is,ie,js,je,ks,ke
	real*8  :: gamma,cv
	real*8  :: q    (5,is:ie,js:je,ks:ke)
	real*8  :: pri_v(7,is:ie,js:je,ks:ke)
	
	integer :: i,j,k
	real*8  :: d,u,v,w,t,p,c
	real*8  :: q2,e
	
	!!get primitive variables
	do k = ks,ke
	do j = js,je
	do i = is,ie
		d  = q(1,i,j,k)
		u  = q(2,i,j,k)/d
		v  = q(3,i,j,k)/d
		w  = q(4,i,j,k)/d
		q2 = 0.5d0*(u**2 + v**2 + w**2)
		e  = q(5,i,j,k)/d
		t  = (e - q2)/cv
		p  = (q(5,i,j,k) - q2*d)*(gamma-1.d0)
		c  = sqrt(gamma*p/d)
				
		pri_v(1,i,j,k) = d
		pri_v(2,i,j,k) = u
		pri_v(3,i,j,k) = v
		pri_v(4,i,j,k) = w
		pri_v(5,i,j,k) = p
		pri_v(6,i,j,k) = t
		pri_v(7,i,j,k) = c

		if (d .lt. 0.d0) then 
			print *, "Negative density found!!!", i,j,k
		end if
	end do
	end do
	end do
	!!*
	
	return
end subroutine