!!=========================================================!!
!! unsteady first order implicit method                    !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.12.24                                      !!
!!=========================================================!!
subroutine implicit1st(q,dt,pri_v,rhs,rhsi,rhsv,dxidx,inv_j, &
					   is,ie,js,je,ks,ke, &
					   is0,ie0,js0,je0,ks0,ke0, &
					   is1,ie1,js1,je1,ks1,ke1, &
					   gamma, cfl, xx, yy, zz)
	use flag_var
	implicit none
	!!************************************************************************************
	integer :: i,j,k
	integer :: n
	!!************************************************************************************
	real*8  :: gamma,cfl
	integer :: is, ie, js, je, ks, ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	integer :: is1,ie1,js1,je1,ks1,ke1
	real*8  :: inv_j (  is:ie,js:je,ks:ke)
	real*8  :: dxidx (9,is:ie,js:je,ks:ke)
	real*8  :: rhs   (5,is:ie,js:je,ks:ke)
	real*8  :: rhsi  (5,is:ie,js:je,ks:ke)
	real*8  :: rhsv  (5,is:ie,js:je,ks:ke)
	real*8  :: pri_v (7,is:ie,js:je,ks:ke)
	real*8  :: dt    (  is:ie,js:je,ks:ke)
	real*8  :: q     (5,is:ie,js:je,ks:ke)

	real*8  :: xx (is:ie,js:je,ks:ke)
	real*8  :: yy (is:ie,js:je,ks:ke)
	real*8  :: zz (is:ie,js:je,ks:ke)
	!!primitive variables:
	!!d  ==> density
	!!u  ==> velocity x, v ==> velocity y, w ==> velocity z
	!!p  ==> pressure
	!!h  ==> enthalpy
	!!c  ==> soundspeed
	!!q2 ==> 0.5*(u^2 + v^2 + w^2)
	real*8 :: d,u,v,w,p,h,c,q2
	real*8 :: qq1,qq2,qq3,qq4,qq5
	real*8 :: dxidx0,dxidy0,dxidz0
	real*8 :: temp
	real*8 :: cu
	real*8 :: phi
	
	real*8,dimension(:,:,:)  ,allocatable :: diag
	real*8,dimension(:,:,:)  ,allocatable :: ra,rb,rc
	real*8,dimension(:,:,:,:),allocatable :: dq,dqstar
	
	!!jacobian matrix and splitted term
	!!jm_a = pf/pq  jm_a = jm_ap + jm_am
	!!jm_b = pg/pq  jm_b = jm_bp + jm_bm
	!!jm_c = ph/pq  jm_c = jm_cp + jm_cm
	real*8 :: jm_ap(5,5),jm_am(5,5)
	real*8 :: jm_bp(5,5),jm_bm(5,5)
	real*8 :: jm_cp(5,5),jm_cm(5,5)
	
	real*8 :: tempa,tempb,tempc
	
	allocate(diag  (  is:ie,js:je,ks:ke))
	allocate(ra    (  is:ie,js:je,ks:ke))
	allocate(rb    (  is:ie,js:je,ks:ke))
	allocate(rc    (  is:ie,js:je,ks:ke))
	
	allocate(dqstar(5,is:ie,js:je,ks:ke))
	allocate(dq    (5,is:ie,js:je,ks:ke))

	!!get right hand side terms
	!!rhs = rhsi + rhsv
	do k = ks1,ke1
	do j = js1,je1
	do i = is1,ie1
		rhs(:,i,j,k) = rhsi(:,i,j,k) + rhsv(:,i,j,k)
	end do
	end do
	end do
	!!*

	if (isDebug .eq. 1) then
		open (99,file = 'result/debug_rhs.dat')
		write(99,*) "variables = x,y,z,rhs1,rhs2,rhs3,rhs4,rhs5"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
					
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xx(i,j,k),yy(i,j,k),zz(i,j,k), &
			            rhs(1,i,j,k),rhs(2,i,j,k),rhs(3,i,j,k),rhs(4,i,j,k),rhs(5,i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	
	!!if cfl > 1000 then performs a newton iteration
	if(cfl .ge. 1.d3)then
		phi = 0.d0
	else 
		phi = 1.d0
	end if
	!!*

	dq     = 0.d0
	dqstar = 0.d0
	
	!!get diag
	!!do k = ks1-1,ke1+1
	k = 1
	do j = js1-1,je1+1
	do i = is1-1,ie1+1
		!!get primitive variables
		d  = pri_v(1,i,j,k)
		u  = pri_v(2,i,j,k)
		v  = pri_v(3,i,j,k)
		w  = pri_v(4,i,j,k)
		p  = pri_v(5,i,j,k)
		c  = pri_v(7,i,j,k)
				
		dxidx0 = dxidx(1,i,j,k)
		dxidy0 = dxidx(2,i,j,k)
		dxidz0 = dxidx(3,i,j,k)
				
		temp  = dxidx0**2 + dxidy0**2 + dxidz0**2
		temp  = sqrt(temp)
		cu    = dxidx0*u + dxidy0*v + dxidz0*w
				
		ra(i,j,k) = abs(cu) + abs(c*temp)
		!!*

		dxidx0 = dxidx(4,i,j,k)
		dxidy0 = dxidx(5,i,j,k)
		dxidz0 = dxidx(6,i,j,k)
				
		temp  = dxidx0**2 + dxidy0**2 + dxidz0**2
		temp  = sqrt(temp)
		cu    = dxidx0*u + dxidy0*v + dxidz0*w
				
		rb(i,j,k) = abs(cu) + abs(c*temp)
		!!*

		dxidx0 = dxidx(7,i,j,k)
		dxidy0 = dxidx(8,i,j,k)
		dxidz0 = dxidx(9,i,j,k)
				
		temp  = dxidx0**2 + dxidy0**2 + dxidz0**2
		temp  = sqrt(temp)
		cu    = dxidx0*u + dxidy0*v + dxidz0*w
				
		rc(i,j,k) = abs(cu) + abs(c*temp)		
		!!*
		diag(i,j,k) = phi*1.d0/dt(i,j,k) + ra(i,j,k) + rb(i,j,k) + rc(i,j,k)
	end do
	end do
    !!end do
	!!*

	if (isDebug .eq. 1) then
		open (99,file = 'result/debug_diag.dat')
		write(99,*) "variables = x,y,z,diag"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1, "k= ",ke0-ks0+1,"datapacking=point"
					
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xx(i,j,k),yy(i,j,k),zz(i,j,k),diag(i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	
	k = 1  
	!!forward sweep
	do j = js1,je1
	do i = is1,ie1
		!!a+(i-1,j,k)
		dxidx0 = dxidx(1,i-1,j,k)
		dxidy0 = dxidx(2,i-1,j,k)
		dxidz0 = dxidx(3,i-1,j,k)
				
		qq1 = q(1,i-1,j,k)
		qq2 = q(2,i-1,j,k)
		qq3 = q(3,i-1,j,k)
		qq4 = q(4,i-1,j,k)
		qq5 = q(5,i-1,j,k)
				
		call calc_jm (jm_ap,qq1,qq2,qq3,qq4,qq5,dxidx0,dxidy0,dxidz0,gamma)
				
		jm_ap(1,1) = jm_ap(1,1) + ra(i-1,j,k)
		jm_ap(2,2) = jm_ap(2,2) + ra(i-1,j,k)
		jm_ap(3,3) = jm_ap(3,3) + ra(i-1,j,k)
		jm_ap(4,4) = jm_ap(4,4) + ra(i-1,j,k)
		jm_ap(5,5) = jm_ap(5,5) + ra(i-1,j,k)
				
		jm_ap = 0.5d0*jm_ap*inv_j(i-1,j,k)

		!!b+(i,j-1,k)
		dxidx0 = dxidx(4,i,j-1,k)
		dxidy0 = dxidx(5,i,j-1,k)
		dxidz0 = dxidx(6,i,j-1,k)
				
		qq1 = q(1,i,j-1,k)
		qq2 = q(2,i,j-1,k)
		qq3 = q(3,i,j-1,k)
		qq4 = q(4,i,j-1,k)
		qq5 = q(5,i,j-1,k)
				
		call calc_jm (jm_bp,qq1,qq2,qq3,qq4,qq5,dxidx0,dxidy0,dxidz0,gamma)				
				
		jm_bp(1,1) = jm_bp(1,1) + rb(i,j-1,k)
		jm_bp(2,2) = jm_bp(2,2) + rb(i,j-1,k)
		jm_bp(3,3) = jm_bp(3,3) + rb(i,j-1,k)
		jm_bp(4,4) = jm_bp(4,4) + rb(i,j-1,k)
		jm_bp(5,5) = jm_bp(5,5) + rb(i,j-1,k)
				
		jm_bp = 0.5d0*jm_bp*inv_j(i,j-1,k)

		do n = 1,5
			tempa = jm_ap(n,1)*dqstar(1,i-1,j,k) &
				  + jm_ap(n,2)*dqstar(2,i-1,j,k) &
				  + jm_ap(n,3)*dqstar(3,i-1,j,k) &
				  + jm_ap(n,4)*dqstar(4,i-1,j,k) &
				  + jm_ap(n,5)*dqstar(5,i-1,j,k)  
						      
			tempb = jm_bp(n,1)*dqstar(1,i,j-1,k) &
				  + jm_bp(n,2)*dqstar(2,i,j-1,k) &
				  + jm_bp(n,3)*dqstar(3,i,j-1,k) &
				  + jm_bp(n,4)*dqstar(4,i,j-1,k) &
				  + jm_bp(n,5)*dqstar(5,i,j-1,k)  
				  
			dqstar(n,i,j,k) = (rhs(n,i,j,k) + (tempa + tempb))/diag(i,j,k)/inv_j(i,j,k)
		end do
	end do
	end do

	if (isDebug .eq. 1) then
		open (99,file = 'result/debug_dqstar.dat')
		write(99,*) "variables = x,y,z,dqstar1,dqstar2,dqstar3,dqstar4,dqstar5"
		write(99,*) "zone i= ", ie0-is0+1, "j= ", je0-js0+1, "k= ", ke0-ks0+1, "datapacking=point"
						
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xx(i,j,k),yy(i,j,k),zz(i,j,k), &
				        dqstar(1,i,j,k),dqstar(2,i,j,k),dqstar(3,i,j,k),dqstar(4,i,j,k),dqstar(5,i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	!!*
	
	k = 1
	!!backward sweep
	do j = je1,js1,-1
	do i = ie1,is1,-1
		!!a-(i+1,j,k)
		dxidx0 = dxidx(1,i+1,j,k)
		dxidy0 = dxidx(2,i+1,j,k)
		dxidz0 = dxidx(3,i+1,j,k)

		qq1 = q(1,i+1,j,k)
		qq2 = q(2,i+1,j,k)
		qq3 = q(3,i+1,j,k)
		qq4 = q(4,i+1,j,k)
		qq5 = q(5,i+1,j,k)
				
		call calc_jm (jm_am,qq1,qq2,qq3,qq4,qq5,dxidx0,dxidy0,dxidz0,gamma)				
				
		jm_am(1,1) = jm_am(1,1) - ra(i+1,j,k)
		jm_am(2,2) = jm_am(2,2) - ra(i+1,j,k)
		jm_am(3,3) = jm_am(3,3) - ra(i+1,j,k)
		jm_am(4,4) = jm_am(4,4) - ra(i+1,j,k)
		jm_am(5,5) = jm_am(5,5) - ra(i+1,j,k)
				
		jm_am = 0.5d0*jm_am*inv_j(i+1,j,k)

		!!b+(i,j+1,k)
		dxidx0 = dxidx(4,i,j+1,k)
		dxidy0 = dxidx(5,i,j+1,k)
		dxidz0 = dxidx(6,i,j+1,k)
				
		qq1 = q(1,i,j+1,k)
		qq2 = q(2,i,j+1,k)
		qq3 = q(3,i,j+1,k)
		qq4 = q(4,i,j+1,k)
		qq5 = q(5,i,j+1,k)
				
		call calc_jm (jm_bm,qq1,qq2,qq3,qq4,qq5,dxidx0,dxidy0,dxidz0,gamma)
				
		jm_bm(1,1) = jm_bm(1,1) - rb(i,j+1,k)
		jm_bm(2,2) = jm_bm(2,2) - rb(i,j+1,k)
		jm_bm(3,3) = jm_bm(3,3) - rb(i,j+1,k)
		jm_bm(4,4) = jm_bm(4,4) - rb(i,j+1,k)
		jm_bm(5,5) = jm_bm(5,5) - rb(i,j+1,k)
				
		jm_bm = 0.5d0*jm_bm*inv_j(i,j+1,k)				
				
		do n = 1, 5
			tempa = jm_am(n,1)*dq(1,i+1,j,k) &
				  + jm_am(n,2)*dq(2,i+1,j,k) &
				  + jm_am(n,3)*dq(3,i+1,j,k) &
				  + jm_am(n,4)*dq(4,i+1,j,k) &
				  + jm_am(n,5)*dq(5,i+1,j,k)
						    
			tempb = jm_bm(n,1)*dq(1,i,j+1,k) &
				  + jm_bm(n,2)*dq(2,i,j+1,k) &
				  + jm_bm(n,3)*dq(3,i,j+1,k) &
				  + jm_bm(n,4)*dq(4,i,j+1,k) &
				  + jm_bm(n,5)*dq(5,i,j+1,k)
                        
			dq(n,i,j,k) = dqstar(n,i,j,k) - (tempa + tempb)/inv_j(i,j,k)/diag(i,j,k)
			!!print *,  i,j,k,dqstar(n,i,j,k), tempa, tempb ,inv_j(i,j,k), diag(i,j,k)
		end do
	end do
	end do
	!!*
	
	if (isDebug .eq. 1) then
		open (99,file = 'result/debug_dq.dat')
		write(99,*) "variables = x,y,z,dq1,dq2,dq3,dq4,dq5"
		write(99,*) "zone i= ",ie0-is0+1, "j= ",je0-js0+1,"k= ", ke0-ks0+1,"datapacking=point"
						
		do k = ks0,ke0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) xx(i,j,k),yy(i,j,k),zz(i,j,k),&
				        dq(1,i,j,k),dq(2,i,j,k),dq(3,i,j,k),dq(4,i,j,k),dq(5,i,j,k)
		end do
		end do
		end do
		close(99)
	end if
	!!pause

	!!update q
	!!q(n+1) = q(n) + dq(n)
	k = 1
	do j = js1, je1
	do i = is1, ie1
		q(:,i,j,k) = q(:,i,j,k) + dq(:,i,j,k)
	end do
	end do
	!!*    
	
	deallocate(diag)
	deallocate(ra)
	deallocate(rb)
	deallocate(rc)
	deallocate(dqstar)
	deallocate(dq)
	
	return
end subroutine implicit1st
