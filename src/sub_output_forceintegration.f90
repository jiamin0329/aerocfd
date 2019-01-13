!!=========================================================!!
!! cl&cd output                                            !!
!!                                                         !!
!!                                                         !!
!! author: jiamin xu                                       !!
!! date:   2012.12.05                                      !!
!!=========================================================!!
subroutine clcd_output
	use flag_var
	use blk_var
	use ns_const
	use mpi_var
	use index_var
	implicit none
	!!**************************************************************************
	real*8  :: cl,cd
	real*8  :: ffp(3),ffv(3)
	real*8  :: fp(3),fv(3)                 !!pressure force and viscous force
	real*8  :: f_prs(3),f_vis(3)           !!pressure force and viscous force of a single element
	real*8  :: uc,vc,wc,pc,miuc            !!variables of a single element
	real*8  :: dxidxc,dxidyc,dxidzc        !!unit vector normal to the surface
	real*8  :: temp00,temp10,temp01,temp11 !!sqrt(dxidx^2 + dxidy^2 + dxidz^2)
	real*8  :: s,dn                        !!s: surface area; dn: normal distance
	real*8  :: dinf,uinf,pinf,div   
	integer :: face
	integer :: i,j,k
	integer :: ip,jp,kp
	integer :: im,jm,km
	integer :: is, ie, js, je, ks, ke !!block dimension
	real*8  :: x0,x1,y0,y1
	real*8  :: t1,t2,n1,n2
	real*8  :: un
	character(len = 180) :: force_out
	!!*************************************************
	dinf = d_0
	uinf = u_0
	pinf = p_0
	div  = 1.d0/(0.5d0*dinf*uinf*uinf*sref)
	!!initialize the force
	fp   = 0.d0
	fv   = 0.d0
	ffp  = 0.d0
	ffv  = 0.d0
	cl   = 0.d0
	cd   = 0.d0
	!!*
	
	if     (iflag_dimension .eq. iflag_2d)then
		k = 1
		!!compute pressure force and viscous force
		do m0 = 1,blk_loop
			do ksub = 1,blk(m0)%num_subface
				if(blk(m0)%bc(ksub)%blk_t .eq. bc_wall)then
					face = blk(m0)%bc(ksub)%face
					!!wall dimention
					is = blk(m0)%bc(ksub)%is
					ie = blk(m0)%bc(ksub)%ie
					js = blk(m0)%bc(ksub)%js
					je = blk(m0)%bc(ksub)%je
					!!*
					
					if      (face .eq. 1)then
						i  = is
						ip = is+1
						do j = js,je-1
							jp = j+1
							
							uc = (blk(m0)%pri_v(2,ip,j,k) + blk(m0)%pri_v(2,ip,jp,k))*0.5d0
							vc = (blk(m0)%pri_v(3,ip,j,k) + blk(m0)%pri_v(3,ip,jp,k))*0.5d0
							wc = 0.d0
							pc = (blk(m0)%pri_v(5,i ,j,k) + blk(m0)%pri_v(5,i ,jp,k))*0.5d0
							
							temp00 = 1.d0/sqrt(blk(m0)%dxidx(1,ip,j ,k)**2 + blk(m0)%dxidx(2,ip,j ,k)**2 + blk(m0)%dxidx(3,ip,j ,k)**2)  
							temp10 = 1.d0/sqrt(blk(m0)%dxidx(1,ip,jp,k)**2 + blk(m0)%dxidx(2,ip,jp,k)**2 + blk(m0)%dxidx(3,ip,jp,k)**2) 
							
							dxidxc =-(blk(m0)%dxidx(1,ip,j,k)*temp00 + blk(m0)%dxidx(1,ip,jp,k)*temp10)*0.5d0    
							dxidyc =-(blk(m0)%dxidx(2,ip,j,k)*temp00 + blk(m0)%dxidx(2,ip,jp,k)*temp10)*0.5d0 
							dxidzc =-(blk(m0)%dxidx(3,ip,j,k)*temp00 + blk(m0)%dxidx(3,ip,jp,k)*temp10)*0.5d0 						         
							  
							miuc = (blk(m0)%amu(3,i,j,k) + blk(m0)%amu(3,i,jp,k))*0.5d0          
							s    =  blk(m0)%spcj(i,j,k)
							dn   =  0.5d0*(blk(m0)%spci(i,j,k)+blk(m0)%spci(i,jp,k))
							
							call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
							!!*
							fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
							fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
							fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
							!!*
						end do
					else if (face .eq. 2)then
						i  = ie
						im = ie-1
						do j = js,je-1
							jp = j+1
							
							uc = (blk(m0)%pri_v(2,im,j,k) + blk(m0)%pri_v(2,im,jp,k))*0.5d0
							vc = (blk(m0)%pri_v(3,im,j,k) + blk(m0)%pri_v(3,im,jp,k))*0.5d0
							wc = 0.d0
							pc = (blk(m0)%pri_v(5,i ,j,k) + blk(m0)%pri_v(5,i ,jp,k))*0.5d0
							
							temp00 = 1.d0/sqrt(blk(m0)%dxidx(1,im,j ,k)**2 + blk(m0)%dxidx(2,im,j ,k)**2 + blk(m0)%dxidx(3,im,j ,k)**2)  
							temp10 = 1.d0/sqrt(blk(m0)%dxidx(1,im,jp,k)**2 + blk(m0)%dxidx(2,im,jp,k)**2 + blk(m0)%dxidx(3,im,jp,k)**2) 
							
							dxidxc = (blk(m0)%dxidx(1,im,j,k)*temp00 + blk(m0)%dxidx(1,im,jp,k)*temp10)*0.5d0    
							dxidyc = (blk(m0)%dxidx(2,im,j,k)*temp00 + blk(m0)%dxidx(2,im,jp,k)*temp10)*0.5d0 
							dxidzc = (blk(m0)%dxidx(3,im,j,k)*temp00 + blk(m0)%dxidx(3,im,jp,k)*temp10)*0.5d0 						         
							
							miuc = (blk(m0)%amu(3,i,j,k) + blk(m0)%amu(3,i,jp,k))*0.5d0          
							s    =  blk(m0)%spcj(i,j,k)
							dn   =  0.5d0*(blk(m0)%spci(i,j,k)+blk(m0)%spci(i,jp,k))
							
							call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
							!!*
							fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
							fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
							fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
							!!*
						end do
					else if (face .eq. 3)then
						j  = js 
						jp = js+1
						do i = is,ie-1
							ip = i+1
							
							uc = (blk(m0)%pri_v(2,i,jp,k ) + blk(m0)%pri_v(2,ip,jp,k ))*0.5d0
							vc = (blk(m0)%pri_v(3,i,jp,k ) + blk(m0)%pri_v(3,ip,jp,k ))*0.5d0
							wc = 0.d0
							pc = (blk(m0)%pri_v(5,i,j ,k ) + blk(m0)%pri_v(5,ip,j ,k ))*0.5d0
							   
							temp00 = 1.d0/sqrt(blk(m0)%dxidx(4,i ,jp,k)**2 + blk(m0)%dxidx(5,i ,jp,k)**2 + blk(m0)%dxidx(6,i ,jp,k)**2)  
							temp10 = 1.d0/sqrt(blk(m0)%dxidx(4,ip,jp,k)**2 + blk(m0)%dxidx(5,ip,jp,k)**2 + blk(m0)%dxidx(6,ip,jp,k)**2)
							   
							dxidxc =-(blk(m0)%dxidx(4,i,jp,k)*temp00 + blk(m0)%dxidx(4,ip,jp,k)*temp10)*0.5d0    
							dxidyc =-(blk(m0)%dxidx(5,i,jp,k)*temp00 + blk(m0)%dxidx(5,ip,jp,k)*temp10)*0.5d0 
							dxidzc =-(blk(m0)%dxidx(6,i,jp,k)*temp00 + blk(m0)%dxidx(6,ip,jp,k)*temp10)*0.5d0 						         
							
							miuc = (blk(m0)%amu(3,i,j,k ) + blk(m0)%amu(3,ip,j,k ))*0.25d0          
							s    =  blk(m0)%spci(i,j,k)
							dn   = (blk(m0)%spcj(i,j,k)+blk(m0)%spcj(ip,j,k))*0.5d0
							        						     						
							call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
							!!*
							fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
							fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
							fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
							!!*
						end do
					else if (face .eq. 4)then
						j  = je
						jm = je-1
						do i = is,ie-1
							ip = i+1
							
							uc = (blk(m0)%pri_v(2,i,jm,k) + blk(m0)%pri_v(2,ip,jm,k))*0.5d0
							vc = (blk(m0)%pri_v(3,i,jm,k) + blk(m0)%pri_v(3,ip,jm,k))*0.5d0
							wc = 0.d0
							pc = (blk(m0)%pri_v(5,i,j ,k) + blk(m0)%pri_v(5,ip,j ,k))*0.5d0
							
							temp00 = 1.d0/sqrt(blk(m0)%dxidx(4,i ,jm,k)**2 + blk(m0)%dxidx(5,i ,jm,k)**2 + blk(m0)%dxidx(6,i ,jm,k)**2)  
							temp10 = 1.d0/sqrt(blk(m0)%dxidx(4,ip,jm,k)**2 + blk(m0)%dxidx(5,ip,jm,k)**2 + blk(m0)%dxidx(6,ip,jm,k)**2) 
							
							dxidxc = (blk(m0)%dxidx(4,i,jm,k)*temp00 + blk(m0)%dxidx(4,ip,jm,k)*temp10)*0.5d0    
							dxidyc = (blk(m0)%dxidx(5,i,jm,k)*temp00 + blk(m0)%dxidx(5,ip,jm,k)*temp10)*0.5d0 
							dxidzc = (blk(m0)%dxidx(6,i,jm,k)*temp00 + blk(m0)%dxidx(6,ip,jm,k)*temp10)*0.5d0 						         
							
							miuc = (blk(m0)%amu(3,i,j,k) + blk(m0)%amu(3,ip,j,k))*0.5d0          
							s    =  blk(m0)%spci(i,j,k)
							dn   = (blk(m0)%spcj(i,j,k)+blk(m0)%spcj(ip,j,k))*0.5d0
							
							call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
							!!*
							fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
							fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
							fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
							!!*
						end do
					end if
				end if
			end do
		end do		
	else if(iflag_dimension .eq. iflag_3d)then
		!!compute pressure force and viscous force
		!!o m0 = 1,blk_loop
		!!	do ksub = 1,blk(m0)%num_subface
		!!		if(blk(m0)%bc(ksub)%blk_t .eq. bc_wall)then
		!!			face = blk(m0)%bc(ksub)%face
		!!			!!wall dimention
		!!			is = blk(m0)%bc(ksub)%is
		!!			ie = blk(m0)%bc(ksub)%ie
		!!			js = blk(m0)%bc(ksub)%js
		!!			je = blk(m0)%bc(ksub)%je
		!!			ks = blk(m0)%bc(ksub)%ks
		!!			ke = blk(m0)%bc(ksub)%ke
		!!			!!*
		!!			if      (face .eq. 1)then
		!!				i  = is
		!!				ip = is+1
		!!				do k = ks,ke-1
		!!					do j = js,je-1
		!!						jp = j+1
		!!						kp = k+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,ip,j,k ) + blk(m0)%pri_v(2,ip,jp,k ) &
		!!					       +blk(m0)%pri_v(2,ip,j,kp) + blk(m0)%pri_v(2,ip,jp,kp))*0.25d0
		!!						vc = (blk(m0)%pri_v(3,ip,j,k ) + blk(m0)%pri_v(3,ip,jp,k ) &
		!!					       +blk(m0)%pri_v(3,ip,j,kp) + blk(m0)%pri_v(3,ip,jp,kp))*0.25d0
		!!						wc = (blk(m0)%pri_v(4,ip,j,k ) + blk(m0)%pri_v(4,ip,jp,k ) &
		!!					       +blk(m0)%pri_v(4,ip,j,kp) + blk(m0)%pri_v(4,ip,jp,kp))*0.25d0
		!!						pc = (blk(m0)%pri_v(5,i ,j,k ) + blk(m0)%pri_v(5,i ,jp,k ) &
		!!					       +blk(m0)%pri_v(5,i ,j,kp) + blk(m0)%pri_v(5,i ,jp,kp))*0.25d0
		!!					     
		!!						temp00 = 1.d0/sqrt(blk(m0)%dxidx2(1,ip,j ,k ))  
		!!						temp10 = 1.d0/sqrt(blk(m0)%dxidx2(1,ip,jp,k ))
		!!						temp01 = 1.d0/sqrt(blk(m0)%dxidx2(1,ip,j ,kp))
		!!						temp11 = 1.d0/sqrt(blk(m0)%dxidx2(1,ip,jp,kp)) 
		!!					     
		!!						dxidxc =-(blk(m0)%dxidx(1,ip,j,k )*temp00 + blk(m0)%dxidx(1,ip,jp,k )*temp10 &
		!!						         +blk(m0)%dxidx(2,ip,j,kp)*temp01 + blk(m0)%dxidx(2,ip,jp,kp)*temp11)*0.25d0    
		!!						dxidyc =-(blk(m0)%dxidx(3,ip,j,k )*temp00 + blk(m0)%dxidx(3,ip,jp,k )*temp10 &
		!!					             +blk(m0)%dxidx(1,ip,j,kp)*temp01 + blk(m0)%dxidx(1,ip,jp,kp)*temp11)*0.25d0 
		!!						dxidzc =-(blk(m0)%dxidx(2,ip,j,k )*temp00 + blk(m0)%dxidx(2,ip,jp,k )*temp10 &
		!!					             +blk(m0)%dxidx(3,ip,j,kp)*temp01 + blk(m0)%dxidx(3,ip,jp,kp)*temp11)*0.25d0 						         
		!!					
		!!						miuc = (blk(m0)%amu(3,i,j,k ) + blk(m0)%amu(3,i,jp,k ) &
		!!						       +blk(m0)%amu(3,i,j,kp) + blk(m0)%amu(3,i,jp,kp))*0.25d0          
		!!						s    =  blk(m0)%spcj(i,j,k)*blk(m0)%spck(i,j,k)
		!!						dn   =  blk(m0)%spci(i,j,k)
		!!					          						     						
		!!						call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!						!!*
		!!						fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!						fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!						fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!						!!*
		!!					end do
		!!				end do
		!!			else if (face .eq. 2)then
		!!				i  = ie
		!!				im = ie-1
		!!				do k = ks,ke-1
		!!					do j = js,je-1
		!!						jp = j+1
		!!						kp = k+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,im,j,k ) + blk(m0)%pri_v(2,im,jp,k ) &
		!!						     +blk(m0)%pri_v(2,im,j,kp) + blk(m0)%pri_v(2,im,jp,kp))*0.25d0
		!!						vc = (blk(m0)%pri_v(3,im,j,k ) + blk(m0)%pri_v(3,im,jp,k ) &
		!!						     +blk(m0)%pri_v(3,im,j,kp) + blk(m0)%pri_v(3,im,jp,kp))*0.25d0
		!!						wc = (blk(m0)%pri_v(4,im,j,k ) + blk(m0)%pri_v(4,im,jp,k ) &
		!!						     +blk(m0)%pri_v(4,im,j,kp) + blk(m0)%pri_v(4,im,jp,kp))*0.25d0
		!!						pc = (blk(m0)%pri_v(5,i ,j,k ) + blk(m0)%pri_v(5,i ,jp,k ) &
		!!						     +blk(m0)%pri_v(5,i ,j,kp) + blk(m0)%pri_v(5,i ,jp,kp))*0.25d0
		!!					     
		!!						temp00 = 1.d0/sqrt(blk(m0)%dxidx2(1,im,j ,k ))  
		!!						temp10 = 1.d0/sqrt(blk(m0)%dxidx2(1,im,jp,k ))
		!!						temp01 = 1.d0/sqrt(blk(m0)%dxidx2(1,im,j ,kp))
		!!						temp11 = 1.d0/sqrt(blk(m0)%dxidx2(1,im,jp,kp)) 
		!!					     
		!!						dxidxc = (blk(m0)%dxidx(1,im,j,k )*temp00 + blk(m0)%dxidx(1,im,jp,k )*temp10 &
		!!						         +blk(m0)%dxidx(2,im,j,kp)*temp01 + blk(m0)%dxidx(2,im,jp,kp)*temp11)*0.25d0    
		!!						dxidyc = (blk(m0)%dxidx(3,im,j,k )*temp00 + blk(m0)%dxidx(3,im,jp,k )*temp10 &
		!!						         +blk(m0)%dxidx(1,im,j,kp)*temp01 + blk(m0)%dxidx(1,im,jp,kp)*temp11)*0.25d0 
		!!						dxidzc = (blk(m0)%dxidx(2,im,j,k )*temp00 + blk(m0)%dxidx(2,im,jp,k )*temp10 &
		!!						         +blk(m0)%dxidx(3,im,j,kp)*temp01 + blk(m0)%dxidx(3,im,jp,kp)*temp11)*0.25d0 						         
		!!					
		!!						miuc = (blk(m0)%amu(3,i,j,k ) + blk(m0)%amu(3,i,jp,k ) &
		!!						       +blk(m0)%amu(3,i,j,kp) + blk(m0)%amu(3,i,jp,kp))*0.25d0          
		!!						s    =  blk(m0)%spcj(i,j,k)*blk(m0)%spck(i,j,k)
		!!						dn   =  blk(m0)%spci(i,j,k) 
		!!					          						     						
		!!						call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!						!!*
		!!						fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!						fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!						fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!						!!*
		!!					end do
		!!				end do
		!!			else if (face .eq. 3)then
		!!				j  = js 
		!!				jp = js+1
		!!				do k = ks,ke-1
		!!					do i = is,ie-1
		!!						ip = i+1
		!!						kp = k+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,i,jp,k ) + blk(m0)%pri_v(2,ip,jp,k ) &
		!!						     +blk(m0)%pri_v(2,i,jp,kp) + blk(m0)%pri_v(2,ip,jp,kp))*0.25d0
		!!						vc = (blk(m0)%pri_v(3,i,jp,k ) + blk(m0)%pri_v(3,ip,jp,k ) &
		!!						     +blk(m0)%pri_v(3,i,jp,kp) + blk(m0)%pri_v(3,ip,jp,kp))*0.25d0
		!!						wc = (blk(m0)%pri_v(4,i,jp,k ) + blk(m0)%pri_v(4,ip,jp,k ) &
		!!						     +blk(m0)%pri_v(4,i,jp,kp) + blk(m0)%pri_v(4,ip,jp,kp))*0.25d0
		!!						pc = (blk(m0)%pri_v(5,i,j ,k ) + blk(m0)%pri_v(5,ip,j ,k ) &
		!!						     +blk(m0)%pri_v(5,i,j ,kp) + blk(m0)%pri_v(5,ip,j ,kp))*0.25d0
		!!					     
		!!						temp00 = 1.d0/sqrt(blk(m0)%dxidx2(2,i ,jp,k ))  
		!!						temp10 = 1.d0/sqrt(blk(m0)%dxidx2(2,ip,jp,k ))
		!!						temp01 = 1.d0/sqrt(blk(m0)%dxidx2(2,i ,jp,kp))
		!!						temp11 = 1.d0/sqrt(blk(m0)%dxidx2(2,ip,jp,kp)) 
		!!					     
		!!						dxidxc =-(blk(m0)%dxidx(4,i,jp,k )*temp00 + blk(m0)%dxidx(4,ip,jp,k )*temp10 &
		!!					             +blk(m0)%dxidx(5,i,jp,kp)*temp01 + blk(m0)%dxidx(5,ip,jp,kp)*temp11)*0.25d0    
		!!						dxidyc =-(blk(m0)%dxidx(6,i,jp,k )*temp00 + blk(m0)%dxidx(6,ip,jp,k )*temp10 &
		!!					             +blk(m0)%dxidx(4,i,jp,kp)*temp01 + blk(m0)%dxidx(4,ip,jp,kp)*temp11)*0.25d0 
		!!						dxidzc =-(blk(m0)%dxidx(5,i,jp,k )*temp00 + blk(m0)%dxidx(5,ip,jp,k )*temp10 &
		!!					             +blk(m0)%dxidx(6,i,jp,kp)*temp01 + blk(m0)%dxidx(6,ip,jp,kp)*temp11)*0.25d0 						         
		!!					
		!!					  miuc = (blk(m0)%amu(3,i,j,k ) + blk(m0)%amu(3,ip,j,k ) &
		!!					         +blk(m0)%amu(3,i,j,kp) + blk(m0)%amu(3,ip,j,kp))*0.25d0          
		!!					  s    =  blk(m0)%spci(i,j,k)*blk(m0)%spck(i,j,k)
		!!					  dn   =  blk(m0)%spcj(i,j,k) 
		!!					          						     						
		!!					  call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!					  !!*
		!!					  fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!					  fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!					  fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!					  !!*
		!!					end do
		!!				end do
		!!			else if (face .eq. 4)then
		!!				j  = je
		!!				jm = je-1
		!!				do k = ks,ke-1
		!!					do i = is,ie-1
		!!						ip = i+1
		!!						kp = k+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,i,jm,k ) + blk(m0)%pri_v(2,ip,jm,k ) &
		!!					       +blk(m0)%pri_v(2,i,jm,kp) + blk(m0)%pri_v(2,ip,jm,kp))*0.25d0
		!!					  vc = (blk(m0)%pri_v(3,i,jm,k ) + blk(m0)%pri_v(3,ip,jm,k ) &
		!!					       +blk(m0)%pri_v(3,i,jm,kp) + blk(m0)%pri_v(3,ip,jm,kp))*0.25d0
		!!					  wc = (blk(m0)%pri_v(4,i,jm,k ) + blk(m0)%pri_v(4,ip,jm,k ) &
		!!					       +blk(m0)%pri_v(4,i,jm,kp) + blk(m0)%pri_v(4,ip,jm,kp))*0.25d0
		!!					  pc = (blk(m0)%pri_v(5,i,j ,k ) + blk(m0)%pri_v(5,ip,j ,k ) &
		!!					       +blk(m0)%pri_v(5,i,j ,kp) + blk(m0)%pri_v(5,ip,j ,kp))*0.25d0
		!!					     
		!!					  temp00 = 1.d0/sqrt(blk(m0)%dxidx2(2,i ,jm,k ))  
		!!					  temp10 = 1.d0/sqrt(blk(m0)%dxidx2(2,ip,jm,k ))
		!!					  temp01 = 1.d0/sqrt(blk(m0)%dxidx2(2,i ,jm,kp))
		!!					  temp11 = 1.d0/sqrt(blk(m0)%dxidx2(2,ip,jm,kp)) 
		!!					     
		!!					  dxidxc = (blk(m0)%dxidx(4,i,jm,k )*temp00 + blk(m0)%dxidx(4,ip,jm,k )*temp10 &
		!!					           +blk(m0)%dxidx(5,i,jm,kp)*temp01 + blk(m0)%dxidx(5,ip,jm,kp)*temp11)*0.25d0    
		!!					  dxidyc = (blk(m0)%dxidx(6,i,jm,k )*temp00 + blk(m0)%dxidx(6,ip,jm,k )*temp10 &
		!!					           +blk(m0)%dxidx(4,i,jm,kp)*temp01 + blk(m0)%dxidx(4,ip,jm,kp)*temp11)*0.25d0 
		!!					  dxidzc = (blk(m0)%dxidx(5,i,jm,k )*temp00 + blk(m0)%dxidx(5,ip,jm,k )*temp10 &
		!!					           +blk(m0)%dxidx(6,i,jm,kp)*temp01 + blk(m0)%dxidx(6,ip,jm,kp)*temp11)*0.25d0 						         
		!!					
		!!					  miuc = (blk(m0)%amu(3,i,j,k ) + blk(m0)%amu(3,ip,j,k ) &
		!!					         +blk(m0)%amu(3,i,j,kp) + blk(m0)%amu(3,ip,j,kp))*0.25d0          
		!!					  s    =  blk(m0)%spci(i,j,k)*blk(m0)%spck(i,j,k)
		!!					  dn   =  blk(m0)%spcj(i,j,k) 
		!!					          						     						
		!!					  call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!					  !!*
		!!					  fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!					  fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!					  fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!					  !!*
		!!					end do
		!!				end do
		!!			else if (face .eq. 5)then
		!!				k  = ks
		!!				kp = ks+1
		!!				do j = js,je-1
		!!					do i = is,ie-1
		!!						ip = i+1
		!!						jp = j+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,i,j ,kp) + blk(m0)%pri_v(2,ip,j ,kp) &
		!!					       +blk(m0)%pri_v(2,i,jp,kp) + blk(m0)%pri_v(2,ip,jp,kp))*0.25d0
		!!					    vc = (blk(m0)%pri_v(3,i,j ,kp) + blk(m0)%pri_v(3,ip,j ,kp) &
		!!					       +blk(m0)%pri_v(3,i,jp,kp) + blk(m0)%pri_v(3,ip,jp,kp))*0.25d0
		!!					    wc = (blk(m0)%pri_v(4,i,j ,kp) + blk(m0)%pri_v(4,ip,j ,kp) &
		!!					       +blk(m0)%pri_v(4,i,jp,kp) + blk(m0)%pri_v(4,ip,jp,kp))*0.25d0
		!!					    pc = (blk(m0)%pri_v(5,i,j ,k ) + blk(m0)%pri_v(5,ip,j ,k ) &
		!!					       +blk(m0)%pri_v(5,i,jp,k ) + blk(m0)%pri_v(5,ip,jp,k ))*0.25d0
		!!					     
		!!					  temp00 = 1.d0/sqrt(blk(m0)%dxidx2(3,i ,j ,kp))  
		!!					  temp10 = 1.d0/sqrt(blk(m0)%dxidx2(3,ip,j ,kp))
		!!					  temp01 = 1.d0/sqrt(blk(m0)%dxidx2(3,i ,jp,kp))
		!!					  temp11 = 1.d0/sqrt(blk(m0)%dxidx2(3,ip,jp,kp)) 
		!!					     
		!!					  dxidxc =-(blk(m0)%dxidx(7,i,j ,kp)*temp00 + blk(m0)%dxidx(7,ip,j ,kp)*temp10 &
		!!					           +blk(m0)%dxidx(8,i,jp,kp)*temp01 + blk(m0)%dxidx(8,ip,jp,kp)*temp11)*0.25d0    
		!!					  dxidyc =-(blk(m0)%dxidx(9,i,j ,kp)*temp00 + blk(m0)%dxidx(9,ip,j ,kp)*temp10 &
		!!					           +blk(m0)%dxidx(7,i,jp,kp)*temp01 + blk(m0)%dxidx(7,ip,jp,kp)*temp11)*0.25d0 
		!!					  dxidzc =-(blk(m0)%dxidx(8,i,j ,kp)*temp00 + blk(m0)%dxidx(8,ip,j ,kp)*temp10 &
		!!					           +blk(m0)%dxidx(9,i,jp,kp)*temp01 + blk(m0)%dxidx(9,ip,jp,kp)*temp11)*0.25d0 						         
		!!					
		!!					  miuc = (blk(m0)%amu(3,i,j ,k) + blk(m0)%amu(3,ip,j ,k) &
		!!					         +blk(m0)%amu(3,i,jp,k) + blk(m0)%amu(3,ip,jp,k))*0.25d0          
		!!					  s    =  blk(m0)%spci(i,j,k)*blk(m0)%spcj(i,j,k)
		!!					  dn   =  blk(m0)%spck(i,j,k) 
		!!					          						     						
		!!					  call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!					  !!*
		!!					  fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!					  fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!					  fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!					  !!*
		!!					end do
		!!				end do
		!!			else if (face .eq. 6)then
		!!				k  = ke     
		!!				km = ke-1
		!!				do j = js,je-1
		!!					do i = is,ie-1
		!!						ip = i+1
		!!						jp = j+1
		!!						
		!!						uc = (blk(m0)%pri_v(2,i,j ,km) + blk(m0)%pri_v(2,ip,j ,km) &
		!!					       +blk(m0)%pri_v(2,i,jp,km) + blk(m0)%pri_v(2,ip,jp,km))*0.25d0
		!!					  vc = (blk(m0)%pri_v(3,i,j ,km) + blk(m0)%pri_v(3,ip,j ,km) &
		!!					       +blk(m0)%pri_v(3,i,jp,km) + blk(m0)%pri_v(3,ip,jp,km))*0.25d0
		!!					  wc = (blk(m0)%pri_v(4,i,j ,km) + blk(m0)%pri_v(4,ip,j ,km) &
		!!					       +blk(m0)%pri_v(4,i,jp,km) + blk(m0)%pri_v(4,ip,jp,km))*0.25d0
		!!					  pc = (blk(m0)%pri_v(5,i,j ,k ) + blk(m0)%pri_v(5,ip,j ,k ) &
		!!					       +blk(m0)%pri_v(5,i,jp,k ) + blk(m0)%pri_v(5,ip,jp,k ))*0.25d0
		!!					     
		!!					  temp00 = 1.d0/sqrt(blk(m0)%dxidx2(3,i ,j ,km))  
		!!					  temp10 = 1.d0/sqrt(blk(m0)%dxidx2(3,ip,j ,km))
		!!					  temp01 = 1.d0/sqrt(blk(m0)%dxidx2(3,i ,jp,km))
		!!					  temp11 = 1.d0/sqrt(blk(m0)%dxidx2(3,ip,jp,km)) 
		!!					     
		!!					  dxidxc = (blk(m0)%dxidx(7,i,j ,km)*temp00 + blk(m0)%dxidx(7,ip,j ,km)*temp10 &
		!!					           +blk(m0)%dxidx(8,i,jp,km)*temp01 + blk(m0)%dxidx(8,ip,jp,km)*temp11)*0.25d0    
		!!					  dxidyc = (blk(m0)%dxidx(9,i,j ,km)*temp00 + blk(m0)%dxidx(9,ip,j ,km)*temp10 &
		!!					           +blk(m0)%dxidx(7,i,jp,km)*temp01 + blk(m0)%dxidx(7,ip,jp,km)*temp11)*0.25d0 
		!!					  dxidzc = (blk(m0)%dxidx(8,i,j ,km)*temp00 + blk(m0)%dxidx(8,ip,j ,km)*temp10 &
		!!					           +blk(m0)%dxidx(9,i,jp,km)*temp01 + blk(m0)%dxidx(9,ip,jp,km)*temp11)*0.25d0 						         
		!!					
		!!					  miuc = (blk(m0)%amu(3,i,j ,k) + blk(m0)%amu(3,ip,j ,k) &
		!!					         +blk(m0)%amu(3,i,jp,k) + blk(m0)%amu(3,ip,jp,k))*0.25d0          
		!!					  s    =  blk(m0)%spci(i,j,k)*blk(m0)%spcj(i,j,k)
		!!					  dn   =  blk(m0)%spck(i,j,k) 
		!!					          						     						
		!!					  call force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
		!!					  !!*
		!!					  fp(1) = fp(1) + f_prs(1); fv(1) = fv(1) + f_vis(1)
		!!					  fp(2) = fp(2) + f_prs(2); fv(2) = fv(2) + f_vis(2)
		!!					  fp(3) = fp(3) + f_prs(3); fv(3) = fv(3) + f_vis(3)
		!!					  !!*
		!!					end do
		!!				end do
		!!			end if
		!!		end if
		!!	end do
		!!end do
	end if
	
	call MPI_Reduce(fp,ffp,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_Reduce(fv,ffv,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	
	if(myid .eq. root)then
		print *, 'force:'
		write(*, 5) ffp(1),ffp(2),ffp(3)
		write(*,10) ffv(1),ffv(2),ffv(3)
5   format("fpx: ",e10.4," fpy: ",e10.4," fpz: ",e10.4)
10  format("fvx: ",e10.4," fvy: ",e10.4," fvz: ",e10.4)
		!!get lift & drag coefficients
		cl = ((ffp(2)+ffv(2))*cos(aoa) - (ffp(1)+ffv(1))*sin(aoa))*div
		cd = ((ffp(2)+ffv(2))*sin(aoa) + (ffp(1)+ffv(1))*cos(aoa))*div
		print *, "cl = ", cl, "cd = ", cd
		!!*
		
		!!output to log file
		write(force_out,"('result/clcd_history.dat')")
		open(99, file = force_out, status = 'old', position = 'append')
		write(99,*) timestep*dt,cl,cd
		close(99)
		!!*
	end if	
	
	return
end subroutine clcd_output

subroutine force_compute(f_prs,f_vis,uc,vc,wc,pc,dxidxc,dxidyc,dxidzc,s,dn,pinf,miuc,re)
	implicit none
	!!*******************************************************************
	real*8 :: f_prs(3),f_vis(3)           !!1: x-direction,2: y-direction,3: z-direction
	real*8 :: s,dn                        !!s: element area, dn: wall dist
	real*8 :: pinf                        !!reference pressure
	real*8 :: re                          !!reynolds number
	!!*******************************************************************
	real*8 :: uc,vc,wc,pc
	real*8 :: miuc
	real*8 :: dxidxc,dxidyc,dxidzc
	real*8 :: un
	!!*******************************************************************
	!!pressure force
	f_prs(1) = (pc - pinf)*s*dxidxc
	f_prs(2) = (pc - pinf)*s*dxidyc
	f_prs(3) = (pc - pinf)*s*dxidzc
	!!*
	
	!!viscous force
	un   = uc*dxidxc + vc*dxidyc + wc*dxidzc
	
	f_vis(1) = miuc*(uc - un*dxidxc)/dn*s/re
	f_vis(2) = miuc*(vc - un*dxidyc)/dn*s/re
	f_vis(3) = miuc*(wc - un*dxidzc)/dn*s/re
	!!*
	
	return
end subroutine
