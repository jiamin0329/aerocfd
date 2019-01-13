!!========================================================= 
!!	initialize the fluid field
!!
!!	purpose: initialize the flow field in two methods:
!!           1. with an initial flow field
!!           2. with a restart file
!!                                                              
!!	author: jiamin xu                                    
!!	date:   2012.11.22                      
!!========================================================= 
subroutine init_flowfield
	use blk_var
	use ns_const
	use mpi_var
	use flag_var
	use index_var
	!!use sa_var
	implicit none
	
	integer :: i,j,k
	real*8  :: q2
	real*8,dimension(:,:,:),allocatable :: tempx,tempy,tempz
	real*8,dimension(:,:,:),allocatable :: tempd,tempu,tempv,tempw,tempt,tempp
	real*8,dimension(:,:,:),allocatable :: temptur1,temptur2,tempamut
	
	character(len = 180) :: restart
	character(len = 180) :: temp1,temp2
	integer:: is,ie,js,je,ks,ke
	!!non-dimensionalization
	d_0   = 1.d0
	u_0   = 1.d0
	t_0   = 1.d0
	p_0   = d_0*t_0/gamma/ma/ma
	c_0   = sqrt(t_0)/ma
	
	cv = 1.d0/gamma/(gamma-1.d0)/ma/ma
	cp = 1.d0/(gamma-1.d0)/ma/ma

	aoa  = aoa*pi/180.0d0
	!!*	

	!!initialized from initial flow field
	if (iflag_init .eq. 0) then
		if(myid .eq. 0)then  
			write(*,*)'***********************************************************'
			write(*,*)'initialize the flowfield from initial value......'
			write(*,*)'***********************************************************'
		end if
		
		do m0 = 1,blk_loop
			is = blk(m0)%is; ie = blk(m0)%ie
			js = blk(m0)%js; je = blk(m0)%je
			ks = blk(m0)%ks; ke = blk(m0)%ke
			
			blk(m0)%dudx    = 0.d0
			blk(m0)%vor     = 0.d0
			blk(m0)%tao     = 0.d0
			blk(m0)%rhsv    = 0.d0
			blk(m0)%rhsi    = 0.d0
			blk(m0)%rhs     = 0.d0
			!!blk(m0)%e       = 0.d0
			!!blk(m0)%f       = 0.d0
			!!blk(m0)%g       = 0.d0
			!!blk(m0)%dedxi   = 0.d0
			!!blk(m0)%dfdeta  = 0.d0
			!!blk(m0)%dgdzeta = 0.d0
		
			do k = ks,ke
			do j = js,je
			do i = is,ie
				!!conservative variables
				blk(m0)%q(1,i,j,k) = d_0
				blk(m0)%q(2,i,j,k) = d_0*u_0*cos(aoa)
				blk(m0)%q(3,i,j,k) = d_0*u_0*sin(aoa)
				blk(m0)%q(4,i,j,k) = 0.d0
				blk(m0)%q(5,i,j,k) = d_0*t_0*cv + 0.5d0*d_0*u_0*u_0
                
				!!blk(m0)%q0(:,1,i,j,k) = blk(m0)%q(1,i,j,k)
				!!blk(m0)%q0(:,2,i,j,k) = blk(m0)%q(2,i,j,k)
				!!blk(m0)%q0(:,3,i,j,k) = blk(m0)%q(3,i,j,k)
				!!blk(m0)%q0(:,4,i,j,k) = blk(m0)%q(4,i,j,k)
				!!blk(m0)%q0(:,5,i,j,k) = blk(m0)%q(5,i,j,k)
						
				!!primitive variables
				blk(m0)%pri_v(1,i,j,k) = d_0
				blk(m0)%pri_v(2,i,j,k) = u_0*cos(aoa)
				blk(m0)%pri_v(3,i,j,k) = u_0*sin(aoa)
				blk(m0)%pri_v(4,i,j,k) = 0.d0
				blk(m0)%pri_v(5,i,j,k) = p_0
				blk(m0)%pri_v(6,i,j,k) = t_0
				blk(m0)%pri_v(7,i,j,k) = c_0
				!!
				!!blk(m0)%nut  (  i,j,k) = 0.d0

						
				!!local time step
				blk(m0)%dt(i,j,k) = dt
				!!initialize the turbulence model variables
				!!if (iflag_turbulence .eq. iflag_sa) then
				!!	blk(m0)%nut      (i,j,k) = nutinf
				!!	blk(m0)%amu    (2,i,j,k) = 0.009d0
				!!	blk(m0)%sa_rhs   (i,j,k) = 0.d0
				!!end if
				
				if (iflag_turbulence .eq. iflag_les .or. iflag_turbulence .eq. iflag_laminar) then
					blk(m0)%amu(2,i,j,k) = 0.d0
				end if
				
				blk(m0)%amu(1,i,j,k) = 1.d0
				blk(m0)%amu(3,i,j,k) = blk(m0)%amu(1,i,j,k) + blk(m0)%amu(2,i,j,k)
			end do
			end do
			end do
		end do
	end if
	!!*
	
	!!initialized from the restart file
	if (iflag_init .eq. 1) then
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		if(myid .eq. 0)then  
			print *,'***********************************************************'
			print *,'initialize the flowfield from restartfile......'
			print *,'***********************************************************'
		end if
				
		!write(restart,"('result/restart_'I4.4'.dat')"),myid
		write(restart,"('result/flowfield_'I10.10'_'I4.4'.dat')"),restart_timestep,myid
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		open (90,file = restart)
		do m0 = 1,blk_loop
			is = 1; ie = blk(m0)%ni
			js = 1; je = blk(m0)%nj
			ks = 1; ke = blk(m0)%nk
			
			blk(m0)%dudx    = 0.d0 
			blk(m0)%vor     = 0.d0 
			blk(m0)%tao     = 0.d0 
			blk(m0)%rhsv    = 0.d0 
			blk(m0)%rhsi    = 0.d0 
			blk(m0)%rhs     = 0.d0 
			!!blk(m0)%e       = 0.d0 
			!!blk(m0)%f       = 0.d0 
			!!blk(m0)%g       = 0.d0 
			!!blk(m0)%dedxi   = 0.d0 
			!!blk(m0)%dfdeta  = 0.d0 
			!!blk(m0)%dgdzeta = 0.d0 
			
			allocate(tempx   (is:ie,js:je,ks:ke))
			allocate(tempy   (is:ie,js:je,ks:ke))
			allocate(tempz   (is:ie,js:je,ks:ke))
			allocate(tempd   (is:ie,js:je,ks:ke))
			allocate(tempu   (is:ie,js:je,ks:ke))
			allocate(tempv   (is:ie,js:je,ks:ke))
			allocate(tempw   (is:ie,js:je,ks:ke))
			allocate(tempp   (is:ie,js:je,ks:ke))
			allocate(tempt   (is:ie,js:je,ks:ke))
			allocate(temptur1(is:ie,js:je,ks:ke))
			allocate(temptur2(is:ie,js:je,ks:ke))
			
			read(90,*)
			read(90,*)
			
			if      (iflag_dimension .eq. iflag_2d) then
				k = 1
				do j = js,je
				do i = is,ie
					read(90,*) tempx(i,j,k),tempy(i,j,k),                           &
						       tempd(i,j,k),tempu(i,j,k),tempv(i,j,k),tempp(i,j,k), &
						       temptur1(i,j,k)
				end do
				end do
				
				do j = js,je
				do i = is,ie
					tempt(i,j,k)       = tempp(i,j,k)/tempd(i,j,k)*gamma*ma*ma
					blk(m0)%q(1,i,j,k) = tempd(i,j,k)
					blk(m0)%q(2,i,j,k) = tempd(i,j,k)*tempu(i,j,k)
					blk(m0)%q(3,i,j,k) = tempd(i,j,k)*tempv(i,j,k)
					blk(m0)%q(4,i,j,k) = tempd(i,j,k)*tempw(i,j,k)
					blk(m0)%q(5,i,j,k) = tempd(i,j,k)*tempt(i,j,k)*cv                   &
							           + 0.5d0*tempd(i,j,k)*(tempu(i,j,k)*tempu(i,j,k)  &
							                               + tempv(i,j,k)*tempv(i,j,k)  &
							                               + tempw(i,j,k)*tempw(i,j,k))
					!!blk(m0)%q0(:,1,i,j,k) = blk(m0)%q(1,i,j,k)
					!!blk(m0)%q0(:,2,i,j,k) = blk(m0)%q(2,i,j,k)
					!!blk(m0)%q0(:,3,i,j,k) = blk(m0)%q(3,i,j,k)
					!!blk(m0)%q0(:,4,i,j,k) = blk(m0)%q(4,i,j,k)
					!!blk(m0)%q0(:,5,i,j,k) = blk(m0)%q(5,i,j,k)
						
					blk(m0)%pri_v(1,i,j,k) = blk(m0)%q(1,i,j,k)
					blk(m0)%pri_v(2,i,j,k) = blk(m0)%q(2,i,j,k)/blk(m0)%pri_v(1,i,j,k)
					blk(m0)%pri_v(3,i,j,k) = blk(m0)%q(3,i,j,k)/blk(m0)%pri_v(1,i,j,k)
					blk(m0)%pri_v(4,i,j,k) = blk(m0)%q(4,i,j,k)/blk(m0)%pri_v(1,i,j,k)
					q2 = blk(m0)%pri_v(1,i,j,k)**2 + blk(m0)%pri_v(2,i,j,k)**2 + blk(m0)%pri_v(3,i,j,k)**2 
					blk(m0)%pri_v(6,i,j,k) =(blk(m0)%q(5,i,j,k)/blk(m0)%pri_v(1,i,j,k) - 0.5d0*q2)/cv
					blk(m0)%pri_v(5,i,j,k) = blk(m0)%pri_v(1,i,j,k)*blk(m0)%pri_v(6,i,j,k)/gamma/ma/ma
					blk(m0)%pri_v(7,i,j,k) = sqrt(blk(m0)%pri_v(5,i,j,k))/ma
						
					blk(m0)%rhs (:,i,j,k) = 0.d0
					blk(m0)%rhsi(:,i,j,k) = 0.d0
					blk(m0)%rhsv(:,i,j,k) = 0.d0
					blk(m0)%dt    (i,j,k) = dt
					
					blk(m0)%amu (1,i,j,k) = 1.d0
						
					!!if      (iflag_turbulence .eq. iflag_sa          ) then
					!!	blk(m0)%nut(i,j,k) = temptur1(i,j,k)
					!!else if (iflag_turbulence .eq. iflag_kwsst_menter) then
					!!end if
				end do
				end do	
			else if (iflag_dimension .eq. iflag_3d) then
				do k = ks,ke
					do j = js,je
						do i = is,ie
							read(90,*) tempx(i,j,k),tempy(i,j,k),tempz(i,j,k),                           &
							           tempd(i,j,k),tempu(i,j,k),tempv(i,j,k),tempw(i,j,k),tempp(i,j,k), &
							           temptur1(i,j,k)
						end do
					end do
				end do
				
				
				do k = ks,ke
					do j = js,je
						do i = is,ie
							tempt(i,j,k)       = tempp(i,j,k)/tempd(i,j,k)*gamma*ma*ma
							blk(m0)%q(1,i,j,k) = tempd(i,j,k)
							blk(m0)%q(2,i,j,k) = tempd(i,j,k)*tempu(i,j,k)
							blk(m0)%q(3,i,j,k) = tempd(i,j,k)*tempv(i,j,k)
							blk(m0)%q(4,i,j,k) = tempd(i,j,k)*tempw(i,j,k)
							blk(m0)%q(5,i,j,k) = tempd(i,j,k)*tempt(i,j,k)*cv                   &
						  	                 + 0.5d0*tempd(i,j,k)*(tempu(i,j,k)*tempu(i,j,k)  &
						  	                                     + tempv(i,j,k)*tempv(i,j,k)  &
						  	                                     + tempw(i,j,k)*tempw(i,j,k))
							!!blk(m0)%q0(:,1,i,j,k) = blk(m0)%q(1,i,j,k)
							!!blk(m0)%q0(:,2,i,j,k) = blk(m0)%q(2,i,j,k)
							!!blk(m0)%q0(:,3,i,j,k) = blk(m0)%q(3,i,j,k)
							!!blk(m0)%q0(:,4,i,j,k) = blk(m0)%q(4,i,j,k)
							!!blk(m0)%q0(:,5,i,j,k) = blk(m0)%q(5,i,j,k)
							
							blk(m0)%pri_v(1,i,j,k) = blk(m0)%q(1,i,j,k)
							blk(m0)%pri_v(2,i,j,k) = blk(m0)%q(2,i,j,k)/blk(m0)%pri_v(1,i,j,k)
							blk(m0)%pri_v(3,i,j,k) = blk(m0)%q(3,i,j,k)/blk(m0)%pri_v(1,i,j,k)
							blk(m0)%pri_v(4,i,j,k) = blk(m0)%q(4,i,j,k)/blk(m0)%pri_v(1,i,j,k)
							q2 = blk(m0)%pri_v(1,i,j,k)**2 + blk(m0)%pri_v(2,i,j,k)**2 + blk(m0)%pri_v(3,i,j,k)**2 
							blk(m0)%pri_v(6,i,j,k) =(blk(m0)%q(5,i,j,k)/blk(m0)%pri_v(1,i,j,k) - 0.5d0*q2)/cv
							blk(m0)%pri_v(5,i,j,k) = blk(m0)%pri_v(1,i,j,k)*blk(m0)%pri_v(6,i,j,k)/gamma/ma/ma
							blk(m0)%pri_v(7,i,j,k) = sqrt(blk(m0)%pri_v(5,i,j,k))/ma
							
							blk(m0)%rhs (:,i,j,k) = 0.d0
							blk(m0)%rhsi(:,i,j,k) = 0.d0
							blk(m0)%rhsv(:,i,j,k) = 0.d0
							blk(m0)%dt    (i,j,k) = dt
							
							blk(m0)%amu (1,i,j,k) = 1.d0
							
							!!if      (iflag_turbulence .eq. iflag_sa          ) then
							!!	blk(m0)%nut(i,j,k) = temptur1(i,j,k)
							!!else if (iflag_turbulence .eq. iflag_kwsst_menter) then
                            !!
							!!end if
							
						end do
					end do
				end do			
			end if

			deallocate (tempx,tempy,tempz,tempd,tempu,tempv,tempw,tempt,tempp,temptur1,temptur2)
		end do
		close(90)
	end if
	!!*end initialization of the flow field
	
	if(myid .eq. 0)then  
		write(*,*)"***********************************************************"
		write(*,*)"finish initialization!"
		write(*,*)"***********************************************************"
		write(*,*)"Press ENTER to continue......"
	end if
	!!*
	
	return
end subroutine

