!!==================================================================!!
!! This is a 2d/3d multi-block structured compressible flow solver	!!
!! finite difference method                                         !!
!! MPI version                                                      !!
!!                                                                  !!
!! Version 2.0 (with overlap block interface implemented)          !!
!! Author: Jiamin Xu                                                !!
!!==================================================================!!
program sjtucfd_mpi
	use blk_var                  !!main data structure: blk,bc and cp
	use flag_var                 !!flags,global variables
	use index_var                !!index variables
	use mpi_var                  !!variables for mpi
	use output_var               !!variables for output
	use ns_const                 !!gas constants in ns equation
	use glbindex_var             !!global index
	implicit none
	
	integer :: i,j,k                   !!index of block dimension
	integer :: sub                     !!index of subiteration in rk or implicit timeadvancement
	integer :: outinterval             !!interval of screen-output
	integer :: n                       !!index of processors
	integer :: is1,ie1,js1,je1,ks1,ke1 !!block dimension
	integer :: is0,ie0,js0,je0,ks0,ke0 !!block dimension
	integer :: is, ie, js, je, ks, ke  !!block dimension

	character(len = 180) :: debugFile
	isDebug = 0

	!!get mpi info
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myid,     ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	call MPI_GET_PROCESSOR_NAME(processor_name,namelen,ierr)
	!!*end get mpi info

	!!print info page
	if (myid .eq. root) then
		write (*,*) '***********************************************************'
		write (*,*) '*  This is a 2d/3d compressible flow solver (mpi version) *'
		write (*,*) '*                                                         *'
		write (*,*) '*                                                         *'
		write (*,*) '*                                                         *'
		write (*,*) '*              Computational AeroAcoustics Research Group *'
		write (*,*) '*                  School of Aeronautics and Astronautics *'
		write (*,*) '*                            Shanghai JiaoTong University *'
		write (*,*) '*                 copyright reserved by prof. Song Wenbin *'
		write (*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*end print info page
	
	!!checking threads
	if(myid .eq. root)then
		write (*,*) "***********************************************************"
		write (*,*) " Checking Threads...... "
		write (*,*) "***********************************************************"
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	do n = 1,numprocs
		if (myid .eq. n-1) then
			write(*,5) myid,numprocs,processor_name
5 		format(" Hey Jude~ from ",I4," of ",I4," on ",A10)
		end if
	end do
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	if(myid .eq. root)then
		write (*,*) "***********************************************************"
		write (*,*) " Finish Checking Threads! "
		write (*,*) "***********************************************************"
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*end checking threads
	
	!!***************************************************!!
	!!                 initialization                    !!
	!!***************************************************!!
	call init_readcontrolfile
	call GetBufferLength(bufferLength, iflag_inviscid)

	if (isDebug .eq. 1) then
		write (*,*) "Debug info: buffer length: ", bufferLength  
	end if

	if (iflag_turbulence .ne. iflag_laminar .and. iflag_turbulence .ne. iflag_sa) then 
		print *, "Wrong turbulence model type is input!!!"
		stop
	end if

	!!read total number of blocks and get myn/lastn
	if(myid .eq. root)then
		write (*,*)"***********************************************************"
		write (*,*)"reading block dimension......                              "
		write (*,*)"***********************************************************"
		write (*,*)"block info:"
		open(99,file = grid)
		read(99,*) nblock
		write (*,*)"total number of blocks = ",nblock
		close(99)
		!!get total block number in the current processor
		if      (mod(nblock,numprocs) .eq. 0) then
			myn   = nblock/numprocs
			lastn = myn
		else if (mod(nblock,numprocs) .ne. 0) then
			myn   = nblock/numprocs+1
			lastn = nblock - myn*(numprocs-1)
			
			if (lastn .le. 0) then
				myn   = nblock/numprocs
				lastn = nblock - myn*(numprocs-1)
			end if
		end if

		if (lastn .le. 0) then
			print *, "Wrong block distribution!"
			stop
		end if
		!!*
	end if

	call MPI_BCAST(nblock,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(myn,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lastn, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	if      (myid .lt. numprocs-1)then
		blk_loop = myn
	else if (myid .eq. numprocs-1)then
		blk_loop = lastn
	end if
	!!*
	
	!!allocate for mod
	allocate (blk (blk_loop)) !!in current processor
	allocate (blk0(nblock  )) !!global index variables in every processor
	
	call init_readbcinfo
	call init_allocatememory
	call init_readbcinfo2
	call UpdateBufferCoordinate
	call init_flowfield

	if(myid .eq. root)then
		call CheckBoundary
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!!grid spacing calculation
	if(myid .eq. root)then
		write(*,*) '***********************************************************'
		write(*,*) 'Computing grid space and wall distance......               '
		write(*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	do m0 = 1,blk_loop
		is = blk(m0)%is;ie = blk(m0)%ie
		js = blk(m0)%js;je = blk(m0)%je
		ks = blk(m0)%ks;ke = blk(m0)%ke

		is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
		js0 = blk(m0)%js0;je0 = blk(m0)%je0
		ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0
		
		call init_spacing (blk(m0)%spci,blk(m0)%spcj,blk(m0)%spck,blk(m0)%spacing,blk(m0)%vol, &
		                   blk(m0)%x,blk(m0)%y,blk(m0)%z,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	end do
	if(myid .eq. root)then
		write(*,*) 'Finish computing grid space!'
	end if

	if(iflag_solver .eq. iflag_nssolver .and. iflag_turbulence .gt. 0) then
		call init_walldist
	end if

	if(myid .eq. root)then
		write(*,*) 'Finish computing walldistance!'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
	!!*
	
	!!jacobian transformation
	if(myid .eq. root)then
		write(*,*) '***********************************************************'
		write(*,*) 'Computing jacobian matics......                            '
		write(*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	do m0 = 1,blk_loop
		is = blk(m0)%is;ie = blk(m0)%ie
		js = blk(m0)%js;je = blk(m0)%je
		ks = blk(m0)%ks;ke = blk(m0)%ke

		is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
		js0 = blk(m0)%js0;je0 = blk(m0)%je0
		ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0
		
		call gridtransformation(blk(m0)%dxidx,blk(m0)%inv_j, &
		                        blk(m0)%aalpha,blk(m0)%bbeta,blk(m0)%ggamma, &
		                        blk(m0)%x,blk(m0)%y,blk(m0)%z,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,m0)
	end do

	if (isDebug .eq. 1) then
	do m0 = 1,blk_loop
		is = blk(m0)%is;ie = blk(m0)%ie
		js = blk(m0)%js;je = blk(m0)%je
		ks = blk(m0)%ks;ke = blk(m0)%ke

		is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
		js0 = blk(m0)%js0;je0 = blk(m0)%je0
		ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0
		
		write(debugFile,"('result/debug_invJac_'I10.10'.dat')"),m0
		write(*,*) "Writing ", debugFile
		!!call DEBUG_OUTPUT_2D_1(debugFile,m0,blk(m0)%inv_j,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		call DEBUG_OUTPUT_3D_1(debugFile,m0,blk(m0)%inv_j,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
		write(debugFile,"('result/debug_alpha_'I10.10'.dat')"),m0
		write(*,*) "Writing ", debugFile
		!!call DEBUG_OUTPUT_2D_3(debugFile,m0,blk(m0)%aalpha,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		call DEBUG_OUTPUT_3D_3(debugFile,m0,blk(m0)%aalpha,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
		write(debugFile,"('result/debug_beta_'I10.10'.dat')"),m0
		write(*,*) "Writing ", debugFile
		!!call DEBUG_OUTPUT_2D_3(debugFile,m0,blk(m0)%bbeta,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		call DEBUG_OUTPUT_3D_3(debugFile,m0,blk(m0)%bbeta,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
		write(debugFile,"('result/debug_gamma_'I10.10'.dat')"),m0
		write(*,*) "Writing ", debugFile
		!!call DEBUG_OUTPUT_2D_3(debugFile,m0,blk(m0)%ggamma,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		call DEBUG_OUTPUT_3D_3(debugFile,m0,blk(m0)%ggamma,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	end do
	!!pause
	end if
	
	if(myid .eq. root)then
		write(*,*) '***********************************************************'
		write(*,*) 'Finish computing jacobian matics!                          '
		write(*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*

	!! implement physical boundary condition
	call physical_bc

	do m0 = 1, blk_loop
		is = blk(m0)%is;ie = blk(m0)%ie
		js = blk(m0)%js;je = blk(m0)%je
		ks = blk(m0)%ks;ke = blk(m0)%ke
		
		call get_primitivevariables(m0,blk(m0)%pri_v,blk(m0)%q,gamma,cv,is,ie,js,je,ks,ke)
	end do
	
	!!***************************************************!!
	!!initial flowfield output
	call jacobian_output
	call geom_output
	call flowfield_output
	call surface_output 
	!!*

	!!creat log file
	if(myid .eq. root)then
		open(99,file = 'result/clcd_history.dat',status = 'replace')
		write (99,*) 'variables="flow time","cl","cd"'
		close(99)
		open(99,file = 'result/resi_history.dat',status = 'replace')
		write (99,*) 'variables="timestep","rou","momentum u","momentum v","momentum w","energy","tur1","tur2"'
		close(99)
              
	end if	

	!!***************************************************!!
	!!               end initialization                  !!
	!!***************************************************!!

	!!***************************************************!!
	!!                  Solver Start                     !!
	!!***************************************************!!  
	if(myid .eq. root)then
		write(*,*) '***********************************************************'
		write(*,*) 'Initialization Procedure Finished!                         '
		write(*,*) 'Press ENTER to start the iteration......                   '
		write(*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
	!!*
	
	!!computation start
	!!get start time
	ts = MPI_WTIME()
	te = MPI_WTIME()
	t1 = MPI_WTIME()
	t2 = MPI_WTIME()
	!!*
	
	if(iflag_timeadvance .ne. iflag_1stimplicit .and. iflag_timeadvance .ne. iflag_2ndcrank) then
		write(*,*) "NOT supported!"
		stop
	end if

	!!***************************************************!!
	!! update main equation in buffer block
	call physical_bc
	call UpdateBuffer
	call average1
	call average2
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

	if (isDebug .eq. 1) then
		do m0 = 1,blk_loop
			is = blk(m0)%is;ie = blk(m0)%ie
			js = blk(m0)%js;je = blk(m0)%je
			ks = blk(m0)%ks;ke = blk(m0)%ke

			is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
			js0 = blk(m0)%js0;je0 = blk(m0)%je0
			ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

			write(debugFile,"('result/debug_init_q_'I10.10'.dat')"),m0
			write(*,*) "Writing ", debugFile
			call DEBUG_OUTPUT_2D_5(debugFile,m0,blk(m0)%q,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
		end do
	end if
	!!***************************************************!!

	!!blk(m0)%rhsv = 0.d0
	!!main loop
	outinterval = 1
	do timestep = 1,ntstep
		if(iflag_timeadvance .eq. iflag_1stimplicit) then
			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1

				call get_primitivevariables(m0,blk(m0)%pri_v,blk(m0)%q,gamma,cv,is,ie,js,je,ks,ke)

				call inviscid_flux_1(blk(m0)%rhsi,blk(m0)%q,blk(m0)%dxidx,blk(m0)%inv_j, &
									 is,ie,js,je,ks,ke, &
									 is0,ie0,js0,je0,ks0,ke0, &
									 is1,ie1,js1,je1,ks1,ke1,gamma)		

				if (isDebug .eq. 1) then
					write(debugFile,"('result/debug_rhsi_'I10.10'.dat')"),m0
					write(*,*) "Writing ", debugFile
					!!call DEBUG_OUTPUT_2D_5(debugFile,m0,blk(m0)%rhsi,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
					call DEBUG_OUTPUT_3D_5(debugFile,m0,blk(m0)%rhsi,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
				end if
				

				if(iflag_solver .eq. iflag_nssolver) then
					call viscous_flux(blk(m0)%rhsv,blk(m0)%pri_v,blk(m0)%amu,blk(m0)%vor,blk(m0)%dudx,blk(m0)%tao, &
									  blk(m0)%dxidx,blk(m0)%inv_j,tinf,prl,prt,gamma,cv,re, &
									  is,ie,js,je,ks,ke, &
									  is0,ie0,js0,je0,ks0,ke0, &
									  is1,ie1,js1,je1,ks1,ke1,m0)

					if (isDebug .eq. 1) then
						write(debugFile,"('result/debug_rhsv_'I10.10'.dat')"), m0
						write(*,*) "Writing ", debugFile
						!!call DEBUG_OUTPUT_2D_5(debugFile,m0,blk(m0)%rhsv,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
						call DEBUG_OUTPUT_3D_5(debugFile,m0,blk(m0)%rhsv,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
					end if
				end if
				!!*

				!!implicit time-advancement
				call find_dt(blk(m0)%dt,dtmax,dtmin,cflmax,cflmin,blk(m0)%pri_v,blk(m0)%amu, &
							 blk(m0)%dxidx,re,cfl,is,ie,js,je,ks,ke,iflag_time,iflag_steady,iflag_unsteady)

				if (isDebug .eq. 1) then
					write(debugFile,"('result/debug_dt_'I10.10'.dat')"),m0
					write(*,*) "Writing ", debugFile
					!!call DEBUG_OUTPUT_2D_1(debugFile,m0,blk(m0)%dt,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
					!!call DEBUG_OUTPUT_3D_1(debugFile,m0,blk(m0)%dt,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
				end if

				call implicit1st(blk(m0)%q,blk(m0)%dt,blk(m0)%pri_v,blk(m0)%rhs,blk(m0)%rhsi,blk(m0)%rhsv, &
								 blk(m0)%dxidx,blk(m0)%inv_j, &
							     is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1, &
								 gamma,cfl,blk(m0)%x,blk(m0)%y,blk(m0)%z)

				!!do k = ks1,ke1
				!!do j = js1,je1
				!!do i = is1,ie1
				!!	blk(m0)%q(:,i,j,k) = blk(m0)%q(:,i,j,k) + blk(m0)%dt(i,j,k) * & 
				!!						(blk(m0)%rhsi(:,i,j,k))/blk(m0)%inv_j(i,j,k)				
				!!end do
				!!end do
				!!end do
				!!*	
				
				if(iflag_turbulence .eq. iflag_sa .and. iflag_solver .eq. iflag_nssolver) then
					call sa_model(blk(m0)%nut,blk(m0)%sa_rhs,blk(m0)%dt,blk(m0)%pri_v,blk(m0)%amu,blk(m0)%dudx,blk(m0)%vor, &
								  blk(m0)%dxidx,blk(m0)%inv_j,blk(m0)%aalpha,blk(m0)%bbeta,blk(m0)%ggamma, &
								  blk(m0)%length,blk(m0)%dist,blk(m0)%spacing,re,cfl, &
								  is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1)	
				end if
			end do

			!!update boundary condition
			call physical_bc
			call UpdateBuffer
			call average1
			call average2

			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1
				call get_primitivevariables(m0,blk(m0)%pri_v,blk(m0)%q,gamma,cv,is,ie,js,je,ks,ke)
			end do
		end if

		if (iflag_timeadvance .eq. iflag_2ndcrank) then 
			do sub = 1,3
			call physical_bc
			call updatebuffer
			call average1
			call average2

			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1

				call get_primitivevariables(m0,blk(m0)%pri_v,blk(m0)%q,gamma,cv,is,ie,js,je,ks,ke)

				call inviscid_flux_1(blk(m0)%rhsi,blk(m0)%q,blk(m0)%dxidx,blk(m0)%inv_j, &
									 is,ie,js,je,ks,ke, &
									 is0,ie0,js0,je0,ks0,ke0, &
									 is1,ie1,js1,je1,ks1,ke1,gamma)		

				if(iflag_solver .eq. iflag_nssolver) then
					call viscous_flux(blk(m0)%rhsv,blk(m0)%pri_v,blk(m0)%amu,blk(m0)%vor,blk(m0)%dudx,blk(m0)%tao, &
									  blk(m0)%dxidx,blk(m0)%inv_j,tinf,prl,prt,gamma,cv,re, &
									  is,ie,js,je,ks,ke, &
									  is0,ie0,js0,je0,ks0,ke0, &
									  is1,ie1,js1,je1,ks1,ke1,m0)
				end if
				!!*

				!!implicit time-advancement
				call find_dt(blk(m0)%dt,dtmax,dtmin,cflmax,cflmin,blk(m0)%pri_v,blk(m0)%amu, &
							 blk(m0)%dxidx,re,cfl,is,ie,js,je,ks,ke,iflag_time,iflag_steady,iflag_unsteady)

				call crank2nd(blk(m0)%q,blk(m0)%q0,blk(m0)%dt,blk(m0)%pri_v,blk(m0)%rhs,blk(m0)%rhsi,blk(m0)%rhsv, &
							  blk(m0)%dxidx,blk(m0)%inv_j, &
							  is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1, &
							  gamma,m0,sub)

				!!write(debugFile,"('result/debug_q_'I10.10'.dat')"),sub
				!!write(*,*) "debugging ", sub
				!!call DEBUG_OUTPUT_2D_5(debugFile,m0,blk(m0)%q,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
			end do !!* end block loop
			end do !!* end sub iteration loop

			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1

				call update_q0(blk(m0)%q0,blk(m0)%q,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
			end do

			!!turbulence
			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1
				
				if(iflag_turbulence .eq. iflag_sa .and. iflag_solver .eq. iflag_nssolver) then
					call sa_model(blk(m0)%nut,blk(m0)%sa_rhs,blk(m0)%dt,blk(m0)%pri_v,blk(m0)%amu,blk(m0)%dudx,blk(m0)%vor, &
								  blk(m0)%dxidx,blk(m0)%inv_j,blk(m0)%aalpha,blk(m0)%bbeta,blk(m0)%ggamma, &
								  blk(m0)%length,blk(m0)%dist,blk(m0)%spacing,re,cfl, &
								  is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0,is1,ie1,js1,je1,ks1,ke1)	
				end if		
				
				!!update boundary condition
				call physical_bc
				call updatebuffer
				call average1
				call average2	
			end do
			!!*  end turbulence part

			do m0 = 1,blk_loop
				is  = blk(m0)%is; ie  = blk(m0)%ie
				js  = blk(m0)%js; je  = blk(m0)%je
				ks  = blk(m0)%ks; ke  = blk(m0)%ke
		
				is0 = blk(m0)%is0;ie0 = blk(m0)%ie0
				js0 = blk(m0)%js0;je0 = blk(m0)%je0
				ks0 = blk(m0)%ks0;ke0 = blk(m0)%ke0

				is1 = blk(m0)%is1;ie1 = blk(m0)%ie1
				js1 = blk(m0)%js1;je1 = blk(m0)%je1
				ks1 = blk(m0)%ks1;ke1 = blk(m0)%ke1
				call get_primitivevariables(m0,blk(m0)%pri_v,blk(m0)%q,gamma,cv,is,ie,js,je,ks,ke)
			end do

		end if

		if(mod(timestep,outinterval) .eq. 0)then  !end in line639
			t1 = t2
			t2 = MPI_WTIME()
			
			t2mt1 = real(t2-t1)
			call MPI_Reduce(t2mt1,tmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
			call MPI_Reduce(t2mt1,tmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,root,MPI_COMM_WORLD,ierr)
			if (myid .eq. root)then
				write(*,*) '***********************************************************'
				write(*,*) "TimeStep =",timestep
				write(*,10) tmax,tmin,outinterval
10			format("Max CPU Time =",f6.3," Min CPU Time =",f6.3, " sec/(" ,I2, " timesteps)")
			end if
			!!*
			call MPI_Reduce(dtmax ,ddtmax ,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
			call MPI_Reduce(dtmin ,ddtmin, 1,MPI_DOUBLE_PRECISION,MPI_MIN,root,MPI_COMM_WORLD,ierr)
			call MPI_Reduce(cflmax,dcflmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
			call MPI_Reduce(cflmin,dcflmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,root,MPI_COMM_WORLD,ierr)
			
			if (myid .eq. root)then
				if     (iflag_time .eq. iflag_steady  )then
					write(*,15) ddtmax,ddtmin
15				    format("Max Timestep =",e10.4," Min Timestep =",e10.4)
				else if(iflag_time .eq. iflag_unsteady)then
					write(*,20) dcflmax,dcflmin
20				    format("Max CFL =",e10.4," Min CFL =",e10.4)				
				end if
			end if    		
			call residual
			call clcd_output
		end if  !end of 10 steps output
		!!*
                
		if(mod(timestep,WriteInterval) .eq. 0)then
			call flowfield_output
			call surface_output
		end if
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	end do
	!!*end timeloop
	
	!!get end time
	if (myid .eq. root) then 
		te = MPI_WTIME()
	end if
	!!*
	
	if(myid .eq. root)then
		print *, 'computation finished!!'
		print *, 'Total CPU time =', real((te-ts)/10000)
	end if
	
	!!mpi end
	call MPI_FINALIZE(ierr)
	!!*
	
	stop
	
end program    


subroutine GetBufferLength(numBuffer, inviscidScheme)
	use flag_var
	implicit none

	integer:: numBuffer, inviscidScheme

	numBuffer = 0;
	if (inviscidScheme .eq. iflag_3rdweno) then
		numBuffer = 2
	end if
		
	if (numBuffer .lt. 1) then
		print *, "Buffer length is not correct, please check the iflag_inviscid!"
		stop
	end if

	return
end subroutine

subroutine CheckBoundary()
	use index_var
	use blk_var
	use glbindex_var
	implicit none
	integer :: iface

	integer :: isPhyBound(6), isBlkInterface(6)

	do m = 1,nblock
	isPhyBound = 0
	isBlkInterface = 0
	do ksub = 1,blk0(m)%num_subface
		if (blk0(m)%bc0(ksub)%blk_t .lt. 0) then
			isPhyBound(blk0(m)%bc0(ksub)%face) = 1
		end if

		if (blk0(m)%bc0(ksub)%blk_t .gt. 0) then
			isBlkInterface(blk0(m)%bc0(ksub)%face) = 1
		end if
	end do

	do iface = 1,6
		if (isPhyBound(iface) * isBlkInterface(iface) .eq. 1) then 
			print *, "Please check boundary"
			print *, "block:", m
			do ksub = 1,blk0(m)%num_subface
				print *, "ksub:", ksub, blk0(m)%bc0(ksub)%blk_t
			end do
			stop
		end if
	end do
	end do

	return 
end subroutine

subroutine update_q0(q0,q,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	implicit none
	
	integer :: i,j,k
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: q0(-1:0,5,is:ie,js:je,ks:ke)
	real*8  :: q (     5,is:ie,js:je,ks:ke)

	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		q0(-1,:,i,j,k) = q0(0,:,i,j,k)
		q0( 0,:,i,j,k) = q (  :,i,j,k)
	end do
	end do
	end do
	
	return
end subroutine