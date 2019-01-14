!!==========================================================
!!	read control file                                       
!!	                                                        
!!	purpose: read global variables in control file          
!!					                                                
!!	root rank is used to read the control file and then     
!!	broadcast to other ranks                                
!!                                                          
!!	author: jiamin xu                                       
!!	date:   2013.09.06                                      
!!==========================================================
subroutine init_readcontrolfile
	use flag_var
	use mpi_var
	use ns_const
	use output_var
	implicit none
	
	character (len=180) :: inp
	character (len=500) :: buff
	
	if (myid .eq. root) then
		!!reading control.in
		write(*,*) '***********************************************************'
		write(*,*) 'reading control file: ./input/Control.in                   '
		write(*,*) '***********************************************************'
		inp = "input/Control.in"
		open(99,file = inp)
		
		read(99,*) buff
		read(99,*) buff
		read(99,*) iflag_dimension,iflag_solver,iflag_time,cfl,dt                                               !!6
		read(99,*) buff
		read(99,*) buff
		read(99,*) ma,re,tinf,aoa                                                                               !!5
		read(99,*) buff
		read(99,*) buff
		read(99,*) scale_x,scale_y,scale_z,lref,sref,fltr                                                       !!7
		read(99,*) buff
		read(99,*) buff
		read(99,*) iflag_splittingtype,iflag_inviscid,iflag_filter,iflag_timeadvance,iflag_turbulence,iflag_des !!6
		read(99,*) buff
		read(99,*) buff
		read(99,*) writeinterval,ntstep,grid,bcinfo,iflag_init,restartfile,nfile,restart_timestep               !!8
        read(99,*) buff
        read(99,*) buff
        read(99,*) nobs,ob_dist,acoustic_length,iflag_acoustic,Start_Acsoutic,End_Acsoutic,iflag_avg,Start_Avg,End_Avg,Avg_interval                    !!6
		close(99)
		
		write (*,*)"Computation Condition:"
		write (*,*)"reynolds number:", re
		write (*,*)"mach number:    ", ma
		write (*,*)"angle of attack:", aoa
		if(iflag_time .eq. iflag_steady) then
			write (*,*)"CFL number:     ", cfl
		else if(iflag_time .eq. iflag_unsteady) then
			write (*,*)"time step:      ", dt
		endif		
		write(*,*) '***********************************************************'
		!!*end output calc condition
		write(*,*) '***********************************************************'
		write(*,*)  "Finish reading control file!"
		write(*,*) '***********************************************************'
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	!!*finish reading control.in
	
	!!broadcast the control parameters to other threads
	!!line 1
	call MPI_BCAST(iflag_dimension,     1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_solver,        1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_time,          1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cfl,                 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dt,                  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	!!line 2
	call MPI_BCAST(ma,                  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(re,                  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tinf,                1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aoa,                 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	!!line 3
	call MPI_BCAST(scale_x,             1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(scale_y,             1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(scale_z,             1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lref,                1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sref,                1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fltr,                1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	!!line 4
	call MPI_BCAST(iflag_splittingtype, 1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_inviscid,      1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_filter,        1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_timeadvance,   1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_turbulence,    1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iflag_des,           1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	!!line5
	call MPI_BCAST(iflag_init,          1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(writeinterval,       1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ntstep,              1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(restart_timestep,    1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
	!!line6
    call MPI_BCAST(nobs,                1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ob_dist,             1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(acoustic_length,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iflag_acoustic      ,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Start_Acsoutic,      1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(End_Acsoutic,        1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iflag_Avg,           1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Start_Avg,           1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(End_Avg,             1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Avg_interval,        1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	 	
	return
end subroutine
