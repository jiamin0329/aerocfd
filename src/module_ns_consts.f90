module ns_const
	implicit none
	!!gas constants (read from input control file)
	real*8, save :: ma                      !!Mach Number
	real*8, save :: re                      !!Reynolds Number
	real*8, save :: aoa                     !!Angle of Attack
	real*8, save :: tinf                    !!reference temperature
	real*8, save :: d_0,u_0,t_0,p_0,c_0     !!density, velocity, temperature, pressure, soundspeed(farfield)
	
	!!geom constants (read from input control file)
	real*8, save :: lref,sref               !!reference length, reference area
	real*8, save :: scale_x,scale_y,scale_z !!scaling factor of three dimensions
	
	real*8, save :: cfl,dt                  !!cfl number, time step scaling
	real*8, save :: dtmax, dtmin            !!maximum timestep, minimum timestep
	real*8, save :: ddtmax,ddtmin           !!maximum timestep, minimum timestep to output
	real*8, save :: cflmax,cflmin           !!maximum timestep, minimum timestep
	real*8, save :: dcflmax,dcflmin         !!maximum timestep, minimum timestep to output
	integer,save :: subiteration            !!subiteration number(needed in unsteady dual time-stepping method)
	integer,save :: timestep                !!current timestep number
	real*8, save :: flow_time               !!current real flow time
	integer,save :: ntstep                  !!total timestep
	real*8, save :: cp,cv                   !!specific heat at constant pressure��specific heat at constant volume(derived in sub_init_flowfield)
	real*8, save :: fltr                    !!numerical filter coefficient
	!!*
	!!acoustic
	integer,save :: nobs
	real*8, save :: ob_dist
	real*8, save :: acoustic_length

	!!
	real*8, parameter :: gamma = 1.4d0      !!gas constant
	real*8, parameter :: prt   = 0.90d0     !!prantl constant(laminar)
	real*8, parameter :: prl   = 0.72d0     !!prantl constant(turbulent)
	real*8, parameter :: pi    = 3.1415926535897932d0
end module ns_const
