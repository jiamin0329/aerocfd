module flag_var
	implicit none
	!!******************************************************************************************************
	integer,save :: iflag_dimension
	integer,save :: iflag_solver
	integer,save :: iflag_time
	integer,save :: iflag_inviscid
	integer,save :: iflag_splittingtype
	integer,save :: iflag_filter
	integer,save :: iflag_timeadvance
	integer,save :: iflag_turbulence
	integer,save :: iflag_des
	integer,save :: iflag_init
	integer,save :: iflag_blockinterface
	integer,save :: iflag_walldisttype
	!!
	integer,save :: iflag_acoustic
	integer,save :: restart_timestep
	integer,save :: Start_Acsoutic
    integer,save :: End_Acsoutic
    integer,save :: iflag_avg
    integer,save :: Start_Avg
    integer,save :: End_Avg
    integer,save :: Avg_interval
	!!
	
	integer :: bufferLength 

	integer :: isDebug = 0

	!!******************************************************************************************************
	integer,parameter :: iflag_2d             =  2 !!(iflag_dimension)dimension : 2d
	integer,parameter :: iflag_3d             =  3 !!(iflag_dimension)dimension : 3d
	!!******************************************************************************************************
	integer,parameter :: iflag_nssolver       =  0 !!(iflag_solver)solver type: ns
	integer,parameter :: iflag_eulersolver    =  1 !!(iflag_solver)solver type: euler
	!!******************************************************************************************************
	integer,parameter :: iflag_steady         =  0 !!(iflag_time)steady simulation
	integer,parameter :: iflag_unsteady       =  1 !!(iflag_time)unsteady simulation
	!!******************************************************************************************************
	integer,parameter :: iflag_1steuler       =  1 !!(iflag_timeadvance)time advance method: 1st euler
	integer,parameter :: iflag_1stimplicit    =  2 !!(iflag_timeadvance)time advance method: 1st implicit method
	integer,parameter :: iflag_2ndcrank       =  3 !!(iflag_timeadvance)time advance method: 2nd crank-nicholson method             
	integer,parameter :: iflag_rk3            =  4 !!(iflag_timeadvance)time advance method: 3rd tvd runge-kutta method
	integer,parameter :: iflag_lusgs          =  5 !!(iflag_timeadvance)time advance method: lusgs method
	integer,parameter :: iflag_lddrk          =  6 !!(iflag_timeadvance)time advance method: 2N-storage ldd rk
	!!******************************************************************************************************
	integer,parameter :: iflag_4thdrp         = -1 !!(iflag_inviscid)central difference method: 4th drp scheme
	integer,parameter :: iflag_4thcompact     = -2 !!(iflag_inviscid)central difference method: 4th compact scheme(lele)
	integer,parameter :: iflag_6thcompact     = -3 !!(iflag_inviscid)central difference method: 6th compact scheme(lele)
	integer,parameter :: iflag_8thcompact     = -4 !!(iflag_inviscid)central difference method: 8th compact scheme(lele)
	integer,parameter :: iflag_10thcompact    = -5 !!(iflag_inviscid)central difference method:10th compact scheme(lele)
	!!******************************************************************************************************
	integer,parameter :: iflag_asd            =  1 !!(iflag_filter)numerical filter: ASD
	integer,parameter :: iflag_6thfilter      =  2 !!(iflag_filter)numerical filter: 4th filter
	integer,parameter :: iflag_8thfilter      =  3 !!(iflag_filter)numerical filter: 6th filter
	integer,parameter :: iflag_10thfilter     =  4 !!(iflag_filter)numerical filter: 8th filter
	!!******************************************************************************************************
	integer,parameter :: iflag_stegerwarming  =  1 !!(iflag_splittingtype)flux splitting method: steger-warming
	integer,parameter :: iflag_vanleer        =  2 !!(iflag_splittingtype)flux splitting method: vanleer
	integer,parameter :: iflag_roe            =  3 !!(iflag_splittingtype)flux splitting method: roe
	integer,parameter :: iflag_ausm           =  4 !!(iflag_splittingtype)flux splitting method: ausm
	integer,parameter :: iflag_hll            =  5 !!(iflag_splittingtype)flux splitting method: hll
	integer,parameter :: iflag_hllc           =  6 !!(iflag_splittingtype)flux splitting method: hllc
	!!******************************************************************************************************
	integer,parameter :: iflag_1stupwind      =  1 !!(iflag_inviscid)1st upwind scheme
	integer,parameter :: iflag_3rdupwind      =  2 !!(iflag_inviscid)3rd upwind scheme
	integer,parameter :: iflag_3rdmuscl       =  3 !!(iflag_inviscid)3rd muscl(van albada limiter)
	integer,parameter :: iflag_3rdweno        =  4 !!(iflag_inviscid)3rd weno scheme
	integer,parameter :: iflag_5thweno        =  5 !!(iflag_inviscid)5th weno scheme
	integer,parameter :: iflag_7thweno        =  6 !!(iflag_inviscid)7th weno scheme
	integer,parameter :: iflag_9thweno        =  7 !!(iflag_inviscid)9th weno scheme
	integer,parameter :: iflag_5thwcns        =  8 !!(iflag_inviscid)5th wcns scheme
	integer,parameter :: iflag_7thwcns        =  9 !!(iflag_inviscid)7th wcns scheme
	!!******************************************************************************************************
	integer,parameter :: iflag_laminar        =  0 !!(iflag_turbulence)laminar flow
	integer,parameter :: iflag_bl             =  1 !!(iflag_turbulence)baldwin-lomax model
	integer,parameter :: iflag_sa             =  2 !!(iflag_turbulence)spalart-allmaras model
	integer,parameter :: iflag_kwsst_menter   =  3 !!(iflag_turbulence)k-w sst model
	integer,parameter :: iflag_les            =  4 !!(iflag_turbulence)large-eddy simulation
	!!******************************************************************************************************
	integer,parameter :: iflag_des97          =  1 !!(iflag_des)detached-eddy simulation 1997
	integer,parameter :: iflag_ddes           =  2 !!(iflag_des)delayed detached-eddy simulation
	!!******************************************************************************************************
	integer,parameter :: iflag_intp           =  1 !!(iflag_blockinterface)average procedure
	integer,parameter :: iflag_gcic           =  2 !!(iflag_blockinterface)generalize characteristic interface conditions
	integer,parameter :: iflag_overlap        =  3 !!(iflag_blockinterface)overlap method
	!!******************************************************************************************************	
	!!boundary condition
	integer,parameter :: bc_user     = -90
	integer,parameter :: bc_periodic = -60
	integer,parameter :: bc_inflow   = -50
	integer,parameter :: bc_outflow  = -40
	integer,parameter :: bc_symm     = -30
	integer,parameter :: bc_wall     = -20
	integer,parameter :: bc_farfield = -10
	!!******************************************************************************************************
end module
