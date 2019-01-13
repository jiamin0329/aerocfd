module output_var
	implicit none
	!!time prediction
	real*8 :: t1    = 0.d0
	real*8 :: t2    = 0.d0
	real*8 :: ts    = 0.d0
	real*8 :: te    = 0.d0
	real*8 :: t2mt1 = 0.d0
	real*8 :: tmax  = 0.d0
	real*8 :: tmin  = 0.d0
	!!*
	
	integer,save :: writeinterval       !!outputinterval
	integer,save :: nfile               !!file number to restart from

	character(len = 180) :: grid        !!grid file name
	character(len = 180) :: bcinfo      !!bcinfo file name
	character(len = 180) :: restartfile !!restart file name
	character(len = 180) :: acoustic
	character(len = 180) :: pts
end module output_var
