module sa_var
	implicit none
	
	real*8,parameter :: dta = 2.d0/3.d0, kar = 0.41d0
	real*8,parameter :: cb1 = 0.1355d0, cb2 = 0.622d0
	real*8,parameter :: cw1 = cb1/kar/kar + (1.d0 + cb2)/dta, cw2 = 0.3d0, cw3 = 2.d0
	real*8,parameter :: ct3 = 1.2d0, ct4 = 0.5d0
	real*8,parameter :: cv1 = 7.1d0
	real*8 :: fv1,fv2,ft2,fw,g,r,s,leix
	
	real*8,save  :: nutinf = 1.341946d0
end module sa_var
