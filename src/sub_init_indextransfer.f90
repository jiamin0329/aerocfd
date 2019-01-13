!!m to rank
subroutine m2node(dest_id,m,myn,lastn,nblock,numprocs)
	implicit none
	
	integer :: dest_id
	integer :: m
	integer :: myn,lastn
	integer :: nblock
	integer :: numprocs
	
	!!id of thread to receive
	if (m .le. myn)then
		dest_id = 0
	else if (m .gt. nblock-lastn)then
		dest_id = numprocs-1
	else
		dest_id = (m+myn-1)/myn-1
	end if
	
	return
end subroutine
    
!!m0 to m
subroutine m02m(m,m0,myid,myn)
	implicit none
	
	integer :: m
	integer :: m0,myid,myn
	
	m = myid*myn + m0
	
	return
end subroutine