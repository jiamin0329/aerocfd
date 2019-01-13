module index_var
	implicit none
	
	integer      :: m                 !!block index(total)
	integer      :: m0                !!block index(current processor)
	integer      :: myn,lastn         !!totalblock number of current processor
	integer,save :: blk_loop          !!totalblock number of current processor

	integer      :: ksub              !!subface index
end module index_var

