module mpi_var
	implicit none
	include "mpif.h"
	!!system parameter
	integer :: ierr
	integer,save :: myid
	integer,save :: numprocs
	integer :: namelen
	character(MPI_MAX_PROCESSOR_NAME):: processor_name
	integer :: status(MPI_STATUS_SIZE)
	!!*
	
	integer,parameter :: root = 0   !!root: rank = 0
	integer           :: dest_id    !!rank which receive the data
	integer           :: source_id  !!rank which send the data 
end module mpi_var