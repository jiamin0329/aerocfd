module glbindex_var
	implicit none
	
	type blockinterface0
		integer :: f_no,face,is,ie,js,je,ks,ke,dim_f,ori,fwh,blk_t,subface_t
		integer :: rank_c,rank_t
		real*8,allocatable,dimension(:,:) :: x,y,z
	end type blockinterface0
	
	type block0
		integer :: blk_no,ni,nj,nk,num_subface
		integer ::  is,  ie,  js,  je,  ks,  ke
		integer :: is0, ie0, js0, je0, ks0, ke0
		integer :: is1, ie1, js1, je1, ks1, ke1
        real*8,pointer,dimension(:,:,:) :: dist
		type(blockinterface0), pointer,dimension(:):: bc0
	end type block0

	type(block0),save,dimension(:),allocatable,target:: blk0
	
end module glbindex_var