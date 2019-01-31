module blk_var
	implicit none

	integer,save :: nblock         !!total number of blocks
	integer,save :: totalnodes     !!total nodes
	!!*
	
	!!subface info
	type blockinterface
		integer :: f_no,face,is,ie,js,je,ks,ke,dim_f,blk_t,subface_t,ori,fwh
		integer :: rank_c,rank_t
		!!buffer block 
		real*8,pointer,dimension(:,:,:,:,:) :: buffer_send
		real*8,pointer,dimension(:,:,:,:,:) :: buffer_recv
		real*8,pointer,dimension(:,:,:,:)   :: buffer
	end type blockinterface
	!!*

	!!block info
	type block_type
		integer :: blk_no
		integer :: ni,nj,nk,num_subface
		integer :: is,  ie,  js,  je,  ks,  ke
		integer :: is0, ie0, js0, je0, ks0, ke0
		integer :: is1, ie1, js1, je1, ks1, ke1
		!!coordinates
		real*8,pointer,dimension(:,:,:):: x,y,z
		!!conservative variables
		!!1: d
		!!2: d*u
		!!3: d*v
		!!4: d*w
		!!5: d*E
		real*8,pointer,dimension(:,:,:,:):: q
		!!primitive variables
		!!1:d  2:u  3:v
		!!4:w  5:p  6:t  7:c
		real*8,pointer,dimension(:,:,:,:):: pri_v 
		!!inviscid right hand side term
		real*8,pointer,dimension(:,:,:,:):: rhsi
		!!viscous right hand side terms
		real*8,pointer,dimension(:,:,:,:):: rhsv
		!!right hand side term
		real*8,pointer,dimension(:,:,:,:):: rhs
		!!for viscous flux
		!!gradient variables
		!!1,i,j,k,1 ==> dudx	  2,i,j,k,1 ==> dvdx   3,i,j,k,1 ==> dwdx
		!!1,i,j,k,2 ==> dudy	  2,i,j,k,2 ==> dvdy   3,i,j,k,2 ==> dwdy
		!!1,i,j,k,3 ==> dudz	  2,i,j,k,3 ==> dvdz   3,i,j,k,3 ==> dwdz
		real*8,pointer,dimension(:,:,:,:):: dudx
		!!stress tensor
		!!1: taoxx	4: taoyy
		!!2: taoxy  5: taoyz
		!!3: taoxz  6: taozz
		real*8,pointer,dimension(:,:,:,:):: tao
		!!vorticity and viscosity coefficients
		real*8,pointer,dimension(:,:,:):: vor
		!!viscosity coefficient
		!!1: amul, 2: amut, 3: amu
		real*8,pointer,dimension(:,:,:,:):: amu
		!!time advancement
		!!non-physical time step
		real*8,pointer,dimension(:,:,:):: dt
		
		!!jacobian related parameters
		!!1,:,:,: ==> dxi/dx  2,:,:,: ==> deta/dx	 3,:,:,: ==> dzeta/dx
		!!4,:,:,: ==> dxi/dy  5,:,:,: ==> deta/dy	 6,:,:,: ==> dzeta/dy
		!!7,:,:,: ==> dxi/dz  8,:,:,: ==> deta/dz	 9,:,:,: ==> dzeta/dz
		real*8,pointer,dimension(:,:,:,:):: dxidx
		!!1/J
		real*8,pointer,dimension(:,:,:):: inv_j  
	
		!!laplacian needed
		real*8,pointer,dimension(:,:,:,:):: aalpha,bbeta,ggamma

		!!geom
		real*8,pointer,dimension(:,:,:):: dist,spci,spcj,spck
		real*8,pointer,dimension(:,:,:):: spacing
		real*8,pointer,dimension(:,:,:):: vol
		real*8,pointer,dimension(:,:,:):: ratio

		!!turbulence model
		!!turbulent length
		real*8,pointer,dimension(:,:,:):: length
		!!spalart-allmaras model
		real*8,pointer,dimension(:,:,:):: nut
		real*8,pointer,dimension(:,:,:):: sa_rhs

		type(blockinterface), pointer,dimension(:):: bc
	end type block_type
	!!*end block type

	type(block_type),save,dimension(:),allocatable,target:: blk
end module blk_var
