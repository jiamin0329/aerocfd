!!=========================================================!!
!! geomtry output                                          !!
!!                                                         !!
!! coordination: x,y,z                                     !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!! geom variables: spci,spcj,spck                          !!
!!                 spacing,dist,vol                        !!
!!primitive variables
!!1:d  2:u  3:v
!!4:w  5:p  6:t  7:c
!!                                                         !!
!! Author: jiamin xu                                       !!
!! Date:   2012.12.05                                      !!
!!=========================================================!!
subroutine flowfield_output
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: i,j,k
	integer:: idim,jdim,kdim
	character(len = 180) :: flow
	integer :: is, ie, js, je, ks, ke !!block dimension
	!!*********************************************************************!!
	!!geom file output
	write(flow,"('result/flowfield_'I10.10'_'I4.4'.dat')"),timestep,myid
	if(myid .eq. root)then
		write(*,*) "writing flowfield......"
	end if
	open (99,file = flow)
	
	if      (iflag_dimension .eq. iflag_2d) then
		k = 1
		do m0 = 1,blk_loop
			is = blk(m0)%is0; ie = blk(m0)%ie0
			js = blk(m0)%js0; je = blk(m0)%je0
			idim = ie-is+1
			jdim = je-js+1

			write(99,*) "variables = x,y,d,u,v,w,p,nut"
			write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
			do j = js,je
			do i = is,ie
				write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),   &
							blk(m0)%pri_v(1,i,j,k),              &
							blk(m0)%pri_v(2,i,j,k),              &
							blk(m0)%pri_v(3,i,j,k),              &
							blk(m0)%pri_v(4,i,j,k),              &
							blk(m0)%pri_v(5,i,j,k),              &
							blk(m0)%nut  (  i,j,k)
			end do
			end do
		end do	
	else if (iflag_dimension .eq. iflag_3d) then
		do m0 = 1,blk_loop
			is = blk(m0)%is0;ie = blk(m0)%ie0
			js = blk(m0)%js0;je = blk(m0)%je0
			ks = blk(m0)%ks0;ke = blk(m0)%ke0

			idim = ie-is+1
			jdim = je-js+1
			kdim = ke-ks+1
			
			write(99,*) "variables = x,y,z,d,u,v,w,p"
			write(99,*) "zone i= ",idim, "j= ",jdim, "k= ",kdim,"datapacking=point"
			
			do k = ks,ke
			do j = js,je
			do i = is,ie
				write(99,*) blk(m0)%x(i,j,k),blk(m0)%y(i,j,k),blk(m0)%z(i,j,k), &
							blk(m0)%pri_v(1,i,j,k),                             &
							blk(m0)%pri_v(2,i,j,k),                             &
							blk(m0)%pri_v(3,i,j,k),                             &
							blk(m0)%pri_v(4,i,j,k),                             &
							blk(m0)%pri_v(5,i,j,k)
			end do
			end do
			end do
		end do
	end if
	close(99)
	return
end subroutine
	!!*********************************************************************!!
!Export FW-H surface(solid wall in the 1.0version,extended to permiable surface later)
subroutine acoustic_output
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none

	integer :: i,j,k,fake_k
	integer :: is, ie, js, je, ks, ke !!block dimension
        integer :: idim,jdim,kdim

if      (iflag_dimension .eq. iflag_2d) then
k = 1  !2D,two dimension
  
		do m0 = 1,blk_loop
                    !write(98,*) "variables = x,y,z,p"
			do ksub = 1,blk(m0)%num_subface
			   if(blk(m0)%bc(ksub)%fwh .eq. 1)then
			           write(acoustic,"('result/acoustic_'I4.4'_'I4.4'.dat')") myn*myid+m0,blk(m0)%bc(ksub)%face
                       write(*,*) acoustic
                       open(98,file=acoustic,status='old',position='append',form='unformatted')
					!face = blk(m0)%bc(ksub)%face
					!!wall dimention
					is = blk(m0)%bc(ksub)%is
					ie = blk(m0)%bc(ksub)%ie
					js = blk(m0)%bc(ksub)%js
					je = blk(m0)%bc(ksub)%je
					!!*
							
  			idim = ie-is+1
  			jdim = je-js+1
  			kdim = 51
  			
  			
				write(98) (timestep-Start_Acsoutic)/10
                                
				if     (blk(m0)%bc(ksub)%face .eq. 1)then
					i=is
						do j = js,je-1
                           do fake_k=1,kdim-1
								write(98) blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
                           end do
                        end do
				else if(blk(m0)%bc(ksub)%face .eq. 2)then
					i=ie
						do j = js,je-1
                           do fake_k=1,kdim-1
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						    end do
						end do
				else if(blk(m0)%bc(ksub)%face .eq. 3)then
					j=js
						do i = is,ie-1	
                            do fake_k=1,kdim-1 !why minus 1�?
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						    end do
                        end do
				else if(blk(m0)%bc(ksub)%face .eq. 4)then
					j=je
						do i = is,ie-1
                           do fake_k=1,kdim-1
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						   end do
                        end do
                 end if
                end if
               end do
              end do
              close(98)
!*********************    3d   ***
else if      (iflag_dimension .eq. iflag_3d) then
!k = 1  !2D,two dimension
  
		do m0 = 1,blk_loop !loop1
                    !write(98,*) "variables = x,y,z,p"
			do ksub = 1,blk(m0)%num_subface !loop2
			   if(blk(m0)%bc(ksub)%fwh .eq. 1)then !if1
			           write(acoustic,"('result/acoustic_'I4.4'_'I4.4'.dat')") myn*myid+m0,blk(m0)%bc(ksub)%face
			if((myn*myid+m0) .eq. 69) then
			write(pts,"('result/pts.dat')")
			open(99,file=pts,status='old',position='append',form='formatted')
			write(99,*) timestep,blk(m0)%pri_v(1:5,17,12,16)
			close(99)
			end if
                       write(*,*) acoustic
                       open(98,file=acoustic,status='old',position='append',form='unformatted')
					!face = blk(m0)%bc(ksub)%face
					!!wall dimention
					is = blk(m0)%bc(ksub)%is
					ie = blk(m0)%bc(ksub)%ie
					js = blk(m0)%bc(ksub)%js
					je = blk(m0)%bc(ksub)%je
					ks = blk(m0)%bc(ksub)%ks
					ke = blk(m0)%bc(ksub)%ke
					!!*
							
  			idim = ie-is+1
  			jdim = je-js+1
  			kdim = ke-ks+1
  			
  			
				write(98) (timestep-Start_Acsoutic)/10
                                
				if     (blk(m0)%bc(ksub)%face .eq. 1)then
					i=is
						do j = js,je-1
                           				 do k=ks,ke-1
								write(98) blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
                          				    end do
                     				   end do
				else if(blk(m0)%bc(ksub)%face .eq. 2)then
					i=ie
						do j = js,je-1
						    do k=ks,ke-1
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						    end do
						end do
				else if(blk(m0)%bc(ksub)%face .eq. 3)then
					j=js
						do i = is,ie-1	
						   do k=ks,ke-1
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						   end do
                      				  end do
				else if(blk(m0)%bc(ksub)%face .eq. 4)then
					j=je
						do i = is,ie-1	
						    do k=ks,ke-1
								write(98)  blk(m0)%pri_v(1,i,j,k),blk(m0)%pri_v(2,i,j,k),&
									blk(m0)%pri_v(3,i,j,k), blk(m0)%pri_v(4,i,j,k),blk(m0)%pri_v(5,i,j,k)
						    end do
                       				 end do
                			end if
               end if !end if 1
           end do !loop2
        end do !loop1
                       close(98)
 end if

!end of export FW-H
!**********************************************
	return
end subroutine

subroutine DEBUG_OUTPUT_2D_1(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: mm
	integer :: i,j,k
	integer:: idim,jdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0
	real*8  :: var(is:ie,js:je,ks:ke)

	!!*********************************************************************!!
	!!geom file output
	open (99,file = filename)
	k = 1
	idim = ie0-is0+1
	jdim = je0-js0+1

	write(99,*) "variables = x,y,var1"
	write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),var(i,j,k)		
	end do
	end do

	!! buffer block i-
	idim = is0-is +1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),var(i,j,k)	
		end do
		end do
	end if

	!! buffer block i+
	idim = ie -ie0+1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),var(i,j,k)	
		end do
		end do
	end if

	!! buffer block j-
	idim = ie0-is0+1
	jdim = js0-js +1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),var(i,j,k)	
		end do
		end do
	end if

	!! buffer block j+
	idim = ie0-is0+1
	jdim = je -je0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = je0,je
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),var(i,j,k)	
		end do
		end do
	end if

	close(99)
	return
end subroutine

subroutine DEBUG_OUTPUT_2D_3(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: i,j,k,mm
	integer:: idim,jdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0
	real*8  :: var(3,is:ie,js:je,ks:ke)

	k = 1
	open (99,file = filename)
	idim = ie0-is0+1
	jdim = je0-js0+1

	write(99,*) "variables = x,y,var1,var2,var3"
	write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
	end do
	end do

	!! buffer block i-
	idim = is0-is +1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)	
		end do
		end do
	end if

	!! buffer block i+
	idim = ie -ie0+1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
	end if

	!! buffer block j-
	idim = ie0-is0+1
	jdim = js0-js +1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
	end if

	!! buffer block j+
	idim = ie0-is0+1
	jdim = je -je0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = je0,je
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &	
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
	end if

	close(99)
	return
end subroutine

subroutine DEBUG_OUTPUT_2D_5(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: mm
	integer :: i,j,k
	integer:: idim,jdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0
	real*8  :: var(5,is:ie,js:je,ks:ke)

	!!*********************************************************************!!
	!!geom file output
	open (99,file = filename)
	k = 1
	idim = ie0-is0+1
	jdim = je0-js0+1

	write(99,*) "variables = x,y,var1,var2,var3,var4,var5"
	write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)		
	end do
	end do

	!! buffer block i-
	idim = is0-is +1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			            var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
	end if

	!! buffer block i+
	idim = ie -ie0+1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
	end if

	!! buffer block j-
	idim = ie0-is0+1
	jdim = js0-js +1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)
		end do
		end do
	end if

	!! buffer block j+
	idim = ie0-is0+1
	jdim = je -je0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "datapacking=point"
			
		do j = je0,je
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
	end if	

	close(99)
	return
end subroutine

subroutine DEBUG_OUTPUT_3D_1(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: i,j,k,mm
	integer:: idim,jdim,kdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: var(is:ie,js:je,ks:ke)

	
	open (99,file = filename)
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ke0-ks0+1

	write(99,*) "variables = x,y,k,var1"
	write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
		
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
	end do
	end do
	end do

	!! buffer block i-
	idim = is0-is +1
	jdim = je0-js0+1
	kdim = ke0-ks0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block i+
	idim = ie -ie0+1
	jdim = je0-js0+1
	kdim = ke0-ks0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block j-
	idim = ie0-is0+1
	jdim = js0-js +1
	kdim = ke0-ks0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block j+
	idim = ie0-is0+1
	jdim = je -je0+1
	kdim = ke0-ks0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = je0,je
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block k-
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ks0-ks +1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks ,ks0
		do j = js0,je0
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block k+
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ke -ke0+1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ke0,ke
		do j = js0,je0
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k),var(i,j,k)
		end do
		end do
		end do
	end if

	close(99)
	return
end subroutine

subroutine DEBUG_OUTPUT_3D_3(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: i,j,k,mm
	integer:: idim,jdim,kdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: var(3,is:ie,js:je,ks:ke)

	
	open (99,file = filename)
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ke0-ks0+1

	write(99,*) "variables = x,y,k,var1,var2,var3"
	write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
		
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
	end do
	end do
	end do

	!! buffer block i-
	idim = is0-is +1
	jdim = je0-js0+1
	kdim = ke0-ks0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)	
		end do
		end do
		end do
	end if

	!! buffer block i+
	idim = ie -ie0+1
	jdim = je0-js0+1
	kdim = ke0-ks0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block j-
	idim = ie0-is0+1
	jdim = js0-js +1
	kdim = ke0-ks0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block j+
	idim = ie0-is0+1
	jdim = je -je0+1
	kdim = ke0-ks0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks0,ke0
		do j = je0,je
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block k-
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ks0-ks +1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ks ,ks0
		do j = js0,je0
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block k+
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ke -ke0+1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,k,var1,var2,var3"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim, "datapacking=point"
			
		do k = ke0,ke
		do j = js0,je0
		do i = is0,ie0
			write(99,*)blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
				        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k)
		end do
		end do
		end do
	end if

	close(99)
	return
end subroutine

subroutine DEBUG_OUTPUT_3D_5(filename,mm,var,is,ie,js,je,ks,ke,is0,ie0,js0,je0,ks0,ke0)
	use blk_var
	use mpi_var
	use index_var
	use output_var
	use ns_const
	use flag_var
	implicit none
	
	integer :: mm
	integer :: i,j,k
	integer:: idim,jdim,kdim
	character(len = 180) :: filename
	integer :: is,ie,js,je,ks,ke
	integer :: is0,ie0,js0,je0,ks0,ke0
	real*8  :: var(5,is:ie,js:je,ks:ke)

	!!*********************************************************************!!
	!!geom file output
	open (99,file = filename)
	kdim = ke0-ks0+1
	idim = ie0-is0+1
	jdim = je0-js0+1

	write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
	write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
		
	do k = ks0,ke0
	do j = js0,je0
	do i = is0,ie0
		write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), blk(mm)%z(i,j,k), &
			        var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)		
	end do
	end do
	end do

	!! buffer block i-
	kdim = ke0-ks0+1
	idim = is0-is +1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = is ,is0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k), blk(mm)%z(i,j,k), &
			            var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
		end do
	end if

	!! buffer block i+
	kdim = ke0-ks0+1
	idim = ie -ie0+1
	jdim = je0-js0+1
	if (idim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ks0,ke0
		do j = js0,je0
		do i = ie0,ie
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			            var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
		end do
	end if

	!! buffer block j-
	kdim = ke0-ks0+1
	idim = ie0-is0+1
	jdim = js0-js +1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ks0,ke0
		do j = js ,js0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)
		end do
		end do
		end do
	end if

	!! buffer block j+
	kdim = ke0-ks0+1
	idim = ie0-is0+1
	jdim = je -je0+1
	if (jdim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ks0,ke0
		do j = je0,je
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
		end do
	end if	

	!! buffer block k-
	idim = ie0-is0+1
	jdim = je0-je0+1
	kdim = ks0-ks +1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ks, ks0
		do j = js0,je0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			            var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
		end do
	end if	

	!! buffer block k+
	idim = ie0-is0+1
	jdim = je0-js0+1
	kdim = ke -ke0+1
	if (kdim .gt. 1) then
		write(99,*) "variables = x,y,z,var1,var2,var3,var4,var5"
		write(99,*) "zone i= ", idim, "j= ", jdim, "k= ", kdim,"datapacking=point"
			
		do k = ke0,ke
		do j = js0,je0
		do i = is0,ie0
			write(99,*) blk(mm)%x(i,j,k),blk(mm)%y(i,j,k),blk(mm)%z(i,j,k), &
			var(1,i,j,k),var(2,i,j,k),var(3,i,j,k),var(4,i,j,k),var(5,i,j,k)	
		end do
		end do
		end do
	end if	

	close(99)
	return
end subroutine