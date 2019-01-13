!!=========================================================!!
!! acoustic files output                                   !!
!! VERSION:1.0 only ouput the single block cylinder                                                       !!
!! coordination: x,y,z                                     !!
!! block dimension: is,ie,js,je,ks,ke                      !!
!!                                                         !!
!! Author: Chi Zhang                                       !!
!! Date:   2016.11.21
!!
!! 2016.12.05 extended to multiblock, not include the case with more than 4 subface!!
!! 2016.
!!=========================================================!!
subroutine acousticfile_output
    use blk_var
    use mpi_var
    use index_var
    use output_var
    use ns_const
    use flag_var
    implicit none

    integer :: i,j,k,kk,i3,j3
    integer :: nodes,z_dim,numbs
    integer :: idim,jdim,kdim
    real*8 :: z_width,dthe,meter
    real*8 :: Ai,Aj, Ak, A2, xc, yc, zc
    real*8 :: xs(1000), ys(1000), zs(1)
    real*8 :: xx(1000,1000,100),yy(1000,1000,100),zz(1000,1000,100)
    character(len=180) :: fwh_surface_xy
    character(len=180) :: fwh_xy
    integer :: is, ie, js, je, ks, ke !!block dimension
    meter=ob_dist
    z_width=acoustic_length
do m0 = 1,blk_loop
   do ksub = 1,blk(m0)%num_subface
    if(blk(m0)%bc(ksub)%fwh .eq. 1)then
! if(blk(m0)%bc(ksub)%fwh .eq. 1)then
    !!**********************************************!!
    !!*******************2d output******************!!
    !!**********************************************!!
     if      (iflag_dimension .eq. iflag_2d) then
                  !if(blk(m0)%bc(ksub)%blk_t .eq. bc_wall)then

                  k = 1  !2D,two dimension
        !face = blk(m0)%bc(ksub)%face
        !!wall dimention

        is = blk(m0)%bc(ksub)%is
        ie = blk(m0)%bc(ksub)%ie
        js = blk(m0)%bc(ksub)%js
        je = blk(m0)%bc(ksub)%je
        !!*
        idim = ie-is+1
        jdim = je-js+1
        kdim = 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if     (blk(m0)%bc(ksub)%face .eq. 1)then
            i=is
            nodes=jdim
            do j=1,nodes
             do kk=1,z_dim
                xx(i,j,kk) = blk(m0)%x(i,j,1)
                yy(i,j,kk) = blk(m0)%y(i,j,1)
                zz(i,j,kk) = z_width*(dble((kk-1))/dble((51-1))-0.50d0)
             end do
            end do
        else if(blk(m0)%bc(ksub)%face .eq. 2)then
                i=ie
                nodes=jdim
                 do j=1,nodes
                  do kk=1,z_dim
                    xx(i,j,kk) = blk(m0)%x(i,j,1)
                    yy(i,j,kk) = blk(m0)%y(i,j,1)
                    zz(i,j,kk) = z_width*(dble((kk-1))/dble((51-1))-0.50d0)
                  end do
                end do
        else if(blk(m0)%bc(ksub)%face .eq. 3)then
                    j=js
                    nodes=idim
                    do i=1,nodes
                        do kk=1,z_dim
                        xx(i,j,kk) = blk(m0)%x(i,j,1)
                        yy(i,j,kk) = blk(m0)%y(i,j,1)
                        zz(i,j,kk)=z_width*(dble((kk-1))/dble((51-1))-0.50d0)
                        end do
                    end do
        else if(blk(m0)%bc(ksub)%face .eq. 4)then
                j=je
                nodes=idim
                do i=1,nodes
                    do kk=1,z_dim
                    xx(i,j,kk) = blk(m0)%x(i,j,1)
                    yy(i,j,kk) = blk(m0)%y(i,j,1)
                    zz(i,j,kk)=z_width*(dble((kk-1))/dble((51-1))-0.50d0)
                    end do
                end do
         end if
                !!!!!!!!!!!!!!!!!!!!!!!!!
        z_dim=51
    !!**********************************************!!
    !!***************end 2d output******************!!
    !!**********************************************!!

    !!**********************************************!!
    !!***************begin 3d output****************!!
    !!**********************************************!!
     else if(iflag_dimension .eq. iflag_3d) then
       !if(blk(m0)%bc(ksub)%fwh .eq. 1)then
        is = blk(m0)%bc(ksub)%is
        ie = blk(m0)%bc(ksub)%ie
        js = blk(m0)%bc(ksub)%js
        je = blk(m0)%bc(ksub)%je
        ks = blk(m0)%bc(ksub)%ks
        ke = blk(m0)%bc(ksub)%ke
        idim = ie-is+1
        jdim = je-js+1
        kdim = ke-ks+1
        if     (blk(m0)%bc(ksub)%face .eq. 1)then
                i=is
                nodes=jdim
                 do j=js,je
                  do k=ks,ke
                xx(i,j,k) = blk(m0)%x(i,j,k)
                yy(i,j,k) = blk(m0)%y(i,j,k)
                zz(i,j,k) = blk(m0)%z(i,j,k)
                  end do
                 end do
        else if(blk(m0)%bc(ksub)%face .eq. 2)then
                i=ie
                nodes=jdim
                 do j=js,je
                  do k=ks,ke
                xx(i,j,k) = blk(m0)%x(i,j,k)
                yy(i,j,k) = blk(m0)%y(i,j,k)
                zz(i,j,k) = blk(m0)%z(i,j,k)
                  end do
                 end do
        else if(blk(m0)%bc(ksub)%face .eq. 3)then
                j=js
                nodes=idim
                 do i=is,ie
                  do k=ks,ke
                xx(i,j,k) = blk(m0)%x(i,j,k)
                yy(i,j,k) = blk(m0)%y(i,j,k)
                zz(i,j,k) = blk(m0)%z(i,j,k)
                  end do
                 end do
        else if(blk(m0)%bc(ksub)%face .eq. 4)then
                j=je
                nodes=idim
                 do i=is,ie
                  do k=ks,ke
                xx(i,j,k) = blk(m0)%x(i,j,k)
                yy(i,j,k) = blk(m0)%y(i,j,k)
                zz(i,j,k) = blk(m0)%z(i,j,k)
                  end do
                 end do
        end if
         kdim=ke-ks+1
    !!**********************************************!!
    !!***************end 3d output******************!!
    !!**********************************************!!
     end if !dimension judge
        !numbs=nodes-1
     !!!!!!!export the planar geometry of the FW-H surface
    write(fwh_surface_xy,"('result/fwh_surface_xy_'I4.4'_'I4.4'.dat')") myn*myid+m0,blk(m0)%bc(ksub)%face
        open(201,file = fwh_surface_xy)!,status='replace')
                if     (blk(m0)%bc(ksub)%face .eq. 1)then
                    i=is
                        do j=js,je
                        k=1
                           write(201,*) j,xx(i,j,k),yy(i,j,k)
                        end do
                else if(blk(m0)%bc(ksub)%face .eq. 2)then
                    i=ie
                        do j=js,je
                        k=1
                           write(201,*) j,xx(i,j,k),yy(i,j,k)
                        end do
                else if(blk(m0)%bc(ksub)%face .eq. 3)then
                    j=js
                        do i=is,ie
                        k=1
                           write(201,*) i,xx(i,j,k),yy(i,j,k)
                        end do
                else if(blk(m0)%bc(ksub)%face .eq. 4)then
                    j=je
                        do i=is,ie
                        k=1
                           write(201,*) i,xx(i,j,k),yy(i,j,k)
                        end do
                end if
                !!!!!!!!!!!!!!!!!!!!!!!!!
        close(201)
    !!!!!!!
        !dtheta = 360.0d0 / numbs
        !ALLOCATE( xx(nodes,1, z_dim))
        !ALLOCATE( yy(nodes,1, z_dim))
        !ALLOCATE( zz(nodes,1, z_dim))
    !j=js
    !all coordinate
    !do i=1,nodes
    !   do kk=1,z_dim
    !    xx(i,1,kk) = blk(m0)%x(i,j,1)
    !    yy(i,1,kk) = blk(m0)%y(i,j,1)
        !    zz(i,1,kk)=z_width*(dble((kk-1))/dble((51-1))-0.50d0)
    !   end do
    !end do
    !!!!!!!!

    !write(6,*) is,ie,js,je,blk(m0)%x(1,1,1),blk(m0)%y(1,1,1)
    !pause
 !!calculate the normal vector,center point and area
    write(fwh_xy,"('result/fwh_xy_'I4.4'_'I4.4'.dat')") myn*myid+m0,blk(m0)%bc(ksub)%face
    OPEN(132,file=fwh_xy)

    WRITE(132,*) (nodes-1)*(kdim-1)    != panels
                if     (blk(m0)%bc(ksub)%face .eq. 1)then
                    i=is
                    i3=is
                        do j=js,je-1
                        j3=j+1
                        DO k = 1,kdim-1    ! azimuthal dir
                              Ai = -0.5D0*( (yy(i3,j3,k)-yy(i,j,k+1)) &
                                        *(zz(i,j,k) -zz(i3,j3,k+1)) &
                                        -(yy(i,j,k) -yy(i3,j3,k+1)) &
                                        *(zz(i3,j3,k)-zz(i,j,k+1)) )
!write(6,*) yy(i3,j3,k),yy(i,j,k+1),zz(i,j,k),zz(i3,j3,k+1),yy(i,j,k),yy(i3,j3,k+1),zz(i3,j3,k),zz(i,j,k+1)
!pause
                              Aj = 0.5D0*( (xx(i3,j3,k)-xx(i,j,k+1))     &
                                        *(zz(i,j,k) -zz(i3,j3,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j3,k+1))     &
                                        *(zz(i3,j3,k)-zz(i,j,k+1)) )
                              Ak = 0.5D0*( (xx(i3,j3,k)-xx(i,j,k+1))     &
                                        *(yy(i,j,k) -yy(i3,j3,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j3,k+1))     &
                                        *(yy(i3,j3,k)-yy(i,j,k+1)) )
                        !
                        !   Unit normal vetor pointing out of the fluid
                        !
                                A2 = DSQRT(Ai*Ai+Aj*Aj+Ak*Ak)       ! surface area
                                Ai =Ai/A2                           ! i component
                                Aj =Aj/A2                           ! j component
                                Ak =Ak/A2                           ! k component
                        !
                        !    get panel centre  position
                        !
                                xc = 0.25d0*(xx(i,j,k) +xx(i,j,k+1)     &
                                    +xx(i3,j,k)+xx(i3,j,k+1))
                                yc = 0.25d0*(yy(i,j,k) +yy(i,j,k+1)     &
                                    +yy(i3,j,k)+yy(i3,j,k+1))
                                zc = 0.25d0*(zz(i,j,k) +zz(i,j,k+1)     &
                                    +zz(i3,j,k)+zz(i3,j,k+1))
                              WRITE(132,'(4f12.8,4f13.7)') Ai,Aj,Ak,A2, xc,yc,zc,1.0d0
                              END DO
                            END DO
                            close(132)
                else if(blk(m0)%bc(ksub)%face .eq. 2)then
                    i=ie
                    i3=ie
                        do j=js,je-1
                        j3=j+1
                        DO k = 1,kdim-1    ! azimuthal dir
                              Ai = -0.5D0*( (yy(i3,j3,k)-yy(i,j,k+1)) &
                                        *(zz(i,j,k) -zz(i3,j3,k+1)) &
                                        -(yy(i,j,k) -yy(i3,j3,k+1)) &
                                        *(zz(i3,j3,k)-zz(i,j,k+1)) )
                              Aj = 0.5D0*( (xx(i3,j3,k)-xx(i,j,k+1))     &
                                        *(zz(i,j,k) -zz(i3,j3,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j3,k+1))     &
                                        *(zz(i3,j3,k)-zz(i,j,k+1)) )
                              Ak = 0.5D0*( (xx(i3,j3,k)-xx(i,j,k+1))     &
                                        *(yy(i,j,k) -yy(i3,j3,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j3,k+1))     &
                                        *(yy(i3,j3,k)-yy(i,j,k+1)) )
                        !
                        !   Unit normal vetor pointing out of the fluid
                        !
                                A2 = DSQRT(Ai*Ai+Aj*Aj+Ak*Ak)       ! surface area
                                Ai =Ai/A2                           ! i component
                                Aj =Aj/A2                           ! j component
                                Ak =Ak/A2                           ! k component
                        !
                        !    get panel centre  position
                        !
                                xc = 0.25d0*(xx(i,j,k) +xx(i,j,k+1)     &
                                    +xx(i3,j,k)+xx(i3,j,k+1))
                                yc = 0.25d0*(yy(i,j,k) +yy(i,j,k+1)     &
                                    +yy(i3,j,k)+yy(i3,j,k+1))
                                zc = 0.25d0*(zz(i,j,k) +zz(i,j,k+1)     &
                                    +zz(i3,j,k)+zz(i3,j,k+1))
                              WRITE(132,'(4f12.8,4f13.7)') Ai,Aj,Ak,A2, xc,yc,zc,1.0d0
                              END DO
                            END DO
                            close(132)
                else if(blk(m0)%bc(ksub)%face .eq. 3)then
                    j=js
                        do i=is,ie-1
                        i3=i+1
                        DO k = 1,kdim-1    ! azimuthal dir
                              Ai = -0.5D0*( (yy(i3,j,k)-yy(i,j,k+1)) &
                                        *(zz(i,j,k) -zz(i3,j,k+1)) &
                                        -(yy(i,j,k) -yy(i3,j,k+1)) &
                                        *(zz(i3,j,k)-zz(i,j,k+1)) )
                              Aj = 0.5D0*( (xx(i3,j,k)-xx(i,j,k+1))     &
                                        *(zz(i,j,k) -zz(i3,j,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j,k+1))     &
                                        *(zz(i3,j,k)-zz(i,j,k+1)) )
                              Ak = 0.5D0*( (xx(i3,j,k)-xx(i,j,k+1))     &
                                        *(yy(i,j,k) -yy(i3,j,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j,k+1))     &
                                        *(yy(i3,j,k)-yy(i,j,k+1)) )
                        !write(6,*) AI,AJ,AK
                        !pause
                        !   Unit normal vetor pointing out of the fluid
                        !
                                A2 = DSQRT(Ai*Ai+Aj*Aj+Ak*Ak)       ! surface area
                                Ai =Ai/A2                           ! i component
                                Aj =Aj/A2                           ! j component
                                Ak =Ak/A2                           ! k component
                        !
                        !    get panel centre  position
                        !
                                xc = 0.25d0*(xx(i,j,k) +xx(i,j,k+1)     &
                                    +xx(i3,j,k)+xx(i3,j,k+1))
                                yc = 0.25d0*(yy(i,j,k) +yy(i,j,k+1)     &
                                    +yy(i3,j,k)+yy(i3,j,k+1))
                                zc = 0.25d0*(zz(i,j,k) +zz(i,j,k+1)     &
                                    +zz(i3,j,k)+zz(i3,j,k+1))
                              WRITE(132,'(4f12.8,4f13.7)') Ai,Aj,Ak,A2, xc,yc,zc,1.0d0
                              END DO
                            END DO
                            close(132)
                else if(blk(m0)%bc(ksub)%face .eq. 4)then
                    j=je
                        do i=is,ie-1
                        i3=i+1
                        DO k = 1,kdim-1    ! azimuthal dir
                              Ai = -0.5D0*( (yy(i3,j,k)-yy(i,j,k+1)) &
                                        *(zz(i,j,k) -zz(i3,j,k+1)) &
                                        -(yy(i,j,k) -yy(i3,j,k+1)) &
                                        *(zz(i3,j,k)-zz(i,j,k+1)) )
                              Aj = 0.5D0*( (xx(i3,j,k)-xx(i,j,k+1))     &
                                        *(zz(i,j,k) -zz(i3,j,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j,k+1))     &
                                        *(zz(i3,j,k)-zz(i,j,k+1)) )
                              Ak = 0.5D0*( (xx(i3,j,k)-xx(i,j,k+1))     &
                                        *(yy(i,j,k) -yy(i3,j,k+1))     &
                                        -(xx(i,j,k) -xx(i3,j,k+1))     &
                                        *(yy(i3,j,k)-yy(i,j,k+1)) )
                        !
                        !   Unit normal vetor pointing out of the fluid
                        !
                                A2 = DSQRT(Ai*Ai+Aj*Aj+Ak*Ak)       ! surface area
                                Ai =Ai/A2                           ! i component
                                Aj =Aj/A2                           ! j component
                                Ak =Ak/A2                           ! k component
                        !
                        !    get panel centre  position
                        !
                                xc = 0.25d0*(xx(i,j,k) +xx(i,j,k+1)     &
                                    +xx(i3,j,k)+xx(i3,j,k+1))
                                yc = 0.25d0*(yy(i,j,k) +yy(i,j,k+1)     &
                                    +yy(i3,j,k)+yy(i3,j,k+1))
                                zc = 0.25d0*(zz(i,j,k) +zz(i,j,k+1)     &
                                    +zz(i3,j,k)+zz(i3,j,k+1))
                              WRITE(132,'(4f12.8,4f13.7)') Ai,Aj,Ak,A2, xc,yc,zc,1.0d0
                              END DO
                            END DO
                            close(132)
                                end if



!----------------------------------------------------------------------
!   Allocate observer positions
!----------------------------------------------------------------------
!

        !deallocate(xx)
        !deallocate(yy)
        !deallocate(zz)

    end if !fwh
   end do !subfacce loop
end do !blk_loop

    if(myid .eq. root)then
    open(111,file='result/obs_xyzdata.dat')

    !ALLOCATE( xs(nobs), ys(nobs), zs(1) )
    zs   = 0.0d0
        write(6,*) 'Observor from 0 to 360 degrees'
        dthe = 360.0d0/float(nobs-1)*pi/180.0d0
        write(6,'("Using a center point at x=0, y=0")')
          DO i = 1,nobs
            xs(i) =meter*dcos(float(i-1)*dthe)+0.0d0
            ys(i) = meter*dsin(float(i-1)*dthe)+0.0d0
          END DO
      write(111,'(i7,i2,"  ! nobs")') nobs
        DO i = 1,nobs
          write(111,'(i5,3f11.5)') i,xs(i),ys(i),zs
          END DO
    CLOSE(111)
    end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !deallocate(xs)
    !deallocate(ys)
    !deallocate(zs)
    return
end subroutine
