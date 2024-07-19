      subroutine bottomboundary(txt_loc)
!______________________________________________________________________
!
! S model
! release 2010.18  - last update: 19-02-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
      character*11 txt_loc
!.....................................................................
! Version date      Description des modifications
! 2010.15 30-12-10  Mise en service
! 2010.15 31-12-10  Evolution de la grille sigma-step (transition continue des marches)
!         12-01-10  suite du point precedent: ajout vitesses de Stokes
! 2010.18 19-02-11  Initialiser kbbc_
!.....................................................................
#ifdef synopsis
       subroutinetitle='bottomboundary'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Under Bottom boundary conditions ensuring a smooth connection of the
! bottom levels in case of grid-steps
      ksecu=0

      if(txt_loc=='velocity   ') then !uvuvuvuvuvuvuvuvuv>
      ksecu=1

! Under the bottom and at the location of a "grid step" the current is taken to
! be that of the first level over the bottom using a detph-averaged
! procedure ensuring that the transport of the modified cells remains
! unchanged:
      do j=1,jmax
      do i=1,imax+1
       sum1=0.
       sum2=0.
       do k=kbbc_u(i,j),kmin_u(i,j)
        sum1=sum1+dz_u(i,j,k,2)
        sum2=sum2+dz_u(i,j,k,2)*vel_u(i,j,k,2)
       enddo
       x1=sum2/sum1
       do k=kbbc_u(i,j),kmin_u(i,j)
        vel_u(i,j,k,2)=x1
       enddo
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       sum1=0.
       sum2=0.
       do k=kbbc_v(i,j),kmin_v(i,j)
        sum1=sum1+dz_v(i,j,k,2)
        sum2=sum2+dz_v(i,j,k,2)*vel_v(i,j,k,2)
       enddo
       x1=sum2/sum1
       do k=kbbc_v(i,j),kmin_v(i,j)
        vel_v(i,j,k,2)=x1
       enddo
      enddo
      enddo

      return
      endif                  !uvuvuvuvuvuvuvuvuv>


      if(txt_loc=='temperature') then !tttttttttttttttttt>
      ksecu=1
! Under the bottom and at the location of a "grid step" T (or S) is taken to
! be that of the first level over the bottom using a detph-averaged
! procedure ensuring that the mean T (or S) of the modified cells remains
! unchanged:

      do j=1,jmax
      do i=1,imax

       sum1=0.
       sum2=0.
       do k=kbbc_w(i,j),kmin_w(i,j)-1
        sum1=sum1+dz_t(i,j,k,2)
        sum2=sum2+dz_t(i,j,k,2)*tridia_in(i,j,k,4)  ! tem(2) sous le fond
       enddo

       x1=(sum2+dz_t(i,j,kmin_w(i,j),2)*tem_t(i,j,kmin_w(i,j),2))  &
         /(sum1+dz_t(i,j,kmin_w(i,j),2)                         )

!      do k=kbbc_w(i,j),kmin_w(i,j)
       do k=1,kmin_w(i,j) ! starts from 1 (and not kbbc) for the tracer advection smoothing procedure
        tem_t(i,j,k,2)=x1
       enddo

      enddo
      enddo

      return
      endif                  !ttttttttttttttttttt>

      if(txt_loc=='salinity   ') then !ssssssssssssssssss>
      ksecu=1

      do j=1,jmax
      do i=1,imax

       sum1=0.
       sum2=0.
       do k=kbbc_w(i,j),kmin_w(i,j)-1
        sum1=sum1+dz_t(i,j,k,2)
        sum2=sum2+dz_t(i,j,k,2)*tridia_in(i,j,k,4)  ! tem(2) sous le fond
       enddo

       x1=(sum2+dz_t(i,j,kmin_w(i,j),2)*sal_t(i,j,kmin_w(i,j),2))  &
         /(sum1+dz_t(i,j,kmin_w(i,j),2)                         )

!      do k=kbbc_w(i,j),kmin_w(i,j)
       do k=1,kmin_w(i,j) ! starts from 1 (and not kbbc) for the tracer advection smoothing procedure
        sal_t(i,j,k,2)=x1
       enddo

      enddo
      enddo

      return
      endif                  !ssssssssssssssssss>



#ifdef stokes
      if(txt_loc=='velstokes  ') then !vsvsvsvsvsvsvs> !12-01-11
      ksecu=1

! Under the bottom and at the location of a "grid step" the current is taken to
! be that of the first level over the bottom using a detph-averaged
! procedure ensuring that the transport of the modified cells remains
! unchanged:
      do j=1,jmax
      do i=1,imax+1
       sum1=0.
       sum2=0.
       do k=kbbc_u(i,j),kmin_u(i,j)
        sum1=sum1+dz_u(i,j,k,1)
        sum2=sum2+dz_u(i,j,k,1)*velstokes_u(i,j,k,1)
       enddo
       x1=sum2/sum1
       do k=kbbc_u(i,j),kmin_u(i,j)
        velstokes_u(i,j,k,1)=x1
       enddo
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       sum1=0.
       sum2=0.
       do k=kbbc_v(i,j),kmin_v(i,j)
        sum1=sum1+dz_v(i,j,k,1)
        sum2=sum2+dz_v(i,j,k,1)*velstokes_v(i,j,k,1)
       enddo
       x1=sum2/sum1
       do k=kbbc_v(i,j),kmin_v(i,j)
        velstokes_v(i,j,k,1)=x1
       enddo
      enddo
      enddo

      return
      endif                  !vsvsvsvsvsvsvsvs>
#endif

      if(txt_loc=='init kbbc  ') then !ininininininininin>     !19-02-11
      ksecu=1

! Initializing kbbc_:

! Debut d'indice vertical pour la couche limite "sous le fond": !31-12-10
      do j=1,jmax
      do i=1,imax+1
      kbbc_u(i,j)=kmin_u(i,j)
       do k=kmin_u(i,j)-1,1,-1
       if(mask_t(i,j,k)/=mask_t(i-1,j,k))kbbc_u(i,j)=k
       enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
      kbbc_v(i,j)=kmin_v(i,j)
       do k=kmin_v(i,j)-1,1,-1
       if(mask_t(i,j,k)/=mask_t(i,j-1,k))kbbc_v(i,j)=k
       enddo
      enddo
      enddo
! PATCH POUR UNE AUTRE VISION DE LA GRILLE C: ouvrir carrement le point
! masqué au calcul complet
!     do j=1,jmax
!     do i=1,imax+1
!     if(mask_u(i,j,kmaxp1)==1) then !------>
!      kmin_u(i,j)=kbbc_u(i,j)
!      do k=kmin_u(i,j),kmaxp1
!       mask_u(i,j,k)=1
!      enddo
!     endif                          !------>
!     enddo
!     enddo
!     do j=1,jmax+1
!     do i=1,imax
!     if(mask_v(i,j,kmaxp1)==1) then !------>
!      kmin_v(i,j)=kbbc_v(i,j)
!      do k=kmin_v(i,j),kmaxp1
!       mask_v(i,j,k)=1
!      enddo
!     endif                          !------>
!     enddo
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,jmax
      do i=1,imax
      kbbc_w(i,j)=kmin_w(i,j)
       do k=kmin_w(i,j)-1,1,-1
        if(mask_t(i,j,k)/=mask_t(i-1,j  ,k))kbbc_w(i,j)=k
        if(mask_t(i,j,k)/=mask_t(i+1,j  ,k))kbbc_w(i,j)=k
        if(mask_t(i,j,k)/=mask_t(i  ,j-1,k))kbbc_w(i,j)=k
        if(mask_t(i,j,k)/=mask_t(i  ,j+1,k))kbbc_w(i,j)=k
       enddo
      enddo
      enddo
      call obc_mixsigstep(1) ! C.L. de type "z1" !13-01-11


      return
      endif                           !ininininininininin>



      write(*,*)
      write(*,*)'Error on the argument passed in routine bottomboundary'
      write(*,'(a,a)') txt_loc,' is not valid'
      stop ' STOP in routine bottomboundary'
      end subroutine bottomboundary
