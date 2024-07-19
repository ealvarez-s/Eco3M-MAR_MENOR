      module module_curvgrdtoolbox
!______________________________________________________________________
! SYMPHONIE ocean model
! release 301 - last update: 26-04-21
!______________________________________________________________________
! Version date      Description des modifications
! S26     25-02-14  creation
!         20-03-14  ajout d'un coef de diffusion dependant de (i,j) dans
!                   les routines smooth...
!         31-07-14  Nouveaux echanges
!         10-06-15  nouvelles fonctionnalites
!         04-06-16  Correction dans bug dans la fonction de lissage
!         15-02-18  correction diag dy_v
! v295    01-02-21  mises A jours diverses
! v301    26-04-21  debug du cas periodique
!______________________________________________________________________
!    _________                    .__                  .__             ! [°v°] Hello!
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................
      use module_principal
      use module_parallele !#mpi
      implicit none
      include 'netcdf.inc'
      double precision , dimension(:,:) , allocatable ::  &
       var_j,var_i,vlx_u,vly_v,pos_x,pos_y,pos_lon,pos_lat
      integer, dimension(:,:) , allocatable :: mask_curvi

contains

      subroutine curvgrdtoolbox_driver(case_) !01-01-14
      implicit none
      double precision dt_,deci_b,decj_b,indice_j_b &
                          ,deci_a,decj_a,indice_j_a &
                          ,vlu_,vlv_,filval_,dlon_di_ &
                          ,dlon_dj_,dist_,lat1_,lat2_,lon2_,lon1_ &
                          ,deci_bb,decj_bb
      integer :: ncid_,dim_x_id_,dim_y_id_,len_,case_,jbase_=-999,ibase_=-999
#ifdef synopsis
       subroutinetitle='curvgrdtoolbox_driver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef bidon
      if(nbdom/=1) then
       write(6,*)'Erreur: Generateur de grille curvgrdtoolbox', &
                 ' ne peut se faire qu''en monoproc'
       stop ' STOP dans curvgrdtoolbox'
      endif

! Numeroter les iles
      do j=0,jmax+1
      do i=0,imax+1
       xy_t(i,j,0)=mask_t(i,j,kmaxp1)
      enddo
      enddo

      k=1
      do j1=0,jmax+1
       write(6,*)'%iles=',100.*real(j1)/real(jmax+1)
      do i1=0,imax+1

      if(nint(xy_t(i1,j1,0))==0)then !>>>>>>>>>>>>

        k=k+1
        xy_t(i1,j1,0)=k

        sum1=0.
 1234   continue

        do j=0,jmax+1
        do i=0,imax+1
         if(nint(xy_t(i,j,0))==k)then !---->
           if(nint(xy_t(i+1,j  ,0))==0)xy_t(i+1,j  ,0)=k
           if(nint(xy_t(i-1,j  ,0))==0)xy_t(i-1,j  ,0)=k
           if(nint(xy_t(i  ,j+1,0))==0)xy_t(i  ,j+1,0)=k
           if(nint(xy_t(i  ,j-1,0))==0)xy_t(i  ,j-1,0)=k
         endif                        !---->
        enddo
        enddo

        sum2=sum1
        sum1=0.
        do j=0,jmax+1
        do i=0,imax+1
         if(nint(xy_t(i,j,0))==k)sum1=sum1+1.
        enddo
        enddo

        if(sum1/=sum2)goto 1234

      endif                          !>>>>>>>>>>>>

!     write(6,*)i1,j1,nint(xy_t(i1,j1,0))
      enddo
      enddo
#endif


!     jbase_=-999
!     jbase_=1
      jbase_=jmax ! phase 1 puis phase 4

!     ibase_=-999
!     ibase_=1    ! phase 2
!     ibase_=imax ! phase 3 puis phase 5

      len_=1
      filval_=-999.
      if(.not.allocated(var_j))then
       allocate(var_j(0:imax+1,0:jmax+1)) ; var_j=0
      endif
      if(.not.allocated(var_i))then
       allocate(var_i(0:imax+1,0:jmax+1)) ; var_i=0
      endif
      if(.not.allocated(mask_curvi)) then
       allocate(mask_curvi(0:imax+1,0:jmax+1)) ; mask_curvi=0
      endif

      if(case_<=1) then !1111111111111>

      if(case_==0)goto 1479 ! Aller directement au fichier netcdf pour voir les numeros des iles

      do j=0,jmax+1
      do i=0,imax+1
            var_j(i,j)=real(j)
            var_i(i,j)=real(i)
       mask_curvi(i,j)=1
!      if(mask_t(i,j,kmax)==0)mask_curvi(i,j)=int(filval_)
      enddo
      enddo

! Step 1:
      if(jbase_==1) then !-------->

       do j=0,jmax+1
       do i=0,imax+1

! Grille 1
!        var_j(i,j)=real(j)+20.*( (1.+tanh(real(i-240)/70. ))*0.5  &
!                                +(1.+tanh(real(670-i)/120.))*0.5) &
!                              *(1-tanh(real(jglb-600-j)/22.))
! Grille 2 plus etendue au large
         var_j(i,j)=real(j)+40.*( (1.+tanh(real(i-200)/70. ))*0.5  &
                                 +(1.+tanh(real(630-i)/120.))*0.5) &
                               *(1-tanh(real(jglb-600-j)/22.))

       enddo
       enddo
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1,3000,'mask')
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1,1000,'invm')
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1, 100,'free')

! Le bord de grille = 1 etape de grd nul pour eviter sortie de domaine +
! surcouche de lissage pour lisser la discontinuite ainsi introduite
!      call curvgrdtoolbox_obc_varj(20,20) ! largeur zone grd nul bords i=1 i=imax
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1,20,'free')


      endif                 !---------->

! Step 4:
      if(ibase_==imax) then !-------->

! xy_t(i,j,1) mesure la distance depuis i=i10=667 pour une valeur de j fixee en j10=454
       i10=667 ; j10=454
       j=j10 
       do i=imax  ,0,-1
        if(i>i10) then 
         xy_t(i,j,1)=0.
        else
         xy_t(i,j,1)=xy_t(i+1,j,1)+dx_u(i+1,j)
        endif
!       if(j==j10)write(666,*)i,xy_t(i,j,1)
       enddo 

       do j=0,jmax+1
       do i=0,imax+1
        mask_curvi(i,j)=1

! STEP 3
! Pour s'assurer que var_i(imax,:) est bien >=imax, on soustrait une constante qui
! permet d'avoir var_i(imax,:)>=imax (c'est A ca que sert i1 en suivant)
!       i1=imax-1
!       if(j<88) then
!        x3=40
!       else
!        x3=80
!       endif
!       var_i(i,j)=real(i) &
!          -30*exp(-(real(j-88)/x3)**2 )*0.5*(1-tanh(real(i -194)/15.))  &
!          +30*exp(-(real(j-88)/x3)**2 )*0.5*(1-tanh(real(i1-194)/15.))  
! STEP 5
!       var_i(i,j)=real(i)+30*( 0.5*(1-tanh(real(i -150)/15.)) &
!                              )*exp(-(real(j-1)/35)**2 )
        if(i>=i10+1) then
          var_i(i,j)=real(i) 
           if(j==j10)write(667,*)i,var_i(i,j),' A'
        else
          var_i(i,j)=real(i10+1)-xy_t(i,j10,1)/6.
          if(j==j10)write(667,*)i,var_i(i,j),real(xy_t(i,j10,1))
        endif


       enddo
       enddo
!      stop 'bravissimo'
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1,1000,'mask')

!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1,3000,'mask')
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1,1000,'invm')
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1, 100,'free')
! Le bord de grille = 1 etape de grd nul pour eviter sortie de domaine +
! surcouche de lissage pour lisser la discontinuite ainsi introduite
!      call curvgrdtoolbox_obc_vari(20,20) ! largeur zone grd nul bords j=1 j=imax
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1, 20,'free')
      endif                 !---------->

! Step 3:
      if(ibase_==1) then !-------->

       do j=0,jmax+1
       do i=0,imax+1


!       mask_curvi(i,j)=1

         var_i(i,j)=real(i)  &
                           +40.*exp(-( real(j-450)/90.)**2) &
                                *(1.-tanh(real(iglb-100-i)/20.))

       enddo
       enddo
!      stop 'bravo'
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1,1000,'mask')
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1, 50,'free')

! Le bord de grille = 1 etape de grd nul pour eviter sortie de domaine +
! surcouche de lissage pour lisser la discontinuite ainsi introduite
!      call curvgrdtoolbox_obc_vari(20,20) ! largeur zone grd nul bords j=1 j=imax
!      call curvgrdtoolbox_smooth_vari(2,imax-1,0,jmax+1, 30,'free')

      endif                 !---------->

! Step 2:
      if(jbase_==jmax) then !-------->
       do j=0,jmax+1
       do i=0,imax+1
!        mask_curvi(i,j)=1

! grille 1
!         var_j(i,j)=real(j)-100.      &
!                                *0.5*(1.-tanh(real(i-160)/50.)) &
!                                *0.5*(1.-tanh(real(j-500)/60.))
! grille 2
          var_j(i,j)=real(j)-10.*exp(-(real(i-713)/90.)**2) &
                                *(1.-tanh(real(j-570)/30.))*0.5

       enddo
       enddo
! STEP 4:
!      call curvgrdtoolbox_smooth_varj(0,imax+1,80,110, 20,'free') ! STEP4 SEULEMENT
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1,1000,'mask')

! Le bord de grille = 1 etape de grd nul pour eviter sortie de domaine +
! surcouche de lissage pour lisser la discontinuite ainsi introduite
!      call curvgrdtoolbox_obc_varj(20,20) ! largeur zone grd nul bords i=1 i=imax
!      call curvgrdtoolbox_smooth_varj(0,imax+1,2,jmax-1,30,'free')
      endif                 !---------->


#ifdef bidon
! Empecher les isolignes de sortir des bords lateraux avec une condition de gradient nul:
      if(jbase_/=-999) then !ooo>
       var_j(0     ,:)=var_j(2     ,:)
       var_j(1     ,:)=var_j(2     ,:)
       var_j(imax+1,:)=var_j(imax-1,:)
       var_j(imax  ,:)=var_j(imax-1,:)
      endif                 !ooo>
      if(ibase_/=-999) then !ooo>
       var_i(:,0)     =var_i(:,2)
       var_i(:,1)     =var_i(:,2)
       var_i(:,jmax+1)=var_i(:,jmax-1)
       var_i(:,jmax  )=var_i(:,jmax-1)
      endif                 !ooo>
#endif


! Calcul des angles des axes de la grille:
      xy_t(:,:,2)=90.
      do j=2,jmax-1
      do i=2,imax-1

! Angle1
       dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
       xy_t(i,j,2)=atan2(lat_t(i+1,j)-lat_t(i-1,j)  &
                        ,dlon_di_*cos(lat_t(i,j)))

! Angle2-angle1
       dlon_dj_=lon_t(i,j+1)-lon_t(i,j-1)
       if(dlon_dj_<-pi)dlon_dj_=dlon_dj_+2.*pi
       if(dlon_dj_> pi)dlon_dj_=dlon_dj_-2.*pi
       xy_t(i,j,2)=-xy_t(i,j,2)                          &
                  +atan2(lat_t(i,j+1)-lat_t(i,j-1)   &
                        ,dlon_dj_*cos(lat_t(i,j)))

! Conversion degres:
       xy_t(i,j,2)=xy_t(i,j,2)*rad2deg
       if(xy_t(i,j,2)<-180.)xy_t(i,j,2)=xy_t(i,j,2)+360.

      enddo
      enddo
! CL
      do j=2,jmax-1
       xy_t(imax,j,2)=xy_t(imax-1,j,2)
       xy_t(1   ,j,2)=xy_t(2     ,j,2)
      enddo
      do i=1,imax
       xy_t(i,jmax,2)=xy_t(i,jmax-1,2)
       xy_t(i,1   ,2)=xy_t(i,2     ,2)
      enddo

 1479 continue

! Ecrire fichier netcdf
      status=nf_create('grille_curvi.nc',nf_clobber,ncid_)
      status=nf_def_dim(ncid_,'ni_t',imax+2,dim_x_id_)
      status=nf_def_dim(ncid_,'nj_t',jmax+2,dim_y_id_)

      vardim(1)=dim_x_id_ ; vardim(2)=dim_y_id_

      status=nf_def_var(ncid_,'longitude',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'latitude',nf_double,2,vardim(1:2),var_id)

      status=nf_def_var(ncid_,'var_j',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'var_i',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'var_j_mask',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'var_i_mask',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'deltavar_j',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'deltavar_i',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_def_var(ncid_,'angle',nf_double,2,vardim(1:2),var_id)

      status=nf_def_var(ncid_,'mask',nf_int,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'mask_curvi',nf_int,2,vardim(1:2),var_id)
      status=nf_put_att_int(ncid_,var_id,'_FillValue' &
                               ,nf_int,len_,int(filval_))

      status=nf_def_var(ncid_,'lands',nf_int,2,vardim(1:2),var_id)

      status=nf_enddef(ncid_)

      status=nf_inq_varid(ncid_,'longitude',var_id)
      if(status/=0)stop 'erreur nf_inq_varid longitude'
      anyv3d(0:imax+1,0:jmax+1,1,1)=lon_t(0:imax+1,0:jmax+1)*rad2deg
      status=nf_put_var_double(ncid_,var_id,anyv3d(0:imax+1,0:jmax+1,1,1))

      status=nf_inq_varid(ncid_,'latitude',var_id)
      if(status/=0)stop 'erreur nf_inq_varid latitude'
      anyv3d(0:imax+1,0:jmax+1,1,1)=lat_t(0:imax+1,0:jmax+1)*rad2deg
      status=nf_put_var_double(ncid_,var_id,anyv3d(0:imax+1,0:jmax+1,1,1))

      status=nf_inq_varid(ncid_,'var_j',var_id)
      if(status/=0)stop 'erreur nf_inq_varid var_j'
      status=nf_put_var_double(ncid_,var_id,var_j(0:imax+1,0:jmax+1))

      status=nf_inq_varid(ncid_,'deltavar_j',var_id)
      if(status/=0)stop 'erreur nf_inq_varid deltavar_j'
      do i=0,imax+1 ; do j=0,jmax+1
       anyv3d(i,j,1,1)=var_j(i,j)-j
      enddo         ; enddo
      status=nf_put_var_double(ncid_,var_id,anyv3d(0:imax+1,0:jmax+1,1,1))

      status=nf_inq_varid(ncid_,'var_i',var_id)
      if(status/=0)stop 'erreur nf_inq_varid var_i'
      status=nf_put_var_double(ncid_,var_id,var_i(0:imax+1,0:jmax+1))

      status=nf_inq_varid(ncid_,'deltavar_i',var_id)
      if(status/=0)stop 'erreur nf_inq_varid deltavar_i'
      do i=0,imax+1 ; do j=0,jmax+1
       anyv3d(i,j,1,1)=var_i(i,j)-i
      enddo         ; enddo
      status=nf_put_var_double(ncid_,var_id,anyv3d(0:imax+1,0:jmax+1,1,1))

      status=nf_inq_varid(ncid_,'var_j_mask',var_id)
      if(status/=0)stop 'erreur nf_inq_varid longitude_t'
      status=nf_put_var_double(ncid_,var_id,                     &
       var_j(0:imax+1,0:jmax+1)*mask_t(0:imax+1,0:jmax+1,kmaxp1) &
       +filval_*(1.-mask_t(0:imax+1,0:jmax+1,kmaxp1)))

      status=nf_inq_varid(ncid_,'var_i_mask',var_id)
      if(status/=0)stop 'erreur nf_inq_varid longitude_t'
      status=nf_put_var_double(ncid_,var_id,                     &
       var_i(0:imax+1,0:jmax+1)*mask_t(0:imax+1,0:jmax+1,kmaxp1) &
       +filval_*(1.-mask_t(0:imax+1,0:jmax+1,kmaxp1)))

      status=nf_inq_varid(ncid_,'angle',var_id)
      if(status/=0)stop 'erreur nf_inq_varid angle'
      status=nf_put_var_double(ncid_,var_id,xy_t(0:imax+1,0:jmax+1,2))

      status=nf_inq_varid(ncid_,'mask',var_id)
      status=nf_put_var_int1(ncid_,var_id,mask_t(0:imax+1,0:jmax+1,kmax))

      status=nf_inq_varid(ncid_,'mask_curvi',var_id)
      status=nf_put_var_int(ncid_,var_id,mask_curvi(0:imax+1,0:jmax+1))

      status=nf_inq_varid(ncid_,'lands',var_id)
      status=nf_put_var_int(ncid_,var_id,nint(xy_t(0:imax+1,0:jmax+1,0)))

      status=nf_close(ncid_)

!     stop 'Miles'
      endif             !1111111111111>

      if(case_==2) then !2222222222222>

       allocate(vlx_u  (1:imax+1,0:jmax+1))
       allocate(vly_v  (0:imax+1,1:jmax+1))
       allocate(pos_x  (imax  ,jmax))
       allocate(pos_y  (imax  ,jmax))
       allocate(pos_lon(imax  ,jmax))
       allocate(pos_lat(imax  ,jmax))

       status=nf_open('grille_curvi.nc',nf_nowrite,ncid_)
       status=nf_inq_varid(ncid_,'var_i',var_id)
       if(status/=0)stop 'erreur 1'
       status=nf_get_var_double(ncid_,var_id,var_i(0:imax+1,0:jmax+1))
       if(status/=0)stop 'erreur 2'
       status=nf_inq_varid(ncid_,'var_j',var_id)
       if(status/=0)stop 'erreur 3'
       status=nf_get_var_double(ncid_,var_id,var_j(0:imax+1,0:jmax+1))
       if(status/=0)stop 'erreur 4'
       status=nf_close(ncid_)

      if(jbase_/=-999) then !2222222>
      do j=0,jmax+1 ; do i=1,imax+1
       j1=min(max(j,1),jmax)
       vlx_u(i,j)=(var_j(i,j)     -var_j(i-1,j))/max(dx_u(i,j1)**2,small1)
!      if(i==1)write(6,*)'i,j,vlx_u(i,j)=',i,j,vlx_u(i,j)
      enddo ; enddo
      do j=1,jmax+1 ; do i=0,imax+1
       i1=min(max(i,1),imax)
       vly_v(i,j)=(var_j(i,j)-var_j(i,j-1))/(dy_v(i1,j)**2)
!      if(mask_v(i,j,kmaxp1)==0) then
!       vly_v(i,j)=max(var_j(i,j)-var_j(i,j-1),0.001d00)/max(dy_v(i1,j)**2,small1)
!      endif
      enddo ; enddo
      endif                 !2222222>
      if(ibase_/=-999) then !3333333>
      do j=0,jmax+1 ; do i=1,imax+1
       j1=min(max(j,1),jmax)
       vlx_u(i,j)=(var_i(i,j)-var_i(i-1,j))/(dx_u(i,j1)**2)
!      if(mask_u(i,j,kmaxp1)==0) then
!       vlx_u(i,j)=max(var_i(i,j)-var_i(i-1,j),0.001d00)/max(dx_u(i,j1)**2,small1)
!      endif
      enddo ; enddo
      do j=1,jmax+1 ; do i=0,imax+1
       i1=min(max(i,1),imax)
       vly_v(i,j)=(var_i(i,j)-var_i(i,j-1))/max(dy_v(i1,j)**2,small1)
      enddo ; enddo

      endif                 !3333333>

      if(jbase_/=-999) then !2222222>
! Lignes iso-i:
      dt_=1.
      do i2=1,imax
       write(6,'(a,i0,a,i0)')'-----imax=',imax,' et ----i2=',i2 !01-02-21
       call curvgrdtoolbox_isoline_i(jbase_,i2*un)
      enddo ! i2
      write(6,*)'Fin de boucle i2'
      endif                 !2222222>

      if(ibase_/=-999) then !3333333>
! Lignes iso-j:
       dt_=1.
       do j2=1,jmax
       write(6,'(a,i0,a,i0)')'-----jmax=',jmax,' et ----j2=',j2 !01-02-21
        call curvgrdtoolbox_isoline_j(ibase_,j2*un)
       enddo
      endif                !3333333>


#ifdef bidon
      open(unit=3,file=trim(tmpdirname)//'trous_grille_curvi')

      i0=-999 ; i1=-999
      read(3,*,end=1532)i2 ; write(6,*)'Trou en i2=',i2 ; i0=i2 ! premier trou
 1531 read(3,*,end=1532)i2 ; write(6,*)'Trou en i2=',i2 ; i1=i2 ! dernier trou
      goto 1531
 1532 close(3)
      if(i0/=-999) then !//////>

! Ici nous choisissons d'elargir un peu plus la zone trouee car on
! on considere que ses bords ont eux memes etes contamines par la zone
! masquee
       write(6,*)'Premier et dernier trou=',i0,i1
       if(jbase_==1   )j2=jmax
       if(jbase_==jmax)j2=1
!      write(6,*)'j2=',j2
       write(6,*)'Positions des fins de trajectoires ',&
                 'valides de part et d''autre:'
       write(6,*)'i et j du point avant ',pos_x(i0-1,j2),pos_y(i0-1,j2)
       write(6,*)'i et j du point apres ',pos_x(i1+1,j2),pos_y(i1+1,j2)
       write(6,*)'Points de depart des trajectoires inverses:'
       open(unit=3,file=trim(tmpdirname)//'departs_inverses')
       do i=i0,i1
       k0=i-i0+2
       k1=i1-i+2
       rap=real(i-i0+1)/real(i1-i0+2)
!      rap=0.5*(1.+sin((rap-0.5)*pi))
!      write(3,*)i,j2,(1.-rap)*pos_x(i0-1,j2)+rap*pos_x(i1+1,j2)
       write(3,*)i,j2,                                          &
           (1.-rap)*( k0*pos_x(i0-1,j2)-(k0-1)*pos_x(i0-2,j2) ) &
          +    rap *( k1*pos_x(i1+1,j2)-(k1-1)*pos_x(i1+2,j2) )
       enddo
       close(3)
! Calcul trajectoires inverses:
       if(jbase_==jmax) then ; jbase_=1 ; else ; jbase_=jmax ; endif
       open(unit=3,file=trim(tmpdirname)//'departs_inverses')
       open(unit=4,file=trim(tmpdirname)//'trous_grille_curvi')
 1557  read(3,*,end=1558)i2,j2,deci
       call curvgrdtoolbox_isoline_i(jbase_,deci)

! Reste � interpoler les points � l'interieur de l'iles qui n'ont pas
! �t� atteints. Ces derniers sont detect�s par differences d'extremites indice j
       read(4,*)i,j1
       write(6,*)'Trajectoire inverse. Manque j ',i,i2,j1,j2
       if(jbase_==jmax) then
        if(j2+1>j1-1) then !>>>
         do j=j1-1,j2+1
          rap=real(j-j1+1)/real(j2-j1+2)
!         rap=0.5*(1.+sin((rap-0.5)*pi))
          pos_x(i,j)=(1.-rap)*pos_x(i,j1-1)+rap*pos_x(i,j2+1)
          pos_y(i,j)=(1.-rap)*pos_y(i,j1-1)+rap*pos_y(i,j2+1)
!         if(i==590)write(6,*)i,j,real(rap),real(pos_x(i,j)),real(pos_y(i,j))
         enddo
        endif              !>>>
       else
        if(j1+1>j2-1) then !>>>
        do j=j2-1,j1+1     ! attention cas non teste
          rap=real(j-j2+1)/real(j1-j2+2)
          pos_x(i,j)=(1.-rap)*pos_x(i,j2-1)+rap*pos_x(i,j1+1)
          pos_y(i,j)=(1.-rap)*pos_y(i,j2-1)+rap*pos_y(i,j1+1)
        enddo
        endif              !>>>
       endif

       goto 1557
 1558  close(3)
       close(4)
       if(jbase_==jmax) then ; jbase_=1 ; else ; jbase_=jmax ; endif


      endif             !//////>
#endif

! Calculer les lon lat correspondant:
      do j=1,jmax
      do i=1,imax

       i1=int(pos_x(i,j))
       j1=int(pos_y(i,j))
       rapi=pos_x(i,j)-i1
       rapj=pos_y(i,j)-j1

       if(i1<0)     then
          write(6,*)'Pb ',i,j,pos_x(i,j),pos_y(i,j)
          stop 'i1<0 pos_lon pos_lat'
       endif
       if(j1<0)     stop 'j1<0 pos_lon pos_lat'
       if(i1>imax+1)stop 'i1>imax+1 pos_lon pos_lat'
       if(j1>jmax+1)stop 'j1>jmax+1 pos_lon pos_lat'

       x1=cos(lat_t(i1  ,j1  ))*cos(lon_t(i1  ,j1  ))
       y1=cos(lat_t(i1  ,j1  ))*sin(lon_t(i1  ,j1  ))
       z1=sin(lat_t(i1  ,j1  ))

       x2=cos(lat_t(i1+1,j1  ))*cos(lon_t(i1+1,j1  ))
       y2=cos(lat_t(i1+1,j1  ))*sin(lon_t(i1+1,j1  ))
       z2=sin(lat_t(i1+1,j1  ))

       x3=cos(lat_t(i1  ,j1+1))*cos(lon_t(i1  ,j1+1))
       y3=cos(lat_t(i1  ,j1+1))*sin(lon_t(i1  ,j1+1))
       z3=sin(lat_t(i1  ,j1+1))

       x4=cos(lat_t(i1+1,j1+1))*cos(lon_t(i1+1,j1+1))
       y4=cos(lat_t(i1+1,j1+1))*sin(lon_t(i1+1,j1+1))
       z4=sin(lat_t(i1+1,j1+1))

       x5=(1.-rapi)*(1.-rapj)*x1     &
         +    rapi *(1.-rapj)*x2     &
         +(1.-rapi)*    rapj *x3     &
         +    rapi *    rapj *x4

       y5=(1.-rapi)*(1.-rapj)*y1     &
         +    rapi *(1.-rapj)*y2     &
         +(1.-rapi)*    rapj *y3     &
         +    rapi *    rapj *y4

       z5=(1.-rapi)*(1.-rapj)*z1     &
         +    rapi *(1.-rapj)*z2     &
         +(1.-rapi)*    rapj *z3     &
         +    rapi *    rapj *z4

       x0=sqrt(x5**2+y5**2+z5**2) ! rayon "terre" diminuE par l'interpolation
       pos_lat(i,j)=asin(z5/x0)
       pos_lon(i,j)=atan2(y5,x5)

!      pos_lat(i,j)=(1.-rapi)*(1.-rapj)*lat_t(i1  ,j1  ) &
!                  +    rapi *(1.-rapj)*lat_t(i1+1,j1  ) &
!                  +(1.-rapi)*    rapj *lat_t(i1  ,j1+1) &
!                  +    rapi *    rapj *lat_t(i1+1,j1+1)
!      pos_lon(i,j)=(1.-rapi)*(1.-rapj)*lon_t(i1  ,j1  ) &
!                  +    rapi *(1.-rapj)*lon_t(i1+1,j1  ) &
!                  +(1.-rapi)*    rapj *lon_t(i1  ,j1+1) &
!                  +    rapi *    rapj *lon_t(i1+1,j1+1)

! First guess du masque terre/mer
       xy_t(i,j,1)=(1.-rapi)*(1.-rapj)*mask_t(i1  ,j1  ,kmaxp1) &
                  +    rapi *(1.-rapj)*mask_t(i1+1,j1  ,kmaxp1) &
                  +(1.-rapi)*    rapj *mask_t(i1  ,j1+1,kmaxp1) &
                  +    rapi *    rapj *mask_t(i1+1,j1+1,kmaxp1)
       if(xy_t(i,j,1)<0.99)xy_t(i,j,1)=0.

! First guess de bathy
       xy_t(i,j,3)=(1.-rapi)*(1.-rapj)*h_w(i1  ,j1  ) &
                  +    rapi *(1.-rapj)*h_w(i1+1,j1  ) &
                  +(1.-rapi)*    rapj *h_w(i1  ,j1+1) &
                  +    rapi *    rapj *h_w(i1+1,j1+1)
      enddo
      enddo

! VERIFIER QU'IL N'Y AIT PAS DE POINTS COLLES:
      do j=1,jmax ; do i=2,imax
       if(abs(pos_lat(i,j)-pos_lat(i-1,j))+abs(pos_lon(i,j)-pos_lon(i-1,j))<1.d-10) then
        write(6,*)'Points collEs en i,j,i-1,j',i,j,i-1,j
        write(6,*)'2 pos_x',pos_x(i-1:i,j)
        write(6,*)'2 pos_y',pos_y(i-1:i,j)
! Sortir la trajectoire buguee dans un fichier trajectoire.bug
         open(unit=3,file='trajectoire1.bug')
          do j1=1,jmax
           write(3,*)j,pos_x(i-1,j1),pos_y(i-1,j1)
          enddo
         close(3)
         open(unit=3,file='trajectoire2.bug')
          do j1=1,jmax
           write(3,*)j,pos_x(i,j1),pos_y(i,j1)
          enddo
         close(3)
        write(6,*)'Voir trajectoire1.bug trajectoire2.bug'
        stop 'Err 575'
       endif
      enddo ; enddo
      do j=2,jmax ; do i=1,imax
       if(abs(pos_lat(i,j)-pos_lat(i,j-1))+abs(pos_lon(i,j)-pos_lon(i,j-1))<1.d-10) then
        write(6,*)'Points collEs en ',i,j,i,j-1
        stop 'Err 576'
       endif
      enddo ; enddo

! Calcul des angles des axes de la nouvelle grille:
      do j=2,jmax-1
      do i=2,imax-1

! Angle1
       dlon_di_=pos_lon(i+1,j)-pos_lon(i-1,j)
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
       xy_t(i,j,2)=atan2(pos_lat(i+1,j)-pos_lat(i-1,j)  &
                        ,dlon_di_*cos(pos_lat(i,j)))
!      xy_t(i,j,2)=atan2(pos_y(i+1,j)-pos_y(i-1,j)  &
!                       ,pos_x(i+1,j)-pos_x(i-1,j))

! Angle2-angle1
       dlon_dj_=pos_lon(i,j+1)-pos_lon(i,j-1)
       if(dlon_dj_<-pi)dlon_dj_=dlon_dj_+2.*pi
       if(dlon_dj_> pi)dlon_dj_=dlon_dj_-2.*pi
       xy_t(i,j,2)=-xy_t(i,j,2)                          &
                  +atan2(pos_lat(i,j+1)-pos_lat(i,j-1)   &
                        ,dlon_dj_*cos(pos_lat(i,j)))
!      xy_t(i,j,2)=-xy_t(i,j,2)                          &
!                 +atan2(pos_y(i,j+1)-pos_y(i,j-1)   &
!                       ,pos_x(i,j+1)-pos_x(i,j-1))

! Conversion degres:
       xy_t(i,j,2)=xy_t(i,j,2)*rad2deg
       if(xy_t(i,j,2)<-180.)xy_t(i,j,2)=xy_t(i,j,2)+360.

! Appliquer le first guess du masque:
!      if(xy_t(i,j,1)<0.99)xy_t(i,j,2)=filval_

!      if(i==250)write(6,*)j,xy_t(i,j,2)

! Rappel des angles de la vieille grille
! Angle1
       dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
       xy_t(i,j,4)=atan2(lat_t(i+1,j)-lat_t(i-1,j)  &
                        ,dlon_di_*cos(lat_t(i,j)))

! Angle2-angle1
       dlon_dj_=lon_t(i,j+1)-lon_t(i,j-1)
       if(dlon_dj_<-pi)dlon_dj_=dlon_dj_+2.*pi
       if(dlon_dj_> pi)dlon_dj_=dlon_dj_-2.*pi
       xy_t(i,j,4)=-xy_t(i,j,4)                          &
                  +atan2(lat_t(i,j+1)-lat_t(i,j-1)   &
                        ,dlon_dj_*cos(lat_t(i,j)))

! Conversion degres:
       xy_t(i,j,4)=xy_t(i,j,4)*rad2deg
       if(xy_t(i,j,4)<-180.)xy_t(i,j,4)=xy_t(i,j,4)+360.

!      if(i==250) then
!        write(6,*)j,xy_t(i,j,4)
!        write(6,*)(lon_t(i+1,j)-lon_t(i-1,j)) &
!                 *(lat_t(i+1,j)-lat_t(i-1,j))/(dx_t(i,j)**2) &
!                 ,(lon_t(i,j+1)-lon_t(i,j-1)) &
!                 *(lat_t(i,j+1)-lat_t(i,j-1))/(dy_t(i,j)**2)
!      endif


      enddo
      enddo
! CL
      do j=2,jmax-1
       xy_t(imax,j,2)=xy_t(imax-1,j,2)
       xy_t(1   ,j,2)=xy_t(2     ,j,2)
      enddo
      do i=1,imax
       xy_t(i,jmax,2)=xy_t(i,jmax-1,2)
       xy_t(i,1   ,2)=xy_t(i,2     ,2)
      enddo

#ifdef bidon
! Dessiner la grille gnuplot dans plan (x,y)
      do j2=1,jmax
       do i2=1,imax
        write(10,*)i2,j2,rad2deg*pos_lon(i2,j2),rad2deg*pos_lat(i2,j2)
       enddo
       do i2=imax,1,-1
        write(10,*)i2,j2,rad2deg*pos_lon(i2,j2),rad2deg*pos_lat(i2,j2)
       enddo
      enddo

      do i2=1,imax
       do j2=1,jmax
        write(11,*)i2,j2,rad2deg*pos_lon(i2,j2),rad2deg*pos_lat(i2,j2)
       enddo
       do j2=jmax,1,-1
        write(11,*)i2,j2,rad2deg*pos_lon(i2,j2),rad2deg*pos_lat(i2,j2)
       enddo
      enddo

      do j=1,jmax
      do i=1,imax
       if(xy_t(i,j,1)<1.)write(12,*)rad2deg*pos_lon(i,j),rad2deg*pos_lat(i,j)
      enddo
      enddo
#endif


! Pour l'ecriture du fichier lonlat calcul de dx et dy:
      do j=1,jmax
      do i=1,imax-1

       x1=pos_lon(i+1,j)
       if(x1-pos_lon(i,j)<-pi)x1=x1+2.*pi
       if(x1-pos_lon(i,j)> pi)x1=x1-2.*pi

! anyvar3d(:,:,1) est dx_u
        anyvar3d(i+1,j,1)=rayonterre*sqrt(                             &
         ((x1-pos_lon(i,j))*cos(0.5*(pos_lat(i,j)+pos_lat(i+1,j))))**2 &
        +(  pos_lat(i+1,j)- pos_lat(i,j)                  )**2 )

      enddo
      enddo

      do j=1,jmax-1
      do i=1,imax

       x1=pos_lon(i,j+1)
       if(x1-pos_lon(i,j)<-pi)x1=x1+2.*pi
       if(x1-pos_lon(i,j)> pi)x1=x1-2.*pi

! anyvar3d(:,:,2) est dy_v
        anyvar3d(i,j+1,2)=rayonterre*sqrt(                         &
         ((x1-pos_lon(i,j))*cos(0.5*(pos_lat(i,j)+pos_lat(i,j+1))))**2 & !15-02-18
        +(  pos_lat(i,j+1)- pos_lat(i,j)                  )**2 )

      enddo
      enddo

! Ecrire fichier netcdf
      status=nf_create('lonlat.nc',nf_clobber,ncid_)
      status=nf_def_dim(ncid_,'ni_t',imax,dim_x_id_)
      status=nf_def_dim(ncid_,'nj_t',jmax,dim_y_id_)

      vardim(1)=dim_x_id_ ; vardim(2)=dim_y_id_

!     status=nf_def_var(ncid_,'vlx_u',nf_double,2,vardim(1:2),var_id)
!     status=nf_def_var(ncid_,'vly_v',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'longitude_t',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'latitude_t',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'i_newgrid',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'j_newgrid',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'i_index',nf_real,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'j_index',nf_real,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'h_w',nf_double,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'dx_t',nf_real,2,vardim(1:2),var_id)
      status=nf_def_var(ncid_,'dy_t',nf_real,2,vardim(1:2),var_id)

      status=nf_def_var(ncid_,'mask_t',nf_int,2,vardim(1:2),var_id)

      status=nf_def_var(ncid_,'axis_angle',nf_double,2,vardim(1:2),var_id)
      status=nf_put_att_double(ncid_,var_id,'_FillValue' &
                               ,nf_double,len_,filval_)

      status=nf_enddef(ncid_)

!     status=nf_inq_varid(ncid_,'vlx_u',var_id)
!     status=nf_put_var_double(ncid_,var_id,vlx_u(1:imax,1:jmax))
!     status=nf_inq_varid(ncid_,'vly_v',var_id)
!     status=nf_put_var_double(ncid_,var_id,vly_v(1:imax,1:jmax))

! conversion en degres:
      pos_lon(1:imax,1:jmax)=pos_lon(1:imax,1:jmax)*rad2deg
      status=nf_inq_varid(ncid_,'longitude_t',var_id)
      status=nf_put_var_double(ncid_,var_id,pos_lon(1:imax,1:jmax))

      pos_lat(1:imax,1:jmax)=pos_lat(1:imax,1:jmax)*rad2deg
      status=nf_inq_varid(ncid_,'latitude_t',var_id)
      status=nf_put_var_double(ncid_,var_id,pos_lat(1:imax,1:jmax))

      status=nf_inq_varid(ncid_,'i_newgrid',var_id)
      status=nf_put_var_double(ncid_,var_id,pos_x(1:imax,1:jmax))

      status=nf_inq_varid(ncid_,'j_newgrid',var_id)
      status=nf_put_var_double(ncid_,var_id,pos_y(1:imax,1:jmax))

      do j=1,jmax ; do i=1,imax ; anyvar2d(i,j)=i ; enddo ; enddo
      status=nf_inq_varid(ncid_,'i_index',var_id)
      status=nf_put_var_real(ncid_,var_id,anyvar2d(1:imax,1:jmax))

      do j=1,jmax ; do i=1,imax ; anyvar2d(i,j)=j ; enddo ; enddo
      status=nf_inq_varid(ncid_,'j_index',var_id)
      status=nf_put_var_real(ncid_,var_id,anyvar2d(1:imax,1:jmax))

      status=nf_inq_varid(ncid_,'h_w',var_id)
      status=nf_put_var_double(ncid_,var_id,xy_t(1:imax,1:jmax,3))

      do j=1,jmax ; do i=2,imax-1
       anyvar2d(i,j)=0.5*(anyvar3d(i,j,1)+anyvar3d(i+1,j,1))*xy_t(i,j,1)
      enddo       ; enddo
      anyvar2d(1   ,:)=0.
      anyvar2d(imax,:)=0.
      status=nf_inq_varid(ncid_,'dx_t',var_id)
      status=nf_put_var_real(ncid_,var_id,anyvar2d(1:imax,1:jmax))

      do j=2,jmax-1 ; do i=1,imax
       anyvar2d(i,j)=0.5*(anyvar3d(i,j,2)+anyvar3d(i,j+1,2))*xy_t(i,j,1)
      enddo       ; enddo
      status=nf_inq_varid(ncid_,'dy_t',var_id)
      status=nf_put_var_real(ncid_,var_id,anyvar2d(1:imax,1:jmax))


      status=nf_inq_varid(ncid_,'mask_t',var_id)
      status=nf_put_var_int(ncid_,var_id,nint(xy_t(1:imax,1:jmax,1)))

      status=nf_inq_varid(ncid_,'axis_angle',var_id)
      status=nf_put_var_double(ncid_,var_id,xy_t(1:imax,1:jmax,2))

      status=nf_close(ncid_)

!     open(unit=3,file=trim(tmpdirname)//'bathy_curvi.dat',recl=100000)
      open(unit=3,file='bathy_curvi.dat',recl=100000)
        do i=1,imax
        write(3,'(4000i1)')(  nint(xy_t(i,j,1)),j=1,jmax)
        enddo
        do i=1,imax
        write(3,'(4000(f10.3,1x))')(xy_t(i,j,3),j=1,jmax)
        enddo
      if(par%rank==0)write(6,*)'bathy_curvi.dat ok'
      close(3)

      open(unit=3,file='bathycote_in.ijhlonlat',recl=100000)
        do i=1,imax
        write(3,'(4000i1)')(  nint(xy_t(i,j,1)),j=1,jmax)
        enddo
        do j=1,jmax
        do i=1,imax
         write(3,'(2i6,1x,f10.4,1x,f15.10,1x,f15.10)') &
         i,j,xy_t(i,j,3),pos_lon(i,j),pos_lat(i,j)
        enddo
        enddo
      if(par%rank==0)write(6,*)'bathycote_in.ijhlonlat ok'
      close(3)

      open(unit=3,file='lonlat_4col.txt',recl=100000)
        do j=1,jmax
        do i=1,imax
         write(3,'(2i6,1x,f15.10,1x,f15.10)') &
         i,j,pos_lon(i,j),pos_lat(i,j)
        enddo
        enddo
      if(par%rank==0)write(6,*)'lonlat_4col.txt ok'
      close(3)

      open(unit=3,file='lonlat_bounds.txt')
       j=1
       do i=1,imax
       write(3,*)pos_lon(i,j),pos_lat(i,j)
       enddo
       i=imax
       do j=1,jmax
       write(3,*)pos_lon(i,j),pos_lat(i,j)
       enddo
       j=jmax
       do i=imax,1,-1
       write(3,*)pos_lon(i,j),pos_lat(i,j)
       enddo
       i=1
       do j=jmax,1,-1
       write(3,*)pos_lon(i,j),pos_lat(i,j)
       enddo
      close(3)

!     stop 'Chet'
      endif             !2222222222222>


      end subroutine curvgrdtoolbox_driver

!...................................................................

      subroutine curvgrdtoolbox_smooth_varj(istart_,istop_,jstart_,jstop_,loopmax_,txt_)
      implicit none
      double precision dt_
      integer jstart_,jstop_,istart_,istop_,loopmax_
      character*4 txt_
#ifdef synopsis
       subroutinetitle='curvgrdtoolbox_smooth_varj'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! DEfinir le coef de diffusion. Par dEfaut 0.5 (valeur max garantissant la stabilite du pas de temps)
      xy_t(:,:,5)=0.5

      dt_=0.25
      do loop1=1,loopmax_
      if(mod(loop1,100)==0.and.par%rank==0)write(6,*)'case1 ',loop1

      if(txt_=='mask') then !mmmmmmmm>
! Lisser la partie marine seulement
        do i=istart_,istop_
        im1=max(i-1,0) ; ip1=min(i+1,imax+1)
        do j=jstart_,jstop_

!        var_j(i,j)=                                       &
         xy_t(i,j,6)=                                      & !04-06-16
         var_j(i,j)+dt_*(                                  &

                 (var_j(i  ,j+1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,j+1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(i  ,j-1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,j-1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(ip1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(ip1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(im1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(im1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                  )*mask_curvi(i  ,j  )

        enddo ; enddo

      endif                 !mmmmmmmm>

      if(txt_=='invm') then !mmmmmmmm> !05-06-16
! Lisser la partie continentale seulement
        do i=istart_,istop_
        im1=max(i-1,0) ; ip1=min(i+1,imax+1)
        do j=jstart_,jstop_

!        var_j(i,j)=                                       &
         xy_t(i,j,6)=                                      & !04-06-16
         var_j(i,j)+dt_*(                                  &

                 (var_j(i  ,j+1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,j+1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(i  ,j-1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,j-1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(ip1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(ip1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(im1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(im1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                  )*(1-mask_curvi(i,j))

        enddo ; enddo

      endif                 !mmmmmmmm>

      if(txt_=='free') then !ffffffff>
! Lisser tout

        do i=istart_,istop_
        im1=max(i-1,0) ; ip1=min(i+1,imax+1)
        do j=jstart_,jstop_
        jm1=max(j-1,0) ; jp1=min(j+1,jmax+1)

!        var_j(i,j)=                                       &
         xy_t(i,j,6)=                                      & !04-06-16
         var_j(i,j)+dt_*(                                  &

                 (var_j(i  ,jp1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,jp1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(i  ,jm1) -var_j(i  ,j  ))          &
                *( xy_t(i  ,jm1,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(ip1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(ip1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                +(var_j(im1,j  ) -var_j(i  ,j  ))          &
                *( xy_t(im1,j  ,5)+xy_t(i  ,j  ,5))        & ! coef dif

                  )

        enddo ; enddo

      endif                 !ffffffff>

      do i=istart_,istop_ ; do j=jstart_,jstop_
       var_j(i,j)=xy_t(i,j,6) !04-06-16
      enddo ; enddo

#ifdef parallele
      call get_type_echange('za','var_j_za',var_j,lbound(var_j),ubound(var_j),i1)
      do loop2=1,subcycle_exchange
        call echange_voisin(var_j,i1,mpi_neighbor_list(loop2)) !31-07-14
      enddo
      call loc_wait()
#endif

      enddo ! loop1

      end subroutine curvgrdtoolbox_smooth_varj

!...................................................................

      subroutine curvgrdtoolbox_smooth_vari(istart_,istop_,jstart_,jstop_,loopmax_,txt_)
      implicit none
      double precision dt_
      integer jstart_,jstop_,istart_,istop_,loopmax_
      character*4 txt_
#ifdef synopsis
       subroutinetitle='curvgrdtoolbox_smooth_vari'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! DEfinir le coef de diffusion. Par dEfaut 0.5 (valeur max garantissant la stabilite du pas de temps)
      xy_t(:,:,5)=0.5

      dt_=0.25
      do loop1=1,loopmax_
      if(mod(loop1,10)==0)write(6,*)'case1 ',loop1

      if(txt_=='mask') then !mmmmmmmm>
        do j=jstart_,jstop_
        jm1=max(j-1,0) ; jp1=min(j+1,jmax+1)
        do i=istart_,istop_

!         var_i(i,j)=                                    &
          xy_t(i,j,6)=                                   & !04-06-16
          var_i(i,j)+dt_*(                               &

                      (var_i(i  ,jp1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jp1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i  ,jm1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jm1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i+1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(i+1,j  ,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i-1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(i-1,j  ,5)+xy_t(i  ,j  ,5)) &

                         )*mask_curvi(i,j)

        enddo ; enddo
      endif                 !mmmmmmmm>

      if(txt_=='invm') then !imimimimim>>
        do j=jstart_,jstop_
        jm1=max(j-1,0) ; jp1=min(j+1,jmax+1)
        do i=istart_,istop_

!         var_i(i,j)=                                    &
          xy_t(i,j,6)=                                   & !04-06-16
          var_i(i,j)+dt_*(                               &

                      (var_i(i  ,jp1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jp1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i  ,jm1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jm1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i+1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(i+1,j  ,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i-1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(i-1,j  ,5)+xy_t(i  ,j  ,5)) &

                         )*(1-mask_curvi(i,j))

        enddo ; enddo
      endif                 !imimimimim>>

      if(txt_=='free') then !ffffffff>
        do j=jstart_,jstop_
        jm1=max(j-1,0) ; jp1=min(j+1,jmax+1)
        do i=istart_,istop_
        im1=max(i-1,0) ; ip1=min(i+1,imax+1)

!         var_i(i,j)=                                    &
          xy_t(i,j,6)=                                   & !04-06-16
          var_i(i,j)+dt_*(                               &

                      (var_i(i  ,jp1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jp1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(i  ,jm1) -var_i(i  ,j  ))   &
                     *( xy_t(i  ,jm1,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(ip1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(ip1,j  ,5)+xy_t(i  ,j  ,5)) &

                     +(var_i(im1,j  ) -var_i(i  ,j  ))   &
                     *( xy_t(im1,j  ,5)+xy_t(i  ,j  ,5)) &

                         )

        enddo ; enddo
      endif                 !ffffffff>

      do j=jstart_,jstop_ ; do i=istart_,istop_
       var_i(i,j)=xy_t(i,j,6) !04-06-16
      enddo ; enddo

#ifdef parallele
      call get_type_echange('za','var_i_za',var_i,lbound(var_i),ubound(var_i),i1)
      do loop2=1,subcycle_exchange
        call echange_voisin(var_i,i1,mpi_neighbor_list(loop2)) !31-07-14
      enddo
      call loc_wait()
#endif

      enddo ! loop1

      end subroutine curvgrdtoolbox_smooth_vari

!...................................................................

!...................................................................

      subroutine curvgrdtoolbox_isoline_i(jbase_,deci_b)
      implicit none
      real signevel_
      integer jbase_,counter_
      double precision dt_,deci_b,decj_b,indice_j_b &
                          ,deci_a,decj_a,indice_j_a &
                          ,vlu_,vlv_,filval_,dlon_di_ &
                          ,dlon_dj_,dist_,lat1_,lat2_,lon2_,lon1_ &
                          ,deci_bb,decj_bb
#ifdef synopsis
       subroutinetitle='curvgrdtoolbox_isoline_i'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      counter_=0
      if(jbase_==jmax)signevel_=-1.  ! Chine:  trajet de jmax vers 1
      if(jbase_==1   )signevel_=+1.  ! Taiwan: trajet de 1 vers jmax

!     deci_b=i2

      decj_b=real(jbase_) ! jbase_=jmax pour Chine:  trajet de jmax vers 1
                          ! jbase_=1    pour Taiwan: trajet de 1 vers jmax
      indice_j_b=decj_b

      pos_x(i2,nint(decj_b))=deci_b
      pos_y(i2,nint(decj_b))=decj_b

      if(jbase_==jmax)j2=jmax-1  ! Chine:  trajet de jmax vers 1
      if(jbase_==1   )j2=2       ! Taiwan: trajet de 1 vers jmax

!      do loop1=1,1
!      do loop1=1,10000
      do loop1=1,100000 ! 01-02-21 Cette boucle doit etre grande

       i1=int(deci_b+0.5)
       j1=int(decj_b)
       rapi=deci_b+0.5-real(i1)
       rapj=decj_b    -real(j1)

       if(i1<1)        then
          write(6,*)'deci_b+0.5 i1=',deci_b+0.5,i1
          stop ' Stop isolinge_i i1<1 vlx_u'
       endif
       if(i1+1>imax+1)stop ' Stop i1+1>imax+1 vlx_u'
       if(j1<0)       stop ' Stop j1<0 vlx_u'
       if(j1+1>jmax+1)stop ' Stop j1+1>jmax+1 vlx_u'

       vlu_=((1.-rapi)*(1.-rapj)*vlx_u(i1  ,j1  ) &
            +(1.-rapi)*    rapj *vlx_u(i1  ,j1+1) &
            +    rapi *(1.-rapj)*vlx_u(i1+1,j1  ) &
            +    rapi *    rapj *vlx_u(i1+1,j1+1))*signevel_

!      if(loop1==2070) then
!      if(jbase_==jmax.and.i2==1) then
!         write(6,*)'vlu',vlu_,i1,j1
!         write(6,*)vlx_u(i1  ,j1  ),vlx_u(i1  ,j1+1),vlx_u(i1+1,j1  ),vlx_u(i1+1,j1+1)
!         stop 'bibi'
!      endif
!      endif

       i1=int(deci_b)
       j1=int(decj_b+0.5)
       rapi=deci_b    -real(i1)
       rapj=decj_b+0.5-real(j1)

       if(i1<0)       stop ' Stop i1<0 vly_v'
       if(j1<1)       stop ' Stop j1<1 vly_v'

       vlv_=((1.-rapi)*(1.-rapj)*vly_v(i1  ,j1  ) &
            +(1.-rapi)*    rapj *vly_v(i1  ,j1+1) &
            +    rapi *(1.-rapj)*vly_v(i1+1,j1  ) &
            +    rapi *    rapj *vly_v(i1+1,j1+1))*signevel_


!      if(jbase_==jmax.and.i2==2)write(666,*)deci_b,decj_b,vlv_

!      write(6,*)'i1,j1,vlv_=',i1,j1,vlv_
!      write(6,*)vly_v(i1  ,j1  ),vly_v(i1  ,j1+1)
!      write(6,*)vly_v(i1+1,j1  ),vly_v(i1+1,j1+1)

!      x1=sqrt(vlu_**2+vlv_**2)
       x1=max( abs(vlu_) , abs(vlv_) )
       if(x1==0.) then
        write(6,*)'vlu vlv=0 donc goto 777'
        goto 777
       endif
!      dt_=0.25/x1
       dt_=0.05/x1
!      dt_=0.01/x1


       deci_bb=deci_b
       decj_bb=decj_b
 2026 continue

       counter_=counter_+1

!     if(i2==1)write(67,*)counter_
!     if(i2==1)write(67,*)vlu_*dt_,vlv_*dt_
!     if(i2==imax)write(67,*)deci_b,decj_b
!     write(10+i2,*)deci_b,decj_b
!     if(i2==3)write(67,*)'deci_b,decj_b',deci_b,decj_b

! Nouvelle position dans le vieux repere (i,j)
       deci_a=deci_b+vlu_*dt_
       decj_a=decj_b+vlv_*dt_

! Bornages et periodicite:
      decj_a=min(max(decj_a,un),un*jmax)

      if(iperiodicboundary) then !------>
       if(deci_a>imax) then !aaaaa>
          deci_a=deci_a-imax+2
          deci_b=deci_b-imax+2 !26-04-21
!         if(deci_b>imax)deci_b=deci_b-imax+2
       endif                !aaaaa>
       if(deci_a<1   ) then !bbbbb>
          deci_a=deci_a+imax-2
          deci_b=deci_b+imax-2 !26-04-21
!         if(deci_b<1   )deci_b=deci_b+imax-2
       endif                !bbbbb>
      else                       !------>

!      deci_a=min(max(deci_a,un),un*imax)
! Nouveau bornage: pour empecher que le bornage n'entraine des points collEs, la borne depend de i2 !01-02-21
       deci_a=min(max(deci_a,1.+i2*0.01) &          !i2*0.01 pour ne pas avoir plusieurs isoligne en 1
                            ,imax-(imax-i2)*0.01)   !(imax-i2)*0.01 pour ne pas avoir plusieurs isoligne en imax


      endif                      !------>

!     if(i2==1)write(67,*)deci_a,decj_a

       i1=int(deci_a)
       j1=int(decj_a)
       rapi=deci_a-real(i1)
       rapj=decj_a-real(j1)

       if(i1<0)       then
         write(6,*)'vlu_ vlv_=',vlu_,vlv_
         write(6,*)'deci_b decj_b=',deci_b,decj_b
         write(6,*)'deci_a decj_a=',deci_a,decj_a
         stop ' Stop i1<0 var_j'
       endif
       if(i1+1>imax+1)stop ' Stop i1+1>imax+1 var_j'
       if(j1<0)       then
         write(6,*)'vlu_ vlv_=',vlu_,vlv_
         write(6,*)'deci_b decj_b=',deci_b,decj_b
         write(6,*)'deci_a decj_a=',deci_a,decj_a
         stop ' Stop j1<0 var_j'
       endif
       if(j1+1>jmax+1)stop ' Stop j1+1>jmax+1 var_j'

! Valeur local de j dans le nouveau repere
       indice_j_a=                            &
         (1.-rapi)*(1.-rapj)*var_j(i1  ,j1  ) &
        +(1.-rapi)*    rapj *var_j(i1  ,j1+1) &
        +    rapi *(1.-rapj)*var_j(i1+1,j1  ) &
        +    rapi *    rapj *var_j(i1+1,j1+1)

       if(jbase_==jmax) then !>>>
         if(var_j(i1,1)>1.or.var_j(i1+1,1)>1) then !ooo>
          write(6,*)'PROBLEME: on dirait que la plus petite valeur' &
      ,' de var_j est >1 ce qui empechera le nouvel axe d atteindre 1' &
         ,' depuis jmax. Revoir calcul de var_i en step 1'
          write(6,*)'var_j(',i1,',1)  ',var_j(i1,1)
          write(6,*)'var_j(',i1+1,',1)',var_j(i1+1,1)
          stop 'var_j(:,1)>1'
         endif                                         !ooo>
       endif              !>>>

!      if(i2==1)write(67,*)'indice_j_a',indice_j_a
!      if(i2==1)write(67,*)'i1 j1',i1,j1
!      if(i2==1)write(67,*)var_j(i1  ,j1  ),var_j(i1  ,j1+1),var_j(i1+1,j1  ),var_j(i1+1,j1+1)

! Si pas de temps trop grand:
      if(abs(indice_j_a-indice_j_b)>0.5) then
!      if(i2==1)write(67,*)'indice_j_a,indice_j_b',indice_j_a,indice_j_b
!      if(counter_==1)stop 'koko'
       deci_b=deci_bb ; decj_b=decj_bb ; dt_=dt_*0.5 ; goto 2026
      endif

!      if(jbase_==jmax.and.i2==2)write(666,*)deci_b,decj_b,vlv_

! Si on vient juste de depasser la valeur entiere de j2 alors jalonner:
       if( (indice_j_a-j2)*(indice_j_b-j2)<=0) then !------>

        if(j2>jmax)stop 'j2>jmax pos_x pos_y'
        if(j2<1   )stop 'j2<1    pos_x pos_y'
        if(i2>imax)stop 'i2>imax pos_x pos_y'
        if(i2<1   )stop 'i2<1    pos_x pos_y'

! Pos_x Pos_y position du point (i2,j2) dans le vieux repere
        rapj=(indice_j_a-real(j2))    &
            /(indice_j_a-indice_j_b)
        pos_x(i2,j2)=(1.-rapj)*deci_a+rapj*deci_b
        pos_y(i2,j2)=(1.-rapj)*decj_a+rapj*decj_b

!      if(jbase_==jmax.and.i2==2)write(667,*)i2,j2,pos_x(i2,j2),pos_y(i2,j2)
!      if(jbase_==jmax.and.i2==2)write(668,*)real(deci_a),real(decj_a),real(deci_b),real(decj_b)

        if(jbase_==jmax)j2=j2-1 ! Chine:  trajet de jmax vers 1
        if(jbase_==1   )j2=j2+1 ! Taiwan: trajet de 1 vers jmax
        if(j2>jmax)goto 777
        if(j2<1   )goto 777


       endif                                     !------>

       indice_j_b=indice_j_a
       deci_b=deci_a
       decj_b=decj_a

!      if(mask(nint(deci_a),nint(decj_a))==0)goto 777

!     if(i2==imax)write(67,*)'j2=',j2
      enddo !loop1
  777 continue
      if(jbase_==jmax.and.j2>1   )stop ' stop car sortie 777 et j2>1'    ! Chine:  trajet jmax vers 1
      if(jbase_==1   .and.j2<jmax) then
!      write(6,*)'loop1=',loop1
!      write(6,*)'jbase_=',jbase_
!      write(6,*)'j2=',j2
!      write(6,*)'jmax=',jmax
!      write(6,*)'indice_j_a=',indice_j_a
!      write(6,*)'indice_j_b=',indice_j_b
!      stop ' stop car sortie 777 et j2<jmax'
! La nouvelle grille va moins loin que jmax et donc on comble simplement la partie manquante de pox,posy
       pos_x(i2,j2:jmax)=pos_x(i2,j2)
       pos_y(i2,j2:jmax)=pos_y(i2,j2)
      endif
  888 continue

      end subroutine curvgrdtoolbox_isoline_i

!...................................................................

      subroutine curvgrdtoolbox_isoline_j(ibase_,decj_b)
      implicit none
      real signevel_
      integer ibase_
      double precision dt_,deci_b,decj_b,indice_i_b &
                          ,deci_a,decj_a,indice_i_a &
                          ,vlu_,vlv_,filval_,dlon_di_ &
                          ,dlon_dj_,dist_,lat1_,lat2_,lon2_,lon1_ &
                          ,deci_bb,decj_bb
#ifdef synopsis
       subroutinetitle='curvgrdtoolbox_isoline_j'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ibase_==imax)signevel_=-1.  ! Chine:  trajet de jmax vers 1
      if(ibase_==1   )signevel_=+1.  ! Taiwan: trajet de 1 vers jmax



      deci_b=real(ibase_)

      indice_i_b=deci_b

      pos_x(nint(deci_b),j2)=deci_b
      pos_y(nint(deci_b),j2)=decj_b

      if(ibase_==imax)i2=imax-1  ! Chine:  trajet de jmax vers 1
      if(ibase_==1   )i2=2       ! Taiwan: trajet de 1 vers jmax

!     do loop1=1,10000
      do loop1=1,100000 ! CETTE BOUCLE DOIT ETRE GRANDE!!!! 01-02-21

       i1=int(deci_b+0.5)
       j1=int(decj_b)
       rapi=deci_b+0.5-real(i1)
       rapj=decj_b    -real(j1)

       if(i1<1)       stop ' Stop i1<1 vlx_u'
       if(i1+1>imax+1)stop ' Stop i1+1>imax+1 vlx_u'
       if(j1<0)       stop ' Stop j1<0 vlx_u'
       if(j1+1>jmax+1)stop ' Stop j1+1>jmax+1 vlx_u'

       vlu_=((1.-rapi)*(1.-rapj)*vlx_u(i1  ,j1  ) &
            +(1.-rapi)*    rapj *vlx_u(i1  ,j1+1) &
            +    rapi *(1.-rapj)*vlx_u(i1+1,j1  ) &
            +    rapi *    rapj *vlx_u(i1+1,j1+1))*signevel_

       i1=int(deci_b)
       j1=int(decj_b+0.5)
       rapi=deci_b    -real(i1)
       rapj=decj_b+0.5-real(j1)

       if(i1<0)       stop ' Stop i1<0 vly_v'
       if(j1<1)       stop ' Stop j1<1 vly_v'

       vlv_=((1.-rapi)*(1.-rapj)*vly_v(i1  ,j1  ) &
            +(1.-rapi)*    rapj *vly_v(i1  ,j1+1) &
            +    rapi *(1.-rapj)*vly_v(i1+1,j1  ) &
            +    rapi *    rapj *vly_v(i1+1,j1+1))*signevel_


!      x1=sqrt(vlu_**2+vlv_**2)
       x1=max( abs(vlu_) , abs(vlv_) )
       if(x1==0.) then
        write(6,*)'ERREUR vlu vlv=0 donc goto 777'
        write(6,*)'i1 j1',i1,j1
        write(6,*)'4vly=',vly_v(i1,j1),vly_v(i1,j1+1),vly_v(i1+1,j1),vly_v(i1+1,j1+1)
        goto 777
       endif
!      dt_=0.25/x1
       dt_=0.05/x1

       deci_bb=deci_b
       decj_bb=decj_b
 2026 continue

! ON AVANCE:
       deci_a=deci_b+vlu_*dt_
       decj_a=decj_b+vlv_*dt_

! Bornages et periodicite:
      deci_a=min(max(deci_a,un),un*imax)
      if(jperiodicboundary) then !------>
       if(decj_a>jmax) then !aaaaa>
          decj_a=decj_a-jmax+2
          decj_b=decj_b-jmax+2 !26-04-21
!         if(decj_b>jmax)decj_b=decj_b-jmax+2
       endif                !aaaaa>
       if(decj_a<1   ) then !bbbbb>
          decj_a=decj_a+jmax-2
          decj_b=decj_b+jmax-2 !26-04-21
!         if(decj_b<1   )decj_b=decj_b+jmax-2
       endif                !bbbbb>
      else                       !------>
!      decj_a=min(max(decj_a,un),un*jmax)
! Nouveau bornage: pour empecher que le bornage n'entraine des points collEs, la borne depend de j2 !01-02-21
       decj_a=min(max(decj_a,1.+j2*0.01) &          !j2*0.01 pour ne pas avoir plusieurs isoligne en 1
                            ,jmax-(jmax-j2)*0.01)   !(jmax-j2)*0.01 pour ne pas avoir plusieurs isoligne en jmax
      endif                      !------>


       i1=int(deci_a)
       j1=int(decj_a)
       rapi=deci_a-real(i1)
       rapj=decj_a-real(j1)

       if(i1<0)       stop ' Stop i1<0 var_j'
       if(i1+1>imax+1)stop ' Stop i1+1>imax+1 var_j'
       if(j1<0)       stop ' Stop j1<0 var_j'
       if(j1+1>jmax+1)stop ' Stop j1+1>jmax+1 var_j'

       indice_i_a=                            &
         (1.-rapi)*(1.-rapj)*var_i(i1  ,j1  ) &
        +(1.-rapi)*    rapj *var_i(i1  ,j1+1) &
        +    rapi *(1.-rapj)*var_i(i1+1,j1  ) &
        +    rapi *    rapj *var_i(i1+1,j1+1)

       if(ibase_==1) then !>>>
         if(var_i(1,j1  )>1.or.var_i(1,j1+1)>1) then !ooo>
          write(6,*)'PROBLEME: on dirait que la plus petite valeur' &
         ,' de var_i est > 1 ce qui empechera le nouvel axe de demarrer' &
         ,' depuis 1. Revoir calcul de var_i en step 1'
          write(6,*)'var_i(1,',j1  ,')=',var_i(1,j1  )
          write(6,*)'var_i(1,',j1+1,')=',var_i(1,j1+1)
         
          stop 'var_i(1,:)/=1'
         endif                                         !ooo>
       endif              !>>>
       if(ibase_==imax) then !>>>
         if(var_i(imax,j1  )<imax.or.var_i(imax,j1+1)<imax) then !ooo>
          write(6,*)'PROBLEME: on dirait que la plus grande valeur' &
      ,' de var_i est <imax ce qui empechera le nouvel axe de demarrer' &
         ,' depuis imax. Revoir calcul de var_i en step 1'
          write(6,*)'var_i(imax,',j1  ,')=',var_i(imax,j1  )
          write(6,*)'var_i(imax,',j1+1,')=',var_i(imax,j1+1)
         
          stop 'var_i(imax,:)/=imax'
         endif                                         !ooo>
       endif              !>>>

! Si pas de temps trop grand:
      if(abs(indice_i_a-indice_i_b)>0.5) then
       deci_b=deci_bb ; decj_b=decj_bb ; dt_=dt_*0.5 ; goto 2026
      endif

!      if( mask_curvi(i1  ,j1  )    &
!         *mask_curvi(i1+1,j1  )    &
!         *mask_curvi(i1  ,j1+1)    &
!         *mask_curvi(i1+1,j1+1)==0) goto 888

!       if(j2==1)write(666,*)i2,indice_i_b,indice_i_a

       if( (indice_i_a-i2)*(indice_i_b-i2)<=0) then !------>

        if(j2>jmax)stop 'j2>jmax pos_x pos_y'
        if(j2<1   )stop 'j2<1    pos_x pos_y'
        if(i2>imax)stop 'i2>imax pos_x pos_y'
        if(i2<1   )stop 'i2<1    pos_x pos_y'

        rapi=(indice_i_a-real(i2))    &
            /(indice_i_a-indice_i_b)
        pos_x(i2,j2)=(1.-rapi)*deci_a+rapi*deci_b
        pos_y(i2,j2)=(1.-rapi)*decj_a+rapi*decj_b

!       if(j2==1)write(666,*)i2,pos_x(i2,j2),pos_y(i2,j2)


        if(ibase_==imax)i2=i2-1 ! Chine:  trajet de jmax vers 1
        if(ibase_==1   )i2=i2+1 ! Taiwan: trajet de 1 vers jmax
        if(i2>imax)goto 777
        if(i2<1   )goto 777

       endif                                     !------>

       indice_i_b=indice_i_a
       deci_b=deci_a
       decj_b=decj_a

      enddo ! loop1?
  777 continue
      if(ibase_==imax.and.i2>1   ) then
        write(6,*)'loop1=',loop1
        write(6,*)'ibase_=',ibase_
        write(6,*)'i2 j2=',i2,j2
        stop ' stop car sortie 777 et i2>1'
      endif
      if(ibase_==1   .and.i2<imax) then
!       write(6,*)'loop1=',loop1
!       write(6,*)'ibase_=',ibase_
!       write(6,*)'i2 j2=',i2,j2
!       write(6,*)'imax=',imax
!       stop ' stop car sortie 777 et i2<imax'
! La nouvelle grille va moins loin que imax et donc on comble simplement la partie manquante de pox,posy
!      pos_x(i2,j2:jmax)=pos_x(i2,j2)
!      pos_y(i2,j2:jmax)=pos_y(i2,j2)
       pos_x(i2:imax,j2)=pos_x(i2,j2) !01-02-21
       pos_y(i2:imax,j2)=pos_y(i2,j2)
      endif
  888 continue

      end subroutine curvgrdtoolbox_isoline_j

!...............................................................
      subroutine curvgrdtoolbox_obc_varj(delta_i1_,delta_i2_)
      implicit none
      integer delta_i1_,delta_i2_
! Appliquer un gradient nul sur les bords de var_j pour garantir que
! les trajectoires "iso-i" ne sortent pas du domaine
       do j=0,jmax+1
! Bord i=1
        do i=0,delta_i1_
         var_j(i,j)=var_j(delta_i1_+1,j)
        enddo
! Bord i=imax
        do i=imax+1,imax+1-delta_i2_,-1
         var_j(i,j)=var_j(imax+1-delta_i2_-1,j)
        enddo
       enddo
      end subroutine curvgrdtoolbox_obc_varj
!...............................................................
      subroutine curvgrdtoolbox_obc_vari(delta_j1_,delta_j2_)
      implicit none
      integer delta_j1_,delta_j2_
! Appliquer un gradient nul sur les bords de var_i pour garantir que
! les trajectoires "iso-j" ne sortent pas du domaine
       do i=0,imax+1
! Bord j=1
        do j=0,delta_j1_
         var_j(i,j)=var_j(i,delta_j1_+1)
        enddo
! Bord j=jmax
        do j=jmax+1,jmax+1-delta_j2_,-1
         var_j(i,j)=var_j(i,jmax+1-delta_j2_-1)
        enddo
       enddo
      end subroutine curvgrdtoolbox_obc_vari
!...............................................................

      end module module_curvgrdtoolbox
