      subroutine dragcoef
!______________________________________________________________________
! SYMPHONIE ocean model
! release 368 - last update: 14-04-23
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_q
      implicit none
! Computes the bottom drag coeffient
#ifdef synopsis
       subroutinetitle='dragcoef'
       subroutinedescription=' Drag coefficient driver subroutine'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         10/03/02: ponderation par la proportion de point de mer entourant
!                   le point CDB_Z
!         04/03/03: passage à la coordonnee hybride.
!         17/01/06: Cas de la prise en compte de la houle
!         26/01/06: finalisation du point precedent
! 2009.3  01-10-09: Suppression de la ponderation qui avait été ajoutée le 10/03/02
!         13-10-09: Attention le cas 2D est maintenant associé à la valeur de nr
!                   Remplacement de 1.+small par 1.1d0
! 2010.8  03-05-10  modif sur test 2d ou 3d. Supression routine modele_wave
! 2010.12 20-09-10  Possibilité de calcul en simple precision
! 2010.25 23-02-12  cas soulsby
!         23-03-12  cas soulsby: blindage bancs decouvrants
! S26     14-11-12  cas z1 proche de z0 revu
!         14-09-13  modif test sur kmax
!                   ajout cdb_2dh
!         02-11-14  retour en arriere par rapport a modif 14-11-12 pour eviter
!                   floating point error en cas de bancs decouvrants
!         07-07-15  spatialisation du z0 initial
!         17-07-16  seuil sur z1/z0 augmentE pour meilleur comportement en zone intertidale
!         13-02-17  fichier z0 aux normes du bathy maker
!         08-04-17  s-z coordinate
!         26-04-17  s-z coordinate suite
!         10-09-17  ajout z0b_land et zlevel_land
!         13-09-17  ajout de vel(:,:,kbottom,1) en projet 
!         09-11-17  ajout subroutine maxbotstress
!         20-11-17  allocate(maxbotstress_w(imax,jmax)) passe dans module_principal
!         05-12-17  amenagement grille verticale
!         23-08-18  suppression flag_nh2d
!         16-10-18  test sur iwve
! v267    16-11-19  cas iwve==2
! v292    07-11-20  - if(flag_merged_levels==1)then !11111111111> !07-11-20
!                   - cas iwve==3
! v298    08-03-21  tauc_ utilise velbot
! v348    27-05-22  cas if(flag_z0_macro==1)
!         28-05-22  la constante 10 remplacEe par dz_over_z0_min
! v349    29-06-22  ajout subroutines dragcoef_no_slip_condition3d et 2d
!         02-07-22  no slip viscosite proportionnelle A la taille de la maille
! v367    07-04-23  friction soulsby si cas 3
! v368    14-04-23  ajout z0_w dans lecture netcdf
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

!     if(kmax>1.and.flag_nh2d==0) then   ! 3d3d3d3d3d> !14-09-13
      if(kmax>1)                  then   ! 3d3d3d3d3d> !14-09-13 !23-08-18
! 3D model:
!       if(iwve==0.or. &
!          iwve==2.or. &
!          iwve==3)call dragcoef_logprofile ! Blumberg et Mellor 1987 !16-11-19
!       if(iwve==1)call dragcoef_soulsby    ! Uchiyama et al, JGR 2009  doi:10.1029/2008JC005135
        if(iwve==0.or. &
           iwve==2)call dragcoef_logprofile ! Blumberg et Mellor 1987 !16-11-19
        if(iwve==1.or. &       
           iwve==3)call dragcoef_soulsby    ! Uchiyama et al, JGR 2009 doi:10.1029/2008JC005135 !07-04-23

      else                               ! 3d3d3d3d3d>
! 2Dh model:

       call dragcoef_2dh

      endif                              ! 3d3d3d3d3d>

      end subroutine dragcoef

!.........................................................................

      subroutine dragcoef_logprofile
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='dragcoef_logprofile'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Si macro rugosite, preparer la fonction de repartition du bottom
! stress sur plusieur niveaux:
      if(flag_z0_macro==1)call dragcoef_macro_logprofile

! Ref: Blumberg et Mellor 1987.
!      const1=1.1
      if(flag_merged_levels==0)then !00000>
       do j=1,jmax
       do i=1,imax

        cdb_t(i,j)=max( (0.4/                        &
        log(max( (h_w(i,j)+depth_t(i,j,kmin_w(i,j))) & !04/03/03
                /z0_w(i,j) , dz_over_z0_min   ) )    & !02-11-14!17-07-16!28-05-22
                          )**2,cdseuil )                               
       enddo
       enddo
      endif                         !00000>

      if(flag_merged_levels==1)then !11111>
       do j=1,jmax
       do i=1,imax
! Methode 1
! Si la couche kmin_w devient infiniment petite le cd
! est maintenu A une valeur "continue" avec la couche kmin_w+1
        cdb_t(i,j)=max( (0.4/                       &
        log(max( 0.5*max(dz_t(i,j,kmin_w(i,j)  ,1)  &
                        ,dz_t(i,j,kmin_w(i,j)+1,1)) &
                /z0_w(i,j) , dz_over_z0_min   ) )   & !28-05-22 
                          )**2,cdseuil )                               
! Methode 2
! La couche kmerge et les couches inferieure ne font qu'une couche. La
! hauteur est la demi epaisseur de la couche totale: 0.5*( zw(kmerged+1)+h_w(1) )
! Note: si kmerged=kmin l'algo revient a l'algo habituel
!       cdb_t(i,j)=max( (0.4/                                    &
!       log(max( 0.5*(depth_w(i,j,kmerged_t(i,j)+1)+h_w(i,j)) &
!               /z0_w(i,j) , dz_over_z0_min    ) )                          & 
!                         )**2,cdseuil )                               
       enddo
       enddo
      endif                         !11111>


      end subroutine dragcoef_logprofile

!.........................................................................

      subroutine dragcoef_soulsby
      use module_principal ; use module_sedw
      implicit none
      real*4 tauw_,cdb_,tauc_
#ifdef synopsis
       subroutinetitle='dragcoef_soulsby'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     const1=1.1      !23-03-12
      if(flag_merged_levels==0)then !00000>
        do j=1,jmax ; do i=1,imax

       tauw_=0.5*rho*fw(i,j)*ubw(i,j)**2 ! Uchiyama et al, JGR 2009, Eq 9

       cdb_=(0.4/log(                                             &
       max( (h_w(i,j)+depth_t(i,j,kmin_w(i,j)))/z0_w(i,j) , dz_over_z0_min ) & !23-03-12!17-07-16!28-05-22
                      )  )**2

       tauc_=rho*cdb_*(0.5*(                                &
                   vel_u(i+1,j  ,kmin_u(i+1,j  ),0)**2      &
                  +vel_u(i  ,j  ,kmin_u(i  ,j  ),0)**2      &
                  +vel_v(i  ,j+1,kmin_v(i  ,j+1),0)**2      &
                  +vel_v(i  ,j  ,kmin_v(i  ,j  ),0)**2 ))

!      taub_=   tauc_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2) ! Uchiyama et al, JGR 2009, Eq 8
!           =rho*cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2)*|u|u
!           =rho*CDequivalent*|u|u with CDequivalent=cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2)
       cdb_t(i,j)=max(cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2),cdseuil)

! In momentum equations taubx/rho=cdb_t(i,j)*|u|u

        enddo ; enddo
      endif                         !00000>

      if(flag_merged_levels==1)then !11111> !05-12-17
        do j=1,jmax ; do i=1,imax

       tauw_=0.5*rho*fw(i,j)*ubw(i,j)**2 ! Uchiyama et al, JGR 2009, Eq 9

       cdb_=(0.4/log(                                          &
        max( 0.5*max(dz_t(i,j,kmin_w(i,j)  ,1)                 &
                    ,dz_t(i,j,kmin_w(i,j)+1,1))/z0_w(i,j),dz_over_z0_min) & !28-05-22
                      )  )**2

!      tauc_=rho*cdb_*(0.5*(                                &
!                  vel_u(i+1,j  ,kmin_u(i+1,j  ),0)**2      &
!                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),0)**2      &
!                 +vel_v(i  ,j+1,kmin_v(i  ,j+1),0)**2      &
!                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),0)**2 ))

! Cet algo favorise la continuite du champs:
       tauc_=rho*cdb_*( & !ooo> !08-03-21
                        (0.5*(velbot_u(i+1,j  )+velbot_u(i  ,j  )))**2 & 
                       +(0.5*(velbot_v(i  ,j+1)+velbot_v(i  ,j  )))**2 &
                      )   !ooo>

!      taub_=   tauc_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2) ! Uchiyama et al, JGR 2009, Eq 8
!           =rho*cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2)*|u|u
!           =rho*CDequivalent*|u|u with CDequivalent=cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2)
       cdb_t(i,j)=max(cdb_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2),cdseuil)

! In momentum equations taubx/rho=cdb_t(i,j)*|u|u

        enddo ; enddo
      endif                         !11111>

      end subroutine dragcoef_soulsby

!............................................................................

      subroutine dragcoef_2dh
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='dragcoef_2dh'
       subroutinedescription=' Drag coefficient 2D case'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       do j=1,jmax
       do i=1,imax
        cdb_t(i,j)=cdb_2dh                                         &    !14-09-13
       *4./max(1,  mask_u(i,j,kmax+1)+  mask_u(i+1,j  ,kmax+1)    &    !10/03/02
                +  mask_v(i,j,kmax+1)+  mask_v(i  ,j+1,kmax+1))        !10/03/02
       enddo
       enddo

      end subroutine dragcoef_2dh

!............................................................................

      subroutine dragcoef_initial
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='dragcoef_initial'
       subroutinedescription='Drag coefficient initial case'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
      if(flag_z0_macro==1) then !pmx> !27-05-22
       if(flag_merged_levels==1) stop &
       'CAS flag_z0_macro=1 & flag_merged_levels=1 EST A FAIRE'
       allocate(sigma_fric_wu(1:imax+1,1:jmax,kmax+1)) ; sigma_fric_wu=1. 
       allocate(sigma_fric_wv(1:imax,1:jmax+1,kmax+1)) ; sigma_fric_wv=1. 
      endif                     !pmx> !27-05-22

      if(kmax>2) then ! 3d3d3d3d3d>
! 3D model:
        call dragcoef_logprofile ! Blumberg et Mellor 1987
      else            ! 3d3d3d3d3d>
! 2Dh model:
       call dragcoef_2dh
      endif           ! 3d3d3d3d3d>


      end subroutine dragcoef_initial

!............................................................................

      subroutine dragcoef_macro_logprofile !27-05-22
      use module_principal ; use module_parallele
      implicit none

! Si macro rugosite, preparer la fonction de repartition du bottom
! stress sur plusieur niveaux:

! calculer (z_w+h)/2h1:
! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
! Note1: en empechant de diviser par plus que hz_w(i,j,1) on garantit que
!        sigma_fric_ww=1 quand k=kmax+1, autrement dit par tres petit
!        fond (encore plus petit que n*z0) on garantit que toute la
!        friction est bien prise en compte et bien repartie entre le fond
!        et la surface. 
! Note2: la multiplication par 2 de dz_over_z0_min vient de ce qu'on
!        considere le niveau "w" soit la facette, or la facette superieure de
!        la couche critique est A une distance du fond egale A 2 fois la
!        hauteur critique
! Note3: inutile de recaluler en k=kmax+1 car la valeur 1 a deja EtE
!        attribuEe A l'allocation
      do k=1,kmax+1 ; do j=1,jmax   ; do i=1,imax+1
       sigma_fric_wu(i,j,k)=min(1.,                                        &
                                     0.5*(depth_w(i  ,j,k)-depth_w(i  ,j,1)  &
                                         +depth_w(i-1,j,k)-depth_w(i-1,j,1)) &
                                  /min(2*dz_over_z0_min                      &
                                        *0.5*(z0_w(i  ,j)                    &
                                             +z0_w(i-1,j)),hz_u(i,j,1)))
!        if(i+par%timax(1)==2308.and.j+par%tjmax(1)==1347) &
!        write(10+par%rank,*)  &
!         real(sigma_fric_wu(i,j,k)) &
!        ,real(0.5*(depth_w(i  ,j,k)-depth_w(i  ,j,1)   &
!                  +depth_w(i-1,j,k)-depth_w(i-1,j,1))) &
!        ,real(hz_u(i,j,1))                          &
!        ,real(2*dz_over_z0_min                      &
!               *0.5*(z0_w(i  ,j)                    &
!                    +z0_w(i-1,j))),' u'
      enddo         ; enddo         ; enddo
      do k=1,kmax+1 ; do j=1,jmax+1 ; do i=1,imax
       sigma_fric_wv(i,j,k)=min(1.,                                        &
                                     0.5*(depth_w(i,j  ,k)-depth_w(i,j  ,1)  &
                                         +depth_w(i,j-1,k)-depth_w(i,j-1,1)) &
                                  /min(2*dz_over_z0_min                      &
                                        *0.5*(z0_w(i,j  )                    &
                                             +z0_w(i,j-1)),hz_v(i,j,1)))
!        if(i+par%timax(1)==2308.and.j+par%tjmax(1)==1347) &
!        write(10+par%rank,*)  &
!         real(sigma_fric_wv(i,j,k)) &
!        ,real(0.5*(depth_w(i,j  ,k)-depth_w(i,j  ,1)   &
!                  +depth_w(i,j-1,k)-depth_w(i,j-1,1))) &
!        ,real(hz_v(i,j,1))                          &
!        ,real(2*dz_over_z0_min                      &
!               *0.5*(z0_w(i  ,j)                    &
!                    +z0_w(i,j-1))),' v'
      enddo         ; enddo         ; enddo

      end subroutine dragcoef_macro_logprofile
!............................................................................

      subroutine dragcoef_initial_z0_file
      use module_principal ; use module_parallele
      implicit none
      integer ncid_
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='dragcoef_initial_z0'
       subroutinedescription='Initial Bottom Roughness'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! The z0 File is produced by xscan sym-tools : !07-06-15
! The z0 File has the same dim than bathycote_in.nc : !13-02-17

      if(par%rank==0)write(6,'(a,a)')'nf_open ',trim(texte250)

      status=nf_open(trim(texte250),nf_nowrite,ncid_)
      if(status/=0)stop 'Err 179 dragcoef_initial_z0 nf_open'

! definir/verifier varstart:
                    status=nf_inq_dimid(ncid_,'nx_t',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'ni_t',k0) !13-02-17   
       if(status/=0)stop 'z0 Err1890'

          status=nf_inq_dimlen(ncid_,k0,i1)             
       if(status/=0)stop 'z0 Err1891'

                    status=nf_inq_dimid(ncid_,'ny_t',k0)
       if(status/=0)status=nf_inq_dimid(ncid_,'nj_t',k0)          
       if(status/=0)stop 'z0 Err1892'

          status=nf_inq_dimlen(ncid_,k0,j1)             
       if(status/=0)stop 'z0 Err1893'

                    status=nf_inq_varid(ncid_,'z0-bottom_diffused_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'z0_roughness_bottom_t',var_id) 
       if(status/=0)status=nf_inq_varid(ncid_,'z0_w',var_id) !14-04-23
       if(status/=0)stop 'Err nf_inq_varid z0'

       status=nf_inq_vartype(ncid_,var_id,k0)                    
       if(k0/=nf_real)stop 'Err uncorrect vartype z0'
       
       ksecu=0
       if(i1==iglb.and.j1==jglb) then     !>>>
        ksecu=1
        varstart(1)=1+par%timax(1) ; varcount(1)=imax
        varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax
        status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                            ,varcount(1:2)   &
                                  ,anyvar2d(1:imax,1:jmax)) 
               z0_w(1:imax,1:jmax)=anyvar2d(1:imax,1:jmax)
       endif                              !>>>

       if(i1==iglb+2.and.j1==jglb+2) then !>>>
        ksecu=1
        varstart(1)=1+par%timax(1) ; varcount(1)=imax+2
        varstart(2)=1+par%tjmax(1) ; varcount(2)=jmax+2
        status=nf_get_vara_real(ncid_,var_id,varstart(1:2)   &
                                            ,varcount(1:2)   &
                               ,anyvar2d(0:imax+1,0:jmax+1))
        z0_w(0:imax+1,0:jmax+1)=anyvar2d(0:imax+1,0:jmax+1)
       endif                              !>>>

       if(ksecu==0)stop 'Err 230 dim dragcoef_initial_z0'

       if(status/=0)stop 'Err 198 nf_get_vara_real z0'

      status=nf_close(ncid_)

      if(status/=0)stop 'Err 2003 nf_close dragcoef_initial_z0'

      end subroutine dragcoef_initial_z0_file

!.........................................................................

      subroutine dragcoef_initial_z0_h !10-09-17
      use module_principal ; use module_parallele
      implicit none

! Par defaut pour ne pas laisser de point non initialisees si dimension
! de h_w ne collent pas a celles de z0_w
       z0_w=z0b
       do j=1,jmax ; do i=1,imax
        if(-h_w(i,j)>zlevel_land)z0_w(i,j)=z0b_land ! valeur des surfaces continentales
       enddo       ; enddo

      end subroutine dragcoef_initial_z0_h

!.........................................................................

      subroutine dragcoef_maxbotstress        !09-11-17
      use module_principal ; use module_sedw
      implicit none
      real*4 tauw_,tauc_,taum_,cdb_

      if(flag_maxbotstress==0)return
      if(iwve==0) stop 'Err 417 iwve=0 and flag_maxbotstress/=0' !16-10-18

! Note: on calcule une integralle dans le temps car module_offline archivera la moyenne
      if(flag_merged_levels==0)then !00000000000>

      do j=1,jmax ; do i=1,imax

       tauw_=0.5*rho*fw(i,j)*ubw(i,j)**2 ! Uchiyama et al, JGR 2009, Eq 9

       cdb_=(0.4/log(                                             &
       max( (h_w(i,j)+depth_t(i,j,kmin_w(i,j)))/z0_w(i,j) ,dz_over_z0_min) & !23-03-12!17-07-16!28-05-22
                      )  )**2

       tauc_=rho*cdb_*(0.5*(                                &
                   vel_u(i+1,j  ,kmin_u(i+1,j  ),0)**2      &
                  +vel_u(i  ,j  ,kmin_u(i  ,j  ),0)**2      &
                  +vel_v(i  ,j+1,kmin_v(i  ,j+1),0)**2      &
                  +vel_v(i  ,j  ,kmin_v(i  ,j  ),0)**2 ))

       taum_=tauc_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2) ! Uchiyama et al, JGR 2009, Eq 8

! Angle (radians) entre le courant et les vagues:
       x1=atan2(  &
                  vel_v(i  ,j+1,kmin_v(i  ,j+1),0)   & ! Angle (u,v) par rapport axe Oi en radians
                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),0)   & 
                , vel_u(i+1,j  ,kmin_u(i+1,j  ),0)   &
                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),0)   &
               )  &
            -dir_wave_t(i  ,j  ,1)                     ! - angle des vagues par rapport axe Oi en radian

! Equation 28, Aurelien Gangloff thesis
! Note: on calcule une integralle dans le temps car module_offline archivera la moyenne
       maxbotstress_w(i,j)=                                         &
       maxbotstress_w(i,j)+dti_now*( (taum_+tauw_*abs(cos(x1)))**2  &
                                    +(      tauw_*    sin(x1) )**2  )**0.5

       stresswave_w(i,j)=    &
       stresswave_w(i,j)+dti_now*tauw_
       stressc_w(i,j)= &
       stressc_w(i,j)+dti_now*tauc_

      enddo ; enddo

      endif                         !00000000000>

      if(flag_merged_levels==1)then !11111111111> !07-11-20

      do j=1,jmax ; do i=1,imax

       tauw_=0.5*rho*fw(i,j)*ubw(i,j)**2 ! Uchiyama et al, JGR 2009, Eq 9

       cdb_=(0.4/log(max(0.5*max(dz_t(i,j,kmin_w(i,j)  ,1)  &
                                ,dz_t(i,j,kmin_w(i,j)+1,1))/z0_w(i,j),dz_over_z0_min)))**2 !28-05-22 

! cet algo favorise les extremas de vitesse
!      tauc_=rho*cdb_*(0.5*(             &
!                  velbot_u(i+1,j  )**2  & 
!                 +velbot_u(i  ,j  )**2  & 
!                 +velbot_v(i  ,j+1)**2  & 
!                 +velbot_v(i  ,j  )**2))
! Cet algo favorise la continuite du champs:
       tauc_=rho*cdb_*( & !ooo>
                        (0.5*(velbot_u(i+1,j  )+velbot_u(i  ,j  )))**2 &
                       +(0.5*(velbot_v(i  ,j+1)+velbot_v(i  ,j  )))**2 &
                      )   !ooo>

       taum_=tauc_*(1.+ 1.2*(tauw_/(tauw_+tauc_+small2))**3.2) ! Uchiyama et al, JGR 2009, Eq 8

! Angle (radians) entre le courant et les vagues:
       x1=atan2( &
                   velbot_v(i  ,j+1)     & 
                  +velbot_v(i  ,j  )     &
                 , velbot_u(i+1,j  )     & 
                  +velbot_u(i  ,j  )     & 
               ) &
            -dir_wave_t(i  ,j  ,1)                     ! - angle des vagues par rapport axe Oi en radian

! Equation 28, Aurelien Gangloff thesis
! Note: on calcule une integralle dans le temps car module_offline archivera la moyenne
       maxbotstress_w(i,j)=                                         &
       maxbotstress_w(i,j)+dti_now*( (taum_+tauw_*abs(cos(x1)))**2  &
                                    +(      tauw_*    sin(x1) )**2  )**0.5

       stresswave_w(i,j)=    &
       stresswave_w(i,j)+dti_now*tauw_
       stressc_w(i,j)= &
       stressc_w(i,j)+dti_now*tauc_

      enddo ; enddo

      endif                         !11111111111>

      end subroutine dragcoef_maxbotstress        !09-11-17

!............................................................................

      subroutine dragcoef_initial_z0_local
      use module_principal ; use module_parallele
      implicit none

! Details dans
! https://docs.google.com/document/d/1cgBEAL0qvT7ernFu7FoxlNq8zwX7vhhwB5NNFeBgsb4/edit

      flag_stop=0

      if(index(texte90,'.ijvaldist')/=0) then ! m[°v°]m >

       ub2=ubound(z0_w) ; lb2=lbound(z0_w)
       if(par%rank==0)write(6,'(a,a)')'opening ',trim(texte90)
       open(unit=3,file=trim(texte90))

 2265    read(3,*,end=2263)i1,j1,i2,j2,x3,i10
         if(abs(i2-i1)>abs(j2-j1)) then !ooo>
          do i=i1,i2,sign(1,i2-i1)
           if(i2/=i1) then !>>>
              rap=real(i-i1)/real(i2-i1)
           else            !>>>
              rap=0.    
           endif           !>>>
           j=nint( rap*real(j2)+(1.-rap)*real(j1) )
           do j3=j-i10,j+i10 ; do i3=i-i10,i+i10
             x0=max((real(i10)-sqrt(real(i3-i)**2+real(j3-j)**2))/real(i10),0.)
             i4=i3-par%timax(1) ; j4=j3-par%tjmax(1)
             if(i4>=lb2(1).and.i4<=ub2(1).and.j4>=lb2(2).and.j4<=ub2(2)) then !pmx>
              z0_w(i4,j4)=max(z0_w(i4,j4),x0*x3+(1.-x0)*z0b)
             endif                                                            !pmx>
           enddo             ; enddo !i3,j3
          enddo ! i
         else                           !ooo>
          do j=j1,j2,sign(1,j2-j1)
           if(j2/=j1) then !>>>
              rap=real(j-j1)/real(j2-j1)
           else            !>>>
              rap=0.    
           endif           !>>>
           i=nint( rap*real(i2)+(1.-rap)*real(i1) )
           do j3=j-i10,j+i10 ; do i3=i-i10,i+i10
             x0=max((real(i10)-sqrt(real(i3-i)**2+real(j3-j)**2))/real(i10),0.)
             i4=i3-par%timax(1) ; j4=j3-par%tjmax(1)
             if(i4>=lb2(1).and.i4<=ub2(1).and.j4>=lb2(2).and.j4<=ub2(2)) then !pmx>
              z0_w(i4,j4)=max(z0_w(i4,j4),x0*x3+(1.-x0)*z0b)
             endif                                                            !pmx>
           enddo             ; enddo !i3,j3
          enddo ! j
         endif                          !ooo>

        goto 2265
 2263  close(3)

      else  ! m[°v°]m >

! Si l'extension du fichier n'est pas correcte STOP
       flag_stop=1

      endif ! m[°v°]m >

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'Err 498 dragcoef_initial_z0_local'

      end subroutine dragcoef_initial_z0_local

!.......................................................

      subroutine dragcoef_no_slip_condition_3d !29-06-22
      use module_principal ; use module_parallele
      implicit none

! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing

! Note:
! Si 2-mask_f(i,j)-mask_f(i,j+1)=0 pas de cote voisine (pas de frottement)
! Si 2-mask_f(i,j)-mask_f(i,j+1)=1 Une cote voisine (frottement)
! Si 2-mask_f(i,j)-mask_f(i,j+1)=0 Deux cotes voisines, c.a.d. canal (double frottement)
      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
! viscosite constante
!      vel_u(i,j,k,2)=vel_u(i,j,k,2)/(1.+coastal_viscosity*(invdy_u(i,j)**2)*(2-mask_f(i,j,kmax)-mask_f(i,j+1,kmax))*dti_lp)
! viscosite proportionnelle A la taille de la maille: ! !02-07-22
       vel_u(i,j,k,2)=vel_u(i,j,k,2)/(1.+coastal_viscosity*invdy_u(i,j)*(2-mask_f(i,j,kmax)-mask_f(i,j+1,kmax))*dti_lp)
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
! viscosite constante
!      vel_v(i,j,k,2)=vel_v(i,j,k,2)/(1.+coastal_viscosity*(invdx_v(i,j)**2)*(2-mask_f(i,j,kmax)-mask_f(i+1,j,kmax))*dti_lp)
! viscosite proportionnelle A la taille de la maille: ! !02-07-22
       vel_v(i,j,k,2)=vel_v(i,j,k,2)/(1.+coastal_viscosity*invdx_v(i,j)*(2-mask_f(i,j,kmax)-mask_f(i+1,j,kmax))*dti_lp)
      enddo
      enddo
      enddo

      end subroutine dragcoef_no_slip_condition_3d

!.......................................................

      subroutine dragcoef_no_slip_condition_2d !29-06-22
      use module_principal ; use module_parallele
      implicit none

! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing

      do j=2,jmax-1
      do i=2,imax
! viscosite constante
!      velbar_u(i,j,2)=velbar_u(i,j,2)/(1.+coastal_viscosity*(invdy_u(i,j)**2)*(2-mask_f(i,j,kmax)-mask_f(i,j+1,kmax))*dte_lp)
! viscosite proportionnelle A la taille de la maille: ! !02-07-22
       velbar_u(i,j,2)=velbar_u(i,j,2)/(1.+coastal_viscosity*invdy_u(i,j)*(2-mask_f(i,j,kmax)-mask_f(i,j+1,kmax))*dte_lp)
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
! viscosite constante
!      velbar_v(i,j,2)=velbar_v(i,j,2)/(1.+coastal_viscosity*(invdx_v(i,j)**2)*(2-mask_f(i,j,kmax)-mask_f(i+1,j,kmax))*dte_lp)
! viscosite proportionnelle A la taille de la maille: ! !02-07-22
       velbar_v(i,j,2)=velbar_v(i,j,2)/(1.+coastal_viscosity*invdx_v(i,j)*(2-mask_f(i,j,kmax)-mask_f(i+1,j,kmax))*dte_lp)
      enddo
      enddo

      end subroutine dragcoef_no_slip_condition_2d
