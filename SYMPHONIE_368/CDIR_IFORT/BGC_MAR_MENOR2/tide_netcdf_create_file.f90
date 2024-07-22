










      subroutine tide_netcdf_create_file                                !16-12-09
!______________________________________________________________________
! SYMPHONIE ocean model
! release 296 - last update: 24-02-21
!______________________________________________________________________
!...............................................................................
! Version date      Description des modifications
! 2010.2  24-12-09  - Version separee de la routine tide_analysis
!                   - boucle sur ktide à l'interieur
!                   - comme les parametres nodaux ne sont plus integrés dans
!                     les amplitudes & phase initiales, il ne faut plus les
!                     enlever dans cette etape.
!         26-12-09  on ajoute le potentiel de marée total dans le fichier
!         27-12-09  modif sur nom du fichier
! 2010.3  04-01-10  - debug attribut "associate"
!                   - des noms de variables exactement comme en sortie de T-UGO
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.8  07-05-10  attribut 'axis' remplacé par 'content'
! 2010.11 04-07-10  amelioration de l'analyse:
!                       - utilisation de tableaux distincts de ceux du forcage
!                       - la chronologie time.dat prend en compte l'evolution
!                         des parametres nodaux
!                       - Possibilite d'une analyse "glissante" melangeant
!                         l'analyse de la periode en cours avec l'analyse
!                         precedente. Permet une reactualisation reguliere
!                         de l'analyse, avec archivage de toutes les echeances
!                         dans un fichier netcdf
!         15-07-10  "status" declaré dans module_principal
! 2010.12 20-09-10  Possibilité de calcul en simple precision
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes
! 2010.24 28-10-11  Ne plus ajouter le nbre d'iteration au nom du fichier d'analyse
! 2010.25 26-02-12  supression mtde et ntde
! 2010.25 21-03-12  bienvenue au tableau ub*, lb* pour remplacer
!                    ubound et lbound qui posent un problème avec pgf 10.6
!         28-03-12  echange compatible avec pgf Stelios
!         03-04-12  routine pour attributs generaux des fichiers netcdf
! S26     08-02-13  ecrire avec la filiere commune
!         15-09-13  un nom correct pour l'archivage du potentiel pour eviter
!                   erreurs si fichiers repris pour initialisation
!         16-11-13  Ecrire LSA dans le fichier de reanalyse pour pouvoir
!                   forcer le modele avec le fichier de reanalyse
!         05-06-15  Message a l'ecran du nom du fichier produit
!         08-12-15  dans le fichier de reanalyse de la maree par S26 on
!                   archive egalement les champs forcants correspondant
!         10-04-16  Message ecran
!         23-07-17  modif nom de champs netcdf
!         17-04-18  ajout call tide_ssh_bias_correction('rmv') 
!         18-11-18  Hg-ObcHg renommE Hg_ObcHg (idem pour Ha) pour compatibilite FERRET
!         28-11-18  texte80(11)='mask_obc_z1' pour les champs S-FES
! v276    04-04-20  analyse harmonique 3d
! v282    08-05-20  fichier netcdf: correction unite densite
!                   ajout de l'analyse de la moyenne temporelle
! v283    15-05-20  Flux energie potentielle
! v285    07-06-20  call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20
! v288    08-09-20  var_addoffset=rho
! v296    24-02-21  signal maree interne dans la SSH
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      use module_principal ; use module_parallele ; use module_s
      use pnetcdf
      implicit none
      integer ichoix,looproc_,looprocmax_
      character*6 type_

! Enlever l'eventuelle Correction de la SSH dans la partie forcage pour bien ecrire
! le vrai forcage FES et non pas le forcage corrige
       do ktide=1,kmaxtide
         if(index(nametide(ktide,6),'none')==0)call tide_ssh_bias_correction('rmv') !17-04-18
       enddo

! Calcul d'erreur RMS sur les amplitudes complexes sur tout le domaine et sur les OBC seulement
! on note que les tableaux de forcage viennent d'etre debarasses de leur correction eventuelle et retrouve
! donc la valeur (interpolee) de FES
       if(par%rank==0)open(unit=3,file='tmp/messages',position='append')
       if(par%rank==0) then
        write(3,'(a)')'----------------------------------------------------'
        write(3,'(a)')'RMS Error (meters) of the complex amplitude of tides'
        write(3,'(a)')'     RMS full domain  RMS obc'
        write(6,'(a)')'RMS Error (meters) of the complex amplitude'
        write(6,'(a)')'     RMS full domain RMS obc'
       endif
       do ktide=1,kmaxtide !KKKKKKKKKK>
        sum1=0. ; sum2=0 ; sum3=0. ; sum4=0.
        do j=1,jmax ; do i=1,imax !IJ-IJ>
         x1=mask_t(i,j,kmax)*mask_i_w(i)*mask_j_w(j)
         x2=(sshtidecos_w(i,j,ktide)-sshtidecosout_w(i,j,ktide))**2 &
           +(sshtidesin_w(i,j,ktide)-sshtidesinout_w(i,j,ktide))**2
         sum1=sum1+x1 ; sum2=sum2+x1*x2
         if(i+par%timax(1)==1.or.i+par%timax(1)==iglb.or.j+par%tjmax(1)==1.or.j+par%tjmax(1)==jglb) then
          sum3=sum3+x1 ; sum4=sum4+x1*x2
         endif
        enddo       ; enddo       !IJ-IJ>
        call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
        if(par%rank==0) then
         if(index(nametide(ktide,1),'Sreanalysis')/=0) then !>>>
       write(6,*)'WARNING RMS is relative to a Sreanalysis forcing file'
       write(3,*)'WARNING RMS is relative to a Sreanalysis forcing file'
         endif                                              !>>>
         write(6,'(a3,1x,2(f7.4,10x))')shortnametide(ktide),real(sqrt(sum2glb/max(sum1glb,1.d0))),real(sqrt(sum4glb/max(sum3glb,1.d0)))
         write(3,'(a3,1x,2(f7.4,10x))')shortnametide(ktide),real(sqrt(sum2glb/max(sum1glb,1.d0))),real(sqrt(sum4glb/max(sum3glb,1.d0)))
        endif
       enddo               !KKKKKKKKKK>
       if(par%rank==0)close(3)

! Archiver la latitude et la longitude en:
!     type_='real'
      type_='double'

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do ktide=1,kmaxtidep1 ! kmaxtide
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(par%rank==0) then !m[°u°]m> 
       k=s_unit(7)
       open(unit=k,file='output_file_extension')
        texte30='' 
        read(k,*,end=126)texte30
  126  close(k)
      endif                !m[°u°]m> 
      call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20

      if(ktide<=kmaxtide) then !ooo>
       texte250=trim(nametide(ktide,4))
      else                     !ooo> !09-05-20
       k1=index(nametide(1,4),'/',back=.true.)
       if(ktide==kmaxtide+1)texte250=nametide(1,4)(1:k1)//'LF0.nc'
       if(ktide==kmaxtide+2)texte250=nametide(1,4)(1:k1)//'LF1.nc'
      endif                    !ooo>

      k0=index(texte250,'.nc')
      texte250=texte250(1:k0-1)//'_2D'//trim(texte30)//'.nc'
      if(par%rank==0)write(6,'(a,a)')'Tidal 2D reanalysis in file ',trim(texte250)

      do loop_netcdf=0,1

      count_netcdfvar=0

      filval=-9999.
      texte80(3)='none'
      texte80(4)='none'

      if(loop_netcdf==0) then  !§§§§§§§>
        status=nfmpi_create(par%comm2d,texte250, nf_clobber + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      else                     !§§§§§§§>
        status=nfmpi_open(par%comm2d,texte250, nf_write + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop erreur tide netcdf 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! write variables dimensions:
      call graph_out_trueaxis

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)   ! 'real'

      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_u'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_v'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_u'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_v'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=h_w(i,j)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3)='undisturbed water depth'                  ! long_name
      texte80(4)='undisturbed_water_depth'                  !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_w')

!Ha S26
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !ooo>
       !Cas standard
        do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*sqrt(sshtidecosout_w(i,j,ktide)**2+sshtidesinout_w(i,j,ktide)**2) &
                  +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
       else                     !ooo>
       !Cas particulier de la moyenne temporelle !09-05-20
        do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*0.5*(sshtidecosout_w(i,j,ktide)+sshtidesinout_w(i,j,ktide)) &
                  +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
       endif                    !ooo>
      endif                  !--------->
      texte80(1)='Ha' ; texte80(2)='m'      ! variable ; units
      texte80(3)='ssh_tide_amplitude'  ! long_name
      if(ktide==kmaxtide+1)texte80(3)='ssh_time_averaged'
      texte80(4)=texte80(3)                 ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_t')

! OiVa S26
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !ooo>
       !Cas standard
        do j=1,jmax ; do i=1,imax+1
         anyvar2d(i,j)=mask_u(i,j,kmax+1)*sqrt(veltidecosout_u(i,j,ktide)**2+veltidesinout_u(i,j,ktide)**2) &
                   +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo
       else                     !ooo>
       !Cas particulier de la moyenne temporelle !09-05-20
        do j=1,jmax ; do i=1,imax+1
         anyvar2d(i,j)=mask_u(i,j,kmax+1)*0.5*(veltidecosout_u(i,j,ktide)+veltidesinout_u(i,j,ktide)) &
                   +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo
       endif                    !ooo>
      endif                  !--------->
      texte80(1)='OiVa' ; texte80(2)='m/s'  ! variable ; units
      texte80(3)='along_i_axis_current_amplitude'     ! long_name
      if(ktide==kmaxtide+1) &
      texte80(3)='along_i_axis_current_time_averaged' ! long_name
      texte80(4)=texte80(3)                 ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_u')

! OjVa S26
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !ooo>
       !Cas standard
        do j=1,jmax+1 ; do i=1,imax
         anyvar2d(i,j)=mask_v(i,j,kmax+1)*sqrt(veltidecosout_v(i,j,ktide)**2+veltidesinout_v(i,j,ktide)**2) &
                   +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo
       else                     !ooo>
       !Cas particulier de la moyenne temporelle !09-05-20
        do j=1,jmax+1 ; do i=1,imax
         anyvar2d(i,j)=mask_v(i,j,kmax+1)*0.5*(veltidecosout_v(i,j,ktide)+veltidesinout_v(i,j,ktide)) &
                   +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo
       endif                    !ooo>
      endif                  !--------->
      texte80(1)='OjVa' ; texte80(2)='m/s'  ! variable ; units
      texte80(3)='along_j_axis_current_amplitude'     ! long_name
      if(ktide==kmaxtide+1) &
      texte80(3)='along_j_axis_current_time_averaged' ! long_name
      texte80(4)=texte80(3)                 ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_v')

! CAS PARTICULIER DE LA MOYENNE: PAS DE PHASE
      if(ktide/=kmaxtide+1) then !PHASE-PHASE>

      if(loop_netcdf==1) then !--------->
! Flux energie potentielle Oi !15-05-20
        x0=0.25*rho*grav/(1000.*10*1000.) ! normalisation par  l'ordre
                                          ! de grandeur de rho*grav*dx
        do j=1,jmax ; do i=1,imax
            anyvar2d(i,j)=mask_t(i,j,kmax+1)*x0                 &
                                          *h_w(i  ,j)           &
                                         *dy_t(i  ,j)*(         &
                               sshtidecosout_w(i  ,j,ktide)     &
                            *( veltidecosout_u(i  ,j,ktide)     &
                              +veltidecosout_u(i+1,j,ktide) )   &

                              +sshtidesinout_w(i  ,j,ktide)     &
                            *( veltidesinout_u(i  ,j,ktide)     &
                              +veltidesinout_u(i+1,j,ktide) ) ) &

                     +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
        anyvar2d(0,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0)=filval ; anyvar2d(:,jmax+1)=filval
      endif                  !--------->
      texte80(1)='flux_oi' ; texte80(2)='10**7*kg*m**2/s**3'      ! variable ; units
      texte80(3)='rho*grav*h*dy*<ssh*u>/10**7'  ! long_name
      texte80(4)=texte80(3)                          ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
! Flux energie potentielle Oj !15-05-20
        x0=0.25*rho*grav/(1000.*10*1000.) ! normalisation par  l'ordre
                                          ! de grandeur de rho*grav*dx
        do j=1,jmax ; do i=1,imax
            anyvar2d(i,j)=mask_t(i,j,kmax+1)*x0                 &
                                          *h_w(i,j  )           &
                                         *dx_t(i,j  )*(         &
                               sshtidecosout_w(i,j  ,ktide)     &
                            *( veltidecosout_v(i,j  ,ktide)     &
                              +veltidecosout_v(i,j+1,ktide) )   &

                              +sshtidesinout_w(i,j  ,ktide)     &
                            *( veltidesinout_v(i,j  ,ktide)     &
                              +veltidesinout_v(i,j+1,ktide) ) ) &

                     +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
        anyvar2d(0,:)=filval ; anyvar2d(imax+1,:)=filval ; anyvar2d(:,0)=filval ; anyvar2d(:,jmax+1)=filval
      endif                  !--------->
      texte80(1)='flux_oj' ; texte80(2)='10**7*kg*m**2/s**3'      ! variable ; units
      texte80(3)='rho*grav*h*dy*<ssh*v>/10**7'  ! long_name
      texte80(4)=texte80(3)             ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_t')

      if(rhp_zavr_xy==1.and.kmax>1.and.ktide/=kmaxtide+1) then !3DCASE> !24-02-21
!Ha SS3D
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*sqrt(ssh3dtidecosout_w(i,j,ktide)**2+ssh3dtidesinout_w(i,j,ktide)**2) &
                  +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
      endif                  !--------->
      texte80(1)='HaSSH3D' ; texte80(2)='m'      ! variable ; units
      texte80(3)='ssh3d_tide_amplitude'  ! long_name
      if(ktide==kmaxtide+1)texte80(3)='ssh3d_time_averaged'
      texte80(4)=texte80(3)                 ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_t')

! Hg SSH3D
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        x1=-atan2(-ssh3dtidesinout_w(i,j,ktide),ssh3dtidecosout_w(i,j,ktide))
        x1=mod( x1*180./pi,360.*un)
        if(x1>180.)x1=x1-360.
        anyvar2d(i,j)=mask_t(i,j,kmax+1)*x1 &
                  +(1-mask_t(i,j,kmax+1))*filval
       enddo ; enddo
      endif                  !--------->
      texte80(1)='HgSSH3D' ; texte80(2)='degree'                   ! variable ; units
      texte80(3)='ssh3d tide phase lag'        ! long_name
      texte80(4)='ssh3d_tide_phase_lag'                         !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

! Ha SSH2D
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*sqrt((sshtidecosout_w(i,j,ktide)-ssh3dtidecosout_w(i,j,ktide))**2  &
                                              +(sshtidesinout_w(i,j,ktide)-ssh3dtidesinout_w(i,j,ktide))**2) &
                  +(1.-mask_t(i,j,kmax+1))*filval
        enddo ; enddo
      endif                  !--------->
      texte80(1)='HaSSH2D' ; texte80(2)='m'      ! variable ; units
      texte80(3)='ssh2d_tide_amplitude'  ! long_name
      if(ktide==kmaxtide+1)texte80(3)='ssh2d_time_averaged'
      texte80(4)=texte80(3)                 ! standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_t')

! Hg SSH2D
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        x1=-atan2(-(sshtidesinout_w(i,j,ktide)-ssh3dtidesinout_w(i,j,ktide)),(sshtidecosout_w(i,j,ktide)-ssh3dtidecosout_w(i,j,ktide)))
        x1=mod( x1*180./pi,360.*un)
        if(x1>180.)x1=x1-360.
        anyvar2d(i,j)=mask_t(i,j,kmax+1)*x1 &
                  +(1-mask_t(i,j,kmax+1))*filval
       enddo ; enddo
      endif                  !--------->
      texte80(1)='HgSSH2D' ; texte80(2)='degree'                   ! variable ; units
      texte80(3)='ssh2d tide phase lag'        ! long_name
      texte80(4)='ssh2d_tide_phase_lag'                         !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')
      endif           !3DCASE>

! Hg S26
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        x1=-atan2(-sshtidesinout_w(i,j,ktide),sshtidecosout_w(i,j,ktide))
        x1=mod( x1*180./pi,360.*un)
        if(x1>180.)x1=x1-360.
        anyvar2d(i,j)=mask_t(i,j,kmax+1)*x1 &
                  +(1-mask_t(i,j,kmax+1))*filval
       enddo ; enddo
      endif                  !--------->
      texte80(1)='Hg' ; texte80(2)='degree'                   ! variable ; units
      texte80(3)='ssh tide phase lag'                         ! long_name
      texte80(4)='ssh_tide_phase_lag'                         !  standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

! OiVg S26
      if(loop_netcdf==1) then !--------->
        do j=1,jmax ; do i=1,imax+1
         x1=-atan2(-veltidesinout_u(i,j,ktide),veltidecosout_u(i,j,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar2d(i,j)=mask_u(i,j,kmax+1)*x1 &
                   +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo
      endif                  !--------->
      texte80(1)='OiVg' ; texte80(2)='degree'                 ! variable ; units
      texte80(3)='along i axis current phase lag'                     ! long_name
      texte80(4)='along_i_axis_current_phase_lag'                     !  standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_u')

! OjVg S26
      if(loop_netcdf==1) then !--------->
        do j=1,jmax+1 ; do i=1,imax
         x1=-atan2(-veltidesinout_v(i,j,ktide),veltidecosout_v(i,j,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar2d(i,j)=mask_v(i,j,kmax+1)*x1 &
                   +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo
      endif                  !--------->
      texte80(1)='OjVg' ; texte80(2)='degree'                 ! variable ; units
      texte80(3)='along j axis current phase lag'                     ! long_name
      texte80(4)='along_j_axis_current_phase_lag'                     !  standard_name
      texte80(5)='YX'  ; texte80(7)='real'
      call netcdf_main('_v')

      endif                      !PHASE-PHASE>

      if(ktide>kmaxtide) goto 409

      if(tideforces==0.or.tideforces==2) then !LSALSALSA> !09-05-20

      if(loop_netcdf==1) then !--------->
! Oter le potentiel astronomique pour n'ecrire que le LSA: !16-11-13
      x0=0.
      x1=0.
      x2=0.
      if(nutide(ktide).eq.0)x0=1.
      if(nutide(ktide).eq.1)x1=1.
      if(nutide(ktide).eq.2)x2=1.

      const1=0.7*equitide(ktide)
      if(tideforces.eq.2)const1=0.
      if(tideforces.eq.3)const1=0.

      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
! Note les signes devant const1 sont coherents avec initial_tide
! x3: LSA cos
         x3=potidecos_w(i,j,ktide)                                    &
         -const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 &
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    cos(nutide(ktide)*lon_t(i,j))

! Note les signes devant const1 sont coherents avec initial_tide
! x4: LSA sin
         x4=potidesin_w(i,j,ktide)                                    &
         +const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 &
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    sin(nutide(ktide)*lon_t(i,j))

!       anyvar2d(i,j)=sqrt(potidecos_w(i,j,ktide)**2+potidesin_w(i,j,ktide)**2)
        anyvar2d(i,j)=sqrt(x3**2+x4**2)
       else
        anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='LSAa' ; texte80(2)='m'                    ! variable ; units
      texte80(3)='LSA amplitude'           ! long_name
      texte80(4)='LSA_amplitude'           ! standard_name
      texte80(5)='YX' ; texte80(7)='real'
      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
! Oter le potentiel astronomique pour n'ecrire que le LSA: !16-11-13
      x0=0.
      x1=0.
      x2=0.
      if(nutide(ktide).eq.0)x0=1.
      if(nutide(ktide).eq.1)x1=1.
      if(nutide(ktide).eq.2)x2=1.

      const1=0.7*equitide(ktide)
      if(tideforces.eq.2)const1=0.
      if(tideforces.eq.3)const1=0.

      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then

! Note les signes devant const1 sont coherents avec initial_tide
! x3: LSA cos
         x3=potidecos_w(i,j,ktide)                                    &
         -const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 &
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    cos(nutide(ktide)*lon_t(i,j))

! Note les signes devant const1 sont coherents avec initial_tide
! x4: LSA sin
         x4=potidesin_w(i,j,ktide)                                    &
         +const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 &
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    sin(nutide(ktide)*lon_t(i,j))

!       x5=-atan2(-potidesin_w(i,j,ktide),potidecos_w(i,j,ktide)) !24-12-09
        x5=-atan2(-x4,x3)
        x5=mod( x5*180./pi,360.*un)
        if(x5>180.)x5=x5-360.
        anyvar2d(i,j)=x5
       else
        anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='LSAg' ; texte80(2)='degree'               ! variable ; units
      texte80(3)='LSA phase lag'           ! long_name
      texte80(4)='LSA_phase_lag'           !  standard_name
      texte80(5)='YX'        ; texte80(7)='real'
      call netcdf_main('_t')

      endif                                   !LSALSALSA>

!Ha FES
       if(loop_netcdf==1) then !---------> !08-12-15
       do j=0,jmax+1
       do i=0,imax+1
        if(mask_t(i,j,kmax+1)==1) then
         anyvar2d(i,j)=sqrt(sshtidecos_w(i,j,ktide)**2+sshtidesin_w(i,j,ktide)**2)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_Ha' ; texte80(2)='m'              ! variable ; units
       texte80(3:4)='forcing_ssh_tide_amplitude'             !  standard_name
       texte80(5)='YX'  ; texte80(7)='real'
       call netcdf_main('_t')

! Ha S26 - Ha FES
       if(loop_netcdf==1) then !---------> !08-12-15
       do j=0,jmax+1
       do i=0,imax+1
        if(mask_t(i,j,kmax+1)==1) then
         anyvar2d(i,j)=sqrt(sshtidecosout_w(i,j,ktide)**2+sshtidesinout_w(i,j,ktide)**2) &
                      -sqrt(sshtidecos_w(i,j,ktide)**2   +sshtidesin_w(i,j,ktide)**2)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='Ha_ObcHa' ; texte80(2)='m'  ! variable ; units
       texte80(3:4)='Ha_ObcHa'             !  standard_name !23-07-17
       texte80(5)='YX'  ; texte80(7)='real'
       texte80(11)='mask_obc_z1' ! l'obc etant traitee de maniere differente pour symphonie et son forcage, la couronne en i=0, i=imax+1 etc...  est masquee !28-11-18
       call netcdf_main('_t')

!Hg FES
       if(loop_netcdf==1) then !---------> !08-12-15
       do j=0,jmax+1
       do i=0,imax+1
        if(mask_t(i,j,kmax+1)==1) then
         x1=-atan2(-sshtidesin_w(i,j,ktide),sshtidecos_w(i,j,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar2d(i,j)=x1
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_Hg' ; texte80(2)='degree'           ! variable ; units
       texte80(3:4)='forcing_ssh_tide_phase_lag'               !  standard_name
       texte80(5)='YX' ; texte80(7)='real'
       call netcdf_main('_t')

! Hg S26 - Hg FES
       if(loop_netcdf==1) then !---------> !08-12-15
       do j=0,jmax+1
       do i=0,imax+1
        if(mask_t(i,j,kmax+1)==1) then
         x1=-atan2(-sshtidesinout_w(i,j,ktide),sshtidecosout_w(i,j,ktide)) &
            +atan2(-sshtidesin_w(i,j,ktide),sshtidecos_w(i,j,ktide))
         x1=mod( x1*180./pi,360.d0)
         if(x1> 180.)x1=x1-360.
         if(x1<-180.)x1=x1+360.
         anyvar2d(i,j)=x1
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='Hg_ObcHg' ; texte80(2)='degree'           ! variable ; units
       texte80(3:4)='Hg_ObcHg'               !  standard_name !23-07-17
       texte80(5)='YX' ; texte80(7)='real'
       texte80(11)='mask_obc_z1' ! l'obc etant traitee de maniere differente pour symphonie et son forcage, la couronne en i=0, i=imax+1 etc...  est masquee !28-11-18
       call netcdf_main('_t')

! OiVa FES
       if(loop_netcdf==1) then !---------> !08-12-15
       do j=1,jmax
       do i=1,imax+1
        if(mask_u(i,j,kmax+1)==1) then
         anyvar2d(i,j)=sqrt(veltidecos_u(i,j,ktide)**2+veltidesin_u(i,j,ktide)**2)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_OiVa' ; texte80(2)='m/s'                    ! variable ; units
       texte80(3:4)='forcing_along_i_axis_current_amplitude'                     !  standard_name
       texte80(5)='YX'  ; texte80(7)='real'
       call netcdf_main('_u')


! OiVg FES
       if(loop_netcdf==1) then !--------->
       do j=1,jmax
       do i=1,imax+1
        if(mask_u(i,j,kmax+1)==1) then
         x1=-atan2(-veltidesin_u(i,j,ktide),veltidecos_u(i,j,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar2d(i,j)=x1
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_OiVg' ; texte80(2)='degree'                 ! variable ; units
       texte80(3:4)='along_i_axis_current_phase_lag'                     !  standard_name
       texte80(5)='YX'  ; texte80(7)='real'
       call netcdf_main('_u')

! OjVa FES
       if(loop_netcdf==1) then !--------->
       do j=1,jmax+1
       do i=1,imax
        if(mask_v(i,j,kmax+1)==1) then
         anyvar2d(i,j)=sqrt(veltidecos_v(i,j,ktide)**2+veltidesin_v(i,j,ktide)**2)
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_OjVa' ; texte80(2)='m/s'                    ! variable ; units
       texte80(3:4)='forcing_along_j_axis_current_amplitude'                     !  standard_name
       texte80(5)='YX'  ; texte80(7)='real'
       call netcdf_main('_v')

! OjVg FES
       if(loop_netcdf==1) then !--------->
       do j=1,jmax+1
       do i=1,imax
        if(mask_v(i,j,kmax+1)==1) then
        x1=-atan2(-veltidesin_v(i,j,ktide),veltidecos_v(i,j,ktide))
        x1=mod( x1*180./pi,360.*un)
        if(x1>180.)x1=x1-360.
        anyvar2d(i,j)=x1
        else
         anyvar2d(i,j)=filval
        endif
       enddo
       enddo
       endif                  !--------->
       texte80(1)='forcing_OjVg' ; texte80(2)='degree'                 ! variable ; units
       texte80(3:4)='forcing_along_j_axis_current_phase_lag'            !  standard_name
       texte80(5)='YX'  ; texte80(7)='real'
       call netcdf_main('_v')

! Hcos S26 - Hcos FES
       if(loop_netcdf==1) then !---------> !17-04-18
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*(sshtidecosout_w(i,j,ktide)-sshtidecos_w(i,j,ktide)) &
                   +(1-mask_t(i,j,kmax+1))*filval
       enddo ; enddo
       endif                   !--------->
       texte80(1)='Hcos-ObcHcos' ; texte80(2)='m' ; texte80(3:4)=texte80(1) ; texte80(5)='YX' ; texte80(7)='real'
       texte80(11)='mask_obc_z1' ! l'obc etant traitee de maniere differente pour symphonie et son forcage, la couronne en i=0, i=imax+1 etc...  est masquee !28-11-18
       call netcdf_main('_t')

! Hsin S26 - Hsin FES
       if(loop_netcdf==1) then !---------> !17-04-18
       do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=mask_t(i,j,kmax+1)*(sshtidesinout_w(i,j,ktide)-sshtidesin_w(i,j,ktide)) &
                   +(1-mask_t(i,j,kmax+1))*filval
       enddo ; enddo
       endif                   !--------->
       texte80(1)='Hsin-ObcHsin' ; texte80(2)='m' ; texte80(3:4)=texte80(1) ; texte80(5)='YX' ; texte80(7)='real'
       texte80(11)='mask_obc_z1' ! l'obc etant traitee de maniere differente pour symphonie et son forcage, la couronne en i=0, i=imax+1 etc...  est masquee !28-11-18
       call netcdf_main('_t')

  409 continue

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

      call netcdf_general_attributes(ncid)

! Definition of variables: done.
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!     endif                         !procprocprocproc>
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!     call barriere(looproc_   ,3)
!      enddo ! fin de boucle sur looproc_

      enddo  ! fin de boucle sur loop_netcdf


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      enddo ! fin de boucke sur ktide
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Analyse du courant 3d et de la densite potentielle:
      if(flag_tide3d_analysis==1)call tide_netcdf_create_3d_file

      end subroutine tide_netcdf_create_file

!...............................................................................

      subroutine tide_netcdf_create_3d_file      !04-04-20                      
      use module_principal ; use module_parallele ; use module_s
      use pnetcdf
      implicit none
      integer ichoix,looproc_,looprocmax_
      character*6 type_

      if(flag_tide3d_analysis/=1)return

! Archiver la latitude et la longitude en:
!     type_='real'
      type_='double'

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do ktide=1,kmaxtidep1 ! kmaxtide
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(par%rank==0) then !m[°u°]m> !05-05-17
       k=s_unit(7)
       open(unit=k,file='output_file_extension')
        texte30='' !22-10-19
        read(k,*,end=471)texte30
  471  close(k)
      endif                !m[°u°]m> 
      call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20

      if(ktide<=kmaxtide) then !ooo>
       texte250=trim(nametide(ktide,4))
      else                     !ooo> !09-05-20
       k1=index(nametide(1,4),'/',back=.true.)
       if(ktide==kmaxtide+1)texte250=nametide(1,4)(1:k1)//'LF0.nc'
       if(ktide==kmaxtide+2)texte250=nametide(1,4)(1:k1)//'LF1.nc'
      endif                    !ooo>

      k0=index(texte250,'.nc')
      texte250=texte250(1:k0-1)//'_3D'//trim(texte30)//'.nc'
      if(par%rank==0)write(6,'(a,a)')'Tidal 3D reanalysis in file ',trim(texte250)

      do loop_netcdf=0,1

      count_netcdfvar=0

      filval=-9999.
      texte80(3)='none'
      texte80(4)='none'

      if(loop_netcdf==0) then  !§§§§§§§>
        status=nfmpi_create(par%comm2d,texte250, nf_clobber + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      else                     !§§§§§§§>
        status=nfmpi_open(par%comm2d,texte250, nf_write + NF_64BIT_OFFSET , MPI_INFO_NULL, ncid) !10-09-12
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop erreur tide netcdf 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! write variables dimensions:
      call graph_out_trueaxis

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)   ! 'real'

      call netcdf_main('_t')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_u'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lon_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lon_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='longitude_v'
      texte80(2)='degree_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_t(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_t(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_t')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_u(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_u(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_u'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_u')


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=lat_v(i,j)*rad2deg
        anyv3d(i,j,1,1)=lat_v(i,j)*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='latitude_v'
      texte80(2)='degree_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX'
      texte80(7)=trim(type_)
      call netcdf_main('_v')


! Oi3dVa 
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !pmx>
       !Cas standard
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         anyvar3d(i,j,k)=mask_u(i,j,kmax+1)*sqrt(vel3dtidecosout_u(i,j,k,ktide)**2+vel3dtidesinout_u(i,j,k,ktide)**2) &
                     +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       else                     !pmx>
       ! Cas particulier de la moyenne temporelle !09-05-20
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         anyvar3d(i,j,k)=mask_u(i,j,kmax+1)*0.5*(vel3dtidecosout_u(i,j,k,ktide)+vel3dtidesinout_u(i,j,k,ktide)) &
                     +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       endif                    !pmx>
      endif                  !--------->
      texte80(1)='Oi3dVa' ; texte80(2)='m/s'               ! variable ; units
      texte80(3)='3Dalong_i_axis_current_amplitude'                     ! long_name
      if(ktide==kmaxtide+1) &
      texte80(3)='3Dalong_i_axis_current_time_averaged'  ! long_name
      texte80(4)=texte80(3)            !  standard_name
      texte80(5)='ZYX'  ; texte80(7)='real'
      call netcdf_main('_u')

! Oj3dVa 
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !pmx>
       !Cas standard
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         anyvar3d(i,j,k)=mask_v(i,j,kmax+1)*sqrt(vel3dtidecosout_v(i,j,k,ktide)**2+vel3dtidesinout_v(i,j,k,ktide)**2) &
                     +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       else                     !pmx>
       ! Cas particulier de la moyenne temporelle !09-05-20
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         anyvar3d(i,j,k)=mask_v(i,j,kmax+1)*0.5*(vel3dtidecosout_v(i,j,k,ktide)+vel3dtidesinout_v(i,j,k,ktide)) &
                     +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       endif                    !pmx>
      endif                  !--------->
      texte80(1)='Oj3dVa' ; texte80(2)='m/s'    ! variable ; units
      texte80(3)='3Dalong_j_axis_current_amplitude'     ! long_name
      if(ktide==kmaxtide+1)  &
      texte80(3)='3Dalong_j_axis_current_time_averaged' ! long_name
      texte80(4)=texte80(3)                     !  standard_name
      texte80(5)='ZYX'  ; texte80(7)='real'
      call netcdf_main('_v')

! RHPa
      if(loop_netcdf==1) then !--------->
       if(ktide/=kmaxtide+1) then !pmx>
       !Cas standard
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         anyvar3d(i,j,k)=mask_t(i,j,kmax+1)*sqrt(rhptidecosout_t(i,j,k,ktide)**2+rhptidesinout_t(i,j,k,ktide)**2) &
                     +(1-mask_t(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       else                     !pmx>
       ! Cas particulier de la moyenne temporelle !09-05-20
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         anyvar3d(i,j,k)=mask_t(i,j,kmax+1)*0.5*(rhptidecosout_t(i,j,k,ktide)+rhptidesinout_t(i,j,k,ktide)) &
                     +(1-mask_t(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
       endif                    !pmx>
      endif                  !--------->
      texte80(1)='RHPa' ; texte80(2)='kg/m**3'      ! variable ; units !08-05-20
      texte80(3)='PotentialDensity_tide_amplitude'  ! long_name
      if(ktide==kmaxtide+1)  &
      texte80(3)='PotentialDensity_tide_time_averaged'  ! long_name
      texte80(4)=texte80(3)                     !  standard_name
      texte80(5)='ZYX'  ; texte80(7)='real'
      var_addoffset=0.                       !08-09-20
      if(ktide==kmaxtide+1)var_addoffset=rho ! Remettre A zero apres call netcdf_main !08-09-20
      call netcdf_main('_t')
      var_addoffset=0.                       ! remise A la valeur par defaut !08-09-20

! CAS PARTICULIER DE LA MOYENNE: PAS DE PHASE
      if(ktide/=kmaxtide+1) then !PHASE-PHASE>

! Oi3dVg 
      if(loop_netcdf==1) then !--------->
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         x1=-atan2(-vel3dtidesinout_u(i,j,k,ktide),vel3dtidecosout_u(i,j,k,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar3d(i,j,k)=mask_u(i,j,kmax+1)*x1 &
                     +(1-mask_u(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
      endif                  !--------->
      texte80(1)='Oi3dVg' ; texte80(2)='degree'                 ! variable ; units
      texte80(3)='3Dalong i axis current phase lag'                     ! long_name
      texte80(4)='3Dalong_i_axis_current_phase_lag'                     !  standard_name
      texte80(5)='ZYX'  ; texte80(7)='real'
      call netcdf_main('_u')

! Oj3dVg
      if(loop_netcdf==1) then !--------->
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         x1=-atan2(-vel3dtidesinout_v(i,j,k,ktide),vel3dtidecosout_v(i,j,k,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar3d(i,j,k)=mask_v(i,j,kmax+1)*x1 &
                     +(1-mask_v(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
      endif                  !--------->
      texte80(1)='Oj3dVg' ; texte80(2)='degree'                 ! variable ; units
      texte80(3)='3Dalong j axis current phase lag'                     ! long_name
      texte80(4)='3Dalong_j_axis_current_phase_lag'                     !  standard_name
      texte80(5)='ZYX'  ; texte80(7)='real'
      call netcdf_main('_v')

! RHPg
      if(loop_netcdf==1) then !--------->
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         x1=-atan2(-rhptidesinout_t(i,j,k,ktide),rhptidecosout_t(i,j,k,ktide))
         x1=mod( x1*180./pi,360.*un)
         if(x1>180.)x1=x1-360.
         anyvar3d(i,j,k)=mask_t(i,j,kmax+1)*x1 &
                     +(1-mask_t(i,j,kmax+1))*filval
        enddo ; enddo ; enddo
      endif                  !--------->
      texte80(1)='RHPg' ; texte80(2)='degree'                   ! variable ; units
      texte80(3)='PotentialDensity tide phase lag'                         ! long_name
      texte80(4)='PotentialDensity_tide_phase_lag'                         !  standard_name
      texte80(5)='ZYX' ; texte80(7)='real'
      call netcdf_main('_t')

      endif                      !PHASE-PHASE>

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

      call netcdf_general_attributes(ncid)

! Definition of variables: done.
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!     endif                         !procprocprocproc>
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!     call barriere(looproc_   ,3)
!      enddo ! fin de boucle sur looproc_

      enddo  ! fin de boucle sur loop_netcdf


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      enddo ! fin de boucke sur ktide
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine tide_netcdf_create_3d_file
