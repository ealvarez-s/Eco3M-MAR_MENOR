      module module_biobc
!______________________________________________________________________
! S model
! release S.26 - last update: 09-12-14
!______________________________________________________________________
! S.26   27-04-13  ajout d'une zone de rappel vers la solution exterieure
!        28-04-13  pas de valeurs negatives (sauf filval)
!        01-05-13  Ne pas appliquer 2 fois le rappel + debug
!        02-05-13  Interpoler sur la zone de rappel
!        01-07-14  modif pour conservation de mpi
!        09-12-14  appliquer les conditions selon la valeur de obcstatus
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
      integer biobc_imax,biobc_jmax,biobc_kmax       &
       ,var_dims,var_type,tabdim(4)                  &
       ,biobczoom_istr,biobczoom_jstr                &
       ,biobczoom_iend,biobczoom_jend
      integer ::                                     &
       biobclistline                                 &
      ,biobc_init_flag=0                             &
      ,biobc_prevfile=1                              &
      ,biobc_nextfile=2                              &
      ,biobc_landvalue_flag=0                        &
      ,biobc_extforcing=0
!     ,bio_relax_size
      double precision,dimension(:,:),allocatable :: &
       biobc_lon                                     &
      ,biobc_lat
       real*4,dimension(:,:,:),allocatable :: &
       biobc_depth                            &
      ,biobc_var
      double precision lat_1_1_gridt,lat_1_1_gridu,lat_1_1_gridv &
                      ,lon_1_1_gridt,lon_1_1_gridu,lon_1_1_gridv &
                      ,lon_2_1_gridt,lat_1_2_gridt
      real*4,dimension(:,:),allocatable :: ij2biobc_i,ij2biobc_j
      real*4,dimension(:),allocatable :: biobc_z1dv
      integer(kind=1) :: flag_biosponge_size=0 !24-07-18

! Elements a inclure dans le fichier restart:
! debut module
!     real*4,dimension(:,:,:,:,:),allocatable :: &
!      bio_relax_north                           &
!     ,bio_relax_south                           &
!     ,bio_relax_east                            &
!     ,bio_relax_west
!     double precision      &
!      biobc_nextfiletime   &
!     ,biobc_prevfiletime
!     character*200         &
!      biobc_nextfilename   &
!     ,biobc_prevfilename
!     real*4                &
!      relax_bio
! fin module

contains

!...........................................................................

      subroutine biobc_driver
      implicit none
#ifdef synopsis
       subroutinetitle='biobc_driver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(initial_main_status==0) then !>>>

         call biobc_initial

      else                            !>>>

       if(elapsedtime_now>biobc_nextfiletime) then !----->
        call biobc_whoisnext
        call biobc_interp_step1(biobc_nextfile)    ! Interpolation sur la grille s du prochain fichier
       endif                                       !----->

       rap_biobc=(   elapsedtime_now-biobc_prevfiletime)    &
                /(biobc_nextfiletime-biobc_prevfiletime)

      endif                           !>>>

      end subroutine biobc_driver

!...........................................................................

      subroutine biobc_initial
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='biobc_initial'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! What kind of initialisation?
! if biobc_type==2 then downscaling case
      do vb=1,vbmax
       if(biobc_type(vb)==2)biobc_extforcing=1
      enddo
      if(biobc_extforcing==0)return

      biobc_init_flag=1                       ! Declarer la phase d'initialisation ouverte

      bio_relax_size=sponge_l+1

!     k0=0
!     if(bio_relax_size>imax-1)k0=1 ! On assure avec imax-1 car imax serait sans doute suffisant
!     if(bio_relax_size>jmax-1)k0=1
!     if(k0==1) then !>>>>>>>!01-07-14
!      k0=1
!      write(6,'(a,a,a,i0)') &
!      ,'Biogeochemical sponge layer can not be larger than'   &
!      ,' one sub-domain. Take a smaller value for sponge_l.'  &
!      ,' A first guess suggestion is:',min(imax-1,jmax-1)
!      stop 'STOP biobc_initial'
!     endif          !>>>>>>>!01-07-14

! Ici on choisit entre la procedure standard (flag_biosponge_size=0) qui cherche A economiser la memoire vive et la
! procedure eponge plus large que la taille des sous-domaines mpi (flag_biosponge_size=1)
      flag_biosponge_size=0 ! cas standard avec des petites tableaux 
      if(bio_relax_size>imax-1)flag_biosponge_size=1 ! On assure avec imax-1 car imax serait sans doute suffisant
      if(bio_relax_size>jmax-1)flag_biosponge_size=1 ! 24-07-18
! Si sur au moins 1 domaine flag_biosponge_size=1 alors tous les domaines doivent avoir flag_biosponge_size=1
      call mpi_allreduce(flag_biosponge_size,k0,1,mpi_integer,mpi_max,par%comm2d ,ierr)
      flag_biosponge_size=k0

      call biobc_allocate_relax

      call biobc_extract_zone                 ! definir les bornes min et max de la zone d'extraction
                                              ! tableaux de passage de la grille s à la grille de forcage
      call biobc_whoisnext                    ! Identification des fichiers correspondant a la date en cours
      call biobc_interp_step1(biobc_prevfile) ! Interpolation sur la grille s du fichier anterieur
      call biobc_interp_step1(biobc_nextfile) ! Interpolation sur la grille s du prochain fichier

      biobc_init_flag=0                       ! Declarer la phase d'initialisation terminee

      end subroutine biobc_initial

!...........................................................................

      subroutine biobc_extract_zone
      implicit none
      integer ncid_
#ifdef synopsis
       subroutinetitle='biobc_extract_zone'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      open(unit=3,file=trim(directory)//txtslash//'biogrdfilelist')
      read(3,'(a)')texte250
      close(3)

      status=nf_open(trim(texte250),nf_nowrite,ncid_)
      if(status/=0)stop 'erreur nf_open biobc_grid'

                    status=nf_inq_dimid(ncid_,'ni_t',dim_x_id)
       if(status/=0)stop 'erreur nf_inq_dimid x biobc'

                    status=nf_inq_dimid(ncid_,'nj_t',dim_y_id)
       if(status/=0)stop 'erreur nf_inq_dimid y biobc'

                    status=nf_inq_dimid(ncid_,'nk_t',dim_z_id)
       if(status/=0)stop 'erreur nf_inq_dimid z biobc'

      status=nf_inq_dimlen(ncid_,dim_x_id,biobc_imax)
      status=nf_inq_dimlen(ncid_,dim_y_id,biobc_jmax)
      status=nf_inq_dimlen(ncid_,dim_z_id,biobc_kmax)

      allocate(biobc_lon(biobc_imax,biobc_jmax))
      allocate(biobc_lat(biobc_imax,biobc_jmax))

! lire les longitudes:
                   status=nf_inq_varid(ncid_,'longitude_t',var_id)
      if(status/=0)stop 'Erreur var_id variable biobc longitude'

      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)


      if(var_type/=nf_double)stop 'erreur type biobc_lon'
      if(var_dims/=2)stop 'biobc lon var_dims/=2'
      i1=biobc_imax ; j1=biobc_jmax
      varstart(1)=1 ; varcount(1)=biobc_imax
      varstart(2)=1 ; varcount(2)=biobc_jmax

       status=nf_get_vara_double(ncid_,var_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,biobc_lon(1:i1,1:j1))
      if(status/=0)stop 'Erreur lecture longitude biobc'

! lire les latitudes:
                   status=nf_inq_varid(ncid_,'latitude_t',var_id)
      if(status/=0)stop 'Erreur var_id variable biobc latitude'

      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type/=nf_double)stop 'Erreur type sur latitude biobc'
      if(var_dims/=2)stop 'biobc lat var_dims/=2'
      i1=biobc_imax ; j1=biobc_jmax
      varstart(1)=1 ; varcount(1)=biobc_imax
      varstart(2)=1 ; varcount(2)=biobc_jmax

       status=nf_get_vara_double(ncid_,var_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,biobc_lat(1:i1,1:j1))
      if(status/=0)stop 'Erreur lecture latitude biobc'

      call biobc_sgrid_to_extgrid             ! tableaux de passage de la grille s à la grille de forcage

       biobczoom_istr=999999 ; biobczoom_iend=-999999 ! first guess
       biobczoom_jstr=999999 ; biobczoom_jend=-999999 ! first guess
       do j=0,jmax+1
       do i=0,imax+1
         i1=nint(ij2biobc_i(i,j))
         j1=nint(ij2biobc_j(i,j))
         biobczoom_istr=min(biobczoom_istr,i1-1)
         biobczoom_jstr=min(biobczoom_jstr,j1-1)
         biobczoom_iend=max(biobczoom_iend,i1+1)
         biobczoom_jend=max(biobczoom_jend,j1+1)
       enddo
       enddo

! on elargit un peu plus pour boucher les trous sans être restreint par la taille
! reduite de la zone d'extraction, et puis aussi pour rattraper l'erreur liee au fait
! que plusieurs grilles biobc (point u, v, t) peuvent être presentes.
       i0=10 ! elargissement à 10 lignes / 10 colonnes supplementaires - repere 1233
       biobczoom_istr=max(biobczoom_istr-i0,1)
       biobczoom_iend=min(biobczoom_iend+i0,biobc_imax)
       biobczoom_jstr=max(biobczoom_jstr-i0,1)
       biobczoom_jend=min(biobczoom_jend+i0,biobc_jmax)

       biobc_imax=biobczoom_iend-biobczoom_istr+1
       biobc_jmax=biobczoom_jend-biobczoom_jstr+1

      deallocate(biobc_lon)
      deallocate(biobc_lat)

! Charger la profondeur
      allocate(biobc_depth(biobc_imax,biobc_jmax,biobc_kmax))
                   status=nf_inq_varid(ncid_,'depth_t',var_id)
      if(status/=0)stop 'Erreur var_id variable biobc depth'
      varstart(1)=biobczoom_istr ; varcount(1)=biobc_imax
      varstart(2)=biobczoom_jstr ; varcount(2)=biobc_jmax
      varstart(3)=1              ; varcount(3)=biobc_kmax
       status=nf_get_vara_real(ncid_,var_id,varstart(1:3)   &
                                           ,varcount(1:3)   &
       ,biobc_depth(1:biobc_imax,1:biobc_jmax,1:biobc_kmax))
      if(status/=0) &
      stop 'Erreur nf_get_vara_real biobc_depth biobc_extract_zone'


! Fermer le fichier netcdf:
      status=nf_close(ncid_)

      end subroutine biobc_extract_zone

!----------------------------------------------------------------------------

      subroutine biobc_sgrid_to_extgrid
      implicit none
      include 'netcdf.inc'
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
#ifdef synopsis
       subroutinetitle='biobc_sgrid_to_extgrid'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       lon_1_1_gridt=biobc_lon(1,1) ; lat_1_1_gridt=biobc_lat(1,1)
       lon_2_1_gridt=biobc_lon(2,1) ; lat_1_2_gridt=biobc_lat(1,2)

      allocate(ij2biobc_i(0:imax+1,0:jmax+1))
      allocate(ij2biobc_j(0:imax+1,0:jmax+1))


      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      x2=real(biobc_imax/2)
      x3=real(biobc_jmax/2)
      deci=x2
      decj=x3

! ETAPE 1: trouver les coordonnées dans la grille ORCA:

! First guess: centre du domaine:
      k10=0
 1456 continue

! Principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*Di+dlon/dj*Dj=Dlon
!      dlat/di*Di+dlat/dj*Dj=Dlat
! On cherche Di et Dj correspondant à Dlon=lon(i,j)-lonbiobc(i0,j0)
!                                et à Dlat=lat(i,j)-latbiobc(i0,j0)

      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1
      i0=i1+i2
      j0=j1+j2

      dlon_di_=biobc_lon(i0+1,j0  )-biobc_lon(i0-1,j0  )
      dlon_dj_=biobc_lon(i0  ,j0+1)-biobc_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)-biobc_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *(biobc_lat(i0  ,j0+1)-biobc_lat(i0  ,j0-1))          &
          -(biobc_lat(i0+1,j0  )-biobc_lat(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

      deci_(i2,j2)=min(max(                                  &
       i0+( dlon_dm_                                            &
          *(biobc_lat(i0  ,j0+1)-biobc_lat(i0,j0-1))            &
          -(rad2deg*lat_t(i,j)  -biobc_lat(i0,j0))              &
          *dlon_dj_)/x1*0.5    &
                   ,2.0001d0),biobc_imax-1.0001d0)

      decj_(i2,j2)=min(max(                                  &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)  -biobc_lat(i0  ,j0))            &
          -(biobc_lat(i0+1,j0  )-biobc_lat(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5    &
                   ,2.0001d0),biobc_jmax-1.0001d0)

      enddo
      enddo

      rapi=deci-i1
      rapj=decj-j1

      deci=(1.-rapi)*(1.-rapj)*deci_(0,0)   &
          +(1.-rapi)*    rapj *deci_(0,1)   &
          +    rapi *    rapj *deci_(1,1)   &
          +    rapi *(1.-rapj)*deci_(1,0)
      decj=(1.-rapi)*(1.-rapj)*decj_(0,0)   &
          +(1.-rapi)*    rapj *decj_(0,1)   &
          +    rapi *    rapj *decj_(1,1)   &
          +    rapi *(1.-rapj)*decj_(1,0)

! Si le point visé est different du first guess refaire le calcul
! avec un first guess donné par le dernier point visé:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       if(k10>20)then !!!!!!>
         if(par%rank==0)write(6,*)'par%rank=',par%rank
         if(par%rank==0)write(6,*)'(i,j)=   ',i,j
         if(par%rank==0)write(6,*)'deci decj',deci,decj
         stop 'hr_to_lr ne converge pas dans le forfait'
       endif          !!!!!!>
       goto 1456
      endif

      ij2biobc_i(i,j)=deci
      ij2biobc_j(i,j)=decj

      enddo
      enddo

!     deallocate(biobc_lon)
!     deallocate(biobc_lat)

      end subroutine biobc_sgrid_to_extgrid

!..........................................................

      subroutine biobc_whoisnext
      implicit none
      integer loop_,year_,month_,day_,hour_,minute_,second_
#ifdef synopsis
       subroutinetitle='biobc_whoisnext'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(biobc_init_flag==1)biobclistline=0 ! initialisation

       open(unit=3,file=trim(directory)//txtslash//'biovarfilelist')

        do loop_=1,biobclistline ! En phase iterative faire defiler
!       read(3,'(a)')texte250    ! (rapidement) le debut du fichier
        read(3,*)
        enddo

  319   read(3,'(a)')texte250
        biobclistline=biobclistline+1
        k=index(texte250,'.nc ')
        read(texte250(k-15:k-12),*)year_
        read(texte250(k-11:k-10),*)month_
        read(texte250(k-9:k-8),*)day_
        read(texte250(k-6:k-5),*)hour_
        read(texte250(k-4:k-3),*)minute_
        read(texte250(k-2:k-1),*)second_


        biobc_prevfilename=biobc_nextfilename
        biobc_nextfilename=texte250

        biobc_prevfiletime=biobc_nextfiletime
        call datetokount(year_,month_,day_,hour_,minute_,second_)
        biobc_nextfiletime=elapsedtime_out

        rap_biobc=(   elapsedtime_now-biobc_prevfiletime)    &
                 /(biobc_nextfiletime-biobc_prevfiletime)

        if(elapsedtime_out<=elapsedtime_now) goto 319
       close(3)

      end subroutine biobc_whoisnext

!..........................................................

      subroutine biobc_interp_step1(time_)
      implicit none
      real*4 filval_
      integer ncid_,time_,vb_
#ifdef synopsis
       subroutinetitle='biobc_interp_step1'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0) then !-------->
      write(6,'(a,a)')'fichier bio previous ',trim(biobc_prevfilename)
      write(6,'(a,a)')'fichier bio next     ',trim(biobc_nextfilename)
      write(6,*)'previous time ',biobc_prevfiletime
      write(6,*)'next     time ',biobc_nextfiletime
      endif                !-------->

      if(.not.allocated(biobc_var))allocate(biobc_var(biobc_imax,biobc_jmax,biobc_kmax))
      if(.not.allocated(biobc_z1dv))allocate(biobc_z1dv(biobc_kmax))

      status=-999
      if(time_==biobc_nextfile)status=nf_open(trim(biobc_nextfilename),nf_nowrite,ncid_)
      if(time_==biobc_prevfile)status=nf_open(trim(biobc_prevfilename),nf_nowrite,ncid_)
      if(status/=0)stop 'erreur nf_open biobc_interp_step1'

      do vb_=1,vbmax ; vb=vb_
!     do vb =1,vbmax

       if(biobc_type(vb)==2) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       status=nf_inq_varid(ncid_,trim(biovarname(vb)),var_id)
       if(status/=0)stop 'erreur nf_inq_varid biobc_interp_step1'

       varstart(1)=biobczoom_istr ; varcount(1)=biobc_imax
       varstart(2)=biobczoom_jstr ; varcount(2)=biobc_jmax
       varstart(3)=1              ; varcount(3)=biobc_kmax
       varstart(4)=1              ; varcount(4)=1

       status=nf_get_vara_real(ncid_,var_id,varstart(1:4)   &
                                           ,varcount(1:4)   &
       ,biobc_var(1:biobc_imax,1:biobc_jmax,1:biobc_kmax))

        if(status/=0) then !.........>
         write(6,'(a,a)')'Pb variable ',trim(biovarname(vb))
         stop 'Erreur nf_get_vara_double biobc_interp_step1'
        endif              !.........>
!.........
! Traitement des valeurs masquées:
       status=nf_get_att_real(ncid_,var_id,'_FillValue',filval_)
       if(status/=0)  &
       stop 'erreur nf_get_att_int _Fillvalue biobc_interp_step1'


! Pas de valeurs negative (sauf filval) ! 28-04-13
! On commence par la surface:
        k=biobc_kmax
        do j=1,biobc_jmax
        do i=1,biobc_imax
         if(biobc_var(i,j,k)/=filval_) &
            biobc_var(i,j,k)=max(biobc_var(i,j,k),0.) !28-04-13
        enddo
        enddo

! Sous le fond extrapolation gradient vertical nul
        do k=biobc_kmax-1,1,-1
        do j=1,biobc_jmax
        do i=1,biobc_imax
         if(biobc_var(i,j,k)==filval_) then !----->
            biobc_var(i,j,k)=biobc_var(i,j,k+1)
         else                               !----->
            biobc_var(i,j,k)=max(biobc_var(i,j,k),0.) !28-04-13
         endif                              !----->
        enddo
        enddo
        enddo

! Dans le continent bouchage avec la valeur non masquee la plus proche:
        if(biobc_landvalue_flag==0)call biobc_landvalue_step1(filval_) ! identifier les trous
                                   call biobc_landvalue_step2          ! les boucher

!.........

! Compute the interpolation:
! En mpi la cpu est indexee sur le plus lent. Inutile donc d'optimiser la cpu
! d'un sous-domaine particulier si les autres vont moins vite (et imposent qu'on
! les attendent). L'interpolation est donc realisee sur tout le domaine. Sans
! compter que les regles de non-chevauchement sont plus complexes en mpi (selon
! qu'un sous-domaine a un voisin ou une frontiere ouverte). On commente donc
! les anciennes lignes pour revenir a une interpolation globale plus simple
! et finalement pas plus couteuse (et plus sure)
         call biobc_interp_step2(0,imax+1,1,0,jmax+1,1,time_) ! produit avnyvar3d
         call biobc_interp_mpiobc                             ! echange z0 sur anyvar3d

         call biobc_interp_step3(time_)

!       if(biobc_init_flag==1) then !-- Initialisation = full grid -->
!        call biobc_interp_step2(0,imax+1,1,0,jmax+1,1,time_)
!        call biobc_interp_step3(time_)
!       endif                       !-- Initialisation = full grid -->

!       if(biobc_init_flag==0) then !-- Iterative phase = Lateral Boundaries Only -->

!        i0=0 ; i1=imax+1 ; j0=0 ; j1=bio_relax_size  !02-05-13
!        call biobc_interp_step2(i0,i1,1,j0,j1,1,time_) ! OBC j=0 j=jmax+1

!        i0=0 ; i1=imax+1 ; j0=jmax-bio_relax_size+1 ; j1=jmax+1
!        call biobc_interp_step2(i0,i1,1,j0,j1,1,time_) ! OBC j=0 j=jmax+1

!        i0=0 ; i1=bio_relax_size ; j0=bio_relax_size+1 ; j1=jmax-bio_relax_size
!        call biobc_interp_step2(i0,i1,1,j0,j1,1,time_) ! OBC j=0 j=jmax+1

!        i0=imax-bio_relax_size+1 ; i1=imax+1 ; j0=bio_relax_size+1 ; j1=jmax-bio_relax_size
!        call biobc_interp_step2(i0,i1,1,j0,j1,1,time_) ! OBC j=0 j=jmax+1

!        call biobc_interp_step3(time_)
!       endif                       !-- Iterative phase = Lateral Boundaries Only -->

       endif                      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      enddo ! fin de boucle sur vb_

      status=nf_close(ncid_)

      end subroutine biobc_interp_step1

!..........................................................

      subroutine biobc_interp_step2(istr_,iend_,idlt_,jstr_,jend_,jdlt_,time_)
      implicit none
      integer istr_,iend_,idlt_,jstr_,jend_,jdlt_,time_
#ifdef synopsis
       subroutinetitle='biobc_interp_step2'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        do j=jstr_,jend_,jdlt_   ! 0,jmax+1
        do i=istr_,iend_,idlt_   ! 0,imax+1

        i1=int(ij2biobc_i(i,j))-biobczoom_istr+1 !01-07-14
        j1=int(ij2biobc_j(i,j))-biobczoom_jstr+1

        rapi=ij2biobc_i(i,j)-int(ij2biobc_i(i,j))
        rapj=ij2biobc_j(i,j)-int(ij2biobc_j(i,j))

         do k2=1,biobc_kmax
          biobc_z1dv(k2)= (1.-rapi)*(1.-rapj)*biobc_depth(i1  ,j1  ,k2) &
                         +(1.-rapi)*    rapj *biobc_depth(i1  ,j1+1,k2) &
                         +    rapi *(1.-rapj)*biobc_depth(i1+1,j1  ,k2) &
                         +    rapi *    rapj *biobc_depth(i1+1,j1+1,k2)
         enddo

         k0=1
         do k=1,kmax

           do k2=k0,biobc_kmax-1
            if(depth_t(i,j,k)>=biobc_z1dv(k2).and. &
               depth_t(i,j,k)< biobc_z1dv(k2+1)) then !>>>>>
               k1=k2
               rapk=(depth_t(i,j,k)-biobc_z1dv(k2))   &
                 /(biobc_z1dv(k2+1)-biobc_z1dv(k2))
               goto 449
            endif                                     !>>>>>
           enddo ! fin de boucle sur k2
           if(depth_t(i,j,k)< biobc_z1dv(1))then          !--->
            k1=1 ; rapk=0.  ! gradient nul
           endif                                          !--->
           if(depth_t(i,j,k)>=biobc_z1dv(biobc_kmax))then !--->
            k1=biobc_kmax-1  ; rapk=1.
           endif                                          !--->
  449      k0=k1

          anyvar3d(i,j,k)=                                          &
          ((1.-rapi)*(1.-rapj)*(1.-rapk)*biobc_var(i1  ,j1  ,k1)    &
          +(1.-rapi)*    rapj *(1.-rapk)*biobc_var(i1  ,j1+1,k1)    &
          +    rapi *(1.-rapj)*(1.-rapk)*biobc_var(i1+1,j1  ,k1)    &
          +    rapi *    rapj *(1.-rapk)*biobc_var(i1+1,j1+1,k1)    &
          +(1.-rapi)*(1.-rapj)*    rapk *biobc_var(i1  ,j1  ,k1+1)  &
          +(1.-rapi)*    rapj *    rapk *biobc_var(i1  ,j1+1,k1+1)  &
          +    rapi *(1.-rapj)*    rapk *biobc_var(i1+1,j1  ,k1+1)  &
          +    rapi *    rapj *    rapk *biobc_var(i1+1,j1+1,k1+1))

         enddo ! fin de boucle sur k
        enddo  ! fin de boucle sur i
        enddo  ! fin de boucle sur j

      end subroutine biobc_interp_step2

!..........................................................

      subroutine biobc_interp_mpiobc
      implicit none
      integer loop_,idi_anyvar3d_z0
#ifdef synopsis
       subroutinetitle='biobc_interp_mpiobc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour mpi conservation il faut echanger anyvar3d
! mpi boundaries:
#ifdef parallele
       call get_type_echange('z0','anyvar3d_z0'  &
                                  ,anyvar3d      &
                           ,lbound(anyvar3d)     &
                           ,ubound(anyvar3d)     &
                              ,idi_anyvar3d_z0)

      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(anyvar3d,idi_anyvar3d_z0,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif


      end subroutine biobc_interp_mpiobc

!..........................................................

      subroutine biobc_interp_step3(time_)
      implicit none
      integer time_
#ifdef synopsis
       subroutinetitle='biobc_interp_step3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Assign tracer arrays with tmp buffer array anyvar3d:

      if(biobc_init_flag==0) then !-- iterative phase -->

! time Shift:

! Tableaux de conditions aux limites du schema d'advection upstream:
       do k=1,kmax
        do i=0,imax+1
         biobc_i_t(i,k,vb,1,1)=biobc_i_t(i,k,vb,1,2)
         biobc_i_t(i,k,vb,2,1)=biobc_i_t(i,k,vb,2,2)
        enddo
        do j=0,jmax+1
         biobc_j_t(j,k,vb,1,1)=biobc_j_t(j,k,vb,1,2)
         biobc_j_t(j,k,vb,2,1)=biobc_j_t(j,k,vb,2,2)
        enddo
       enddo

! Tableaux pour la couche de rappel vers la solution externe
      if(flag_biosponge_size==0) then !m°v°m> !24-07-18
       do k=1,kmax
       do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
        do i=1,imax
         bio_relax_north(i,j0,k,vb,1)=bio_relax_north(i,j0,k,vb,2)
         bio_relax_south(i,j0,k,vb,1)=bio_relax_south(i,j0,k,vb,2)
        enddo
       enddo
       enddo
       do k=1,kmax
       do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
        do j=1,jmax
         bio_relax_east(j,i0,k,vb,1)=bio_relax_east(j,i0,k,vb,2)
         bio_relax_west(j,i0,k,vb,1)=bio_relax_west(j,i0,k,vb,2)
        enddo
       enddo
       enddo
      endif                           !m°v°m> !24-07-18

      if(flag_biosponge_size==1) then !mOvOm> !24-07-18
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
         bio_relax_full(i,j,k,vb,1) &
        =bio_relax_full(i,j,k,vb,2)
       enddo       ; enddo       ; enddo
      endif                           !mOvOm> !24-07-18

      endif                       !-- iterative phase -->

! Tableaux de conditions aux limites du schema d'advection upstream:
      i1=0 ; i2=imax+1 ; j1=0 ; j2=jmax+1
      do k=1,kmax
       do i=0,imax+1
          biobc_i_t(i,k,vb,1,time_)=anyvar3d(i,j1,k)
          biobc_i_t(i,k,vb,2,time_)=anyvar3d(i,j2,k)
       enddo
       do j=0,jmax+1
         biobc_j_t(j,k,vb,1,time_)=anyvar3d(i1,j,k)
         biobc_j_t(j,k,vb,2,time_)=anyvar3d(i2,j,k)
       enddo
      enddo

! Tableaux pour la couche de rappel vers la solution externe
      if(flag_biosponge_size==0) then !m°v°m> !24-07-18
       do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
       j1=j0 ; j2=jmax+1-j0
       do k=1,kmax
        do i=1,imax
         bio_relax_north(i,j0,k,vb,time_)=anyvar3d(i,j2,k)
         bio_relax_south(i,j0,k,vb,time_)=anyvar3d(i,j1,k)
        enddo
       enddo
       enddo
       do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
       i1=i0 ; i2=imax+1-i0
       do k=1,kmax
        do j=1,jmax

         bio_relax_east(j,i0,k,vb,time_)=anyvar3d(i2,j,k)
         bio_relax_west(j,i0,k,vb,time_)=anyvar3d(i1,j,k)

        enddo
       enddo
       enddo
      endif                           !m°v°m> !24-07-18
      if(flag_biosponge_size==1) then !mOvOm> !24-07-18
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
         bio_relax_full(i,j,k,vb,time_) &
              =anyvar3d(i,j,k)
       enddo       ; enddo       ; enddo
      endif                           !mOvOm> !24-07-18

      if(biobc_init_flag==1) then !--initialisation-->
       x1=rap_biobc
       if(time_==1)x1=1.
       do k=1,kmax
       do j=0,jmax+1
       do i=0,imax+1
        bio_t(i,j,k,vb)=x1*anyvar3d(i,j,k)+(1.-x1)*bio_t(i,j,k,vb)
       enddo
       enddo
       enddo
      endif                       !--initialisation-->

      end subroutine biobc_interp_step3

!..........................................................

      subroutine biobc_landvalue_step1(filval_)
      implicit none
      real*4 filval_
#ifdef synopsis
       subroutinetitle='biobc_landvalue_step1'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      k=biobc_kmax
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !----->

        i1=int(ij2biobc_i(i,j))-biobczoom_istr+1 !01-07-14
        j1=int(ij2biobc_j(i,j))-biobczoom_jstr+1

        rapi=ij2biobc_i(i,j)-int(ij2biobc_i(i,j))
        rapj=ij2biobc_j(i,j)-int(ij2biobc_j(i,j))

      if(biobc_var(i1  ,j1  ,k)==filval_)biobc_var(i1  ,j1  ,k)=-filval_
      if(biobc_var(i1  ,j1+1,k)==filval_)biobc_var(i1  ,j1+1,k)=-filval_
      if(biobc_var(i1+1,j1+1,k)==filval_)biobc_var(i1+1,j1+1,k)=-filval_
      if(biobc_var(i1+1,j1  ,k)==filval_)biobc_var(i1+1,j1  ,k)=-filval_

      endif                          !----->
      enddo
      enddo

      write(texte30,'(a,i0)')trim(tmpdirname)//'bouchetrou_biobc',par%rank
      open(unit=4,file=trim(texte30))

      k=biobc_kmax

      do 1862 j=1,biobc_jmax
      do 1862 i=1,biobc_imax

       if(biobc_var(i,j,k)==-filval_)then !§§§§§§§§§§§§§§§>

        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=max(j-i10,1)
         j2=min(j+i10,biobc_jmax)
         j3=k1*(j2-j0)+(1-k1)

         i0=max(i-i10+k1,1)
         i2=min(i+i10-k1,biobc_imax)
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3

         if(i1<1         )stop 'faute biobc fp5'  !01-07-10
         if(i1>biobc_imax)stop 'faute biobc fp6'
         if(j1<1         )stop 'faute biobc fp7'
         if(j1>biobc_jmax)stop 'faute biobc fp8'
! si le modele plante sur fautes 5 6 7 ou 8 c'est qu'il n'y a pas de bouchage
! suffisament proche du trou et que l'algo cherche au delà de la zone d'extraction,
! ce qu'on ne permet pas pour avoir une parallelisation parfaite.
! Ce qu'on peu faire: on peut augmenter la taille de la
! zone d'extraction en augmentant i0 (repere 1498). Mais cela revele surtout
! que les masques oceano et biobc sont trop differents et qu'il faut sans doute
! ameliorer le masque oceano

          if(     biobc_var(i1,j1,k)/=-filval_                      &
             .and.biobc_var(i1,j1,k)/= filval_)then !%%%%%%%%%%%%%>
           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then                    !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              dist1=dist2
             endif                                      !>>>>>>>>>>>>>
          endif                                      !%%%%%%%%%%%%%>
 1861    continue

 1863   continue
        i10=i10+1
        if(ksecu.eq.0)goto 1864
                 write(4,'(4(i4,1x))')i,j  ,i4,j4

       endif                           !§§§§§§§§§§§§§§§>
 1862 continue

      close(4)

      biobc_landvalue_flag=1


      end subroutine biobc_landvalue_step1

!..........................................................

      subroutine biobc_landvalue_step2
      implicit none
#ifdef synopsis
       subroutinetitle='biobc_landvalue_step2'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     open(unit=3,file=trim(tmpdirname)//'bouchetrou_biobc'//dom_c//'.out')
      write(texte30,'(a,i0)')trim(tmpdirname)//'bouchetrou_biobc',par%rank
      open(unit=3,file=trim(texte30))

      do j1=1,biobc_jmax
      do i1=1,biobc_imax
       read(3,*,end=1760)i,j,i4,j4
       do k=1,biobc_kmax
        biobc_var(i,j,k)=biobc_var(i4,j4,k)
       enddo
      enddo
      enddo
 1760 continue

      close(3)

      end subroutine biobc_landvalue_step2

!...........................................................

      subroutine biobc_allocate_relax
      implicit none
#ifdef synopsis
       subroutinetitle='biobc_allocate_relax'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(bio_relax_size==0)return
      if(vbmax==0)return

      if(flag_biosponge_size==0) then !24-07-18 ! eponges economes >
       allocate(bio_relax_north(imax,bio_relax_size,kmax,vbmax,2))
       allocate(bio_relax_south(imax,bio_relax_size,kmax,vbmax,2))
       allocate(bio_relax_east (jmax,bio_relax_size,kmax,vbmax,2))
       allocate(bio_relax_west (jmax,bio_relax_size,kmax,vbmax,2))
      endif                           !24-07-18 ! eponges economes >

      if(flag_biosponge_size==1) then !24-07-18 ! eponges larges >
       allocate(bio_relax_full(imax,jmax,kmax,vbmax,2))
      endif                           !24-07-18 ! eponges larges >

      end subroutine biobc_allocate_relax

!...........................................................

      subroutine biobc_nudging_area
      implicit none
      integer istr_,iend_,jstr_,jend_

      if(bio_relax_size==0)return

      if(flag_biosponge_size==1) then !pmx> !24-07-18> 
       call biobc_nudging_area_full
       return
      endif                           !pmx> !24-07-18> 

! Tableaux de melange de la solution avec la solution externe:
!      relax_bio=86400. ! temps de rappel = un jour sur la frontiere
       relax_bio=dti_fw ! Rappel 100% sur la frontiere
       x0=min(dti_fw/max(real(relax_bio,kind=8),small1),1.d0)
       x1=(1.-rap_biobc)
       x2=    rap_biobc

! Si frontière ouest (i=1) ouverte:
       if(obcstatus(ouest) == 1) then !--------->

        do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
        i=i0 !01-05-13
          do j=1,jmax
           xy_t(i,j,0)= 1.-sponge_t(i,j,1)*x0
           xy_t(i,j,1)=    sponge_t(i,j,1)*x0*x1
           xy_t(i,j,2)=    sponge_t(i,j,1)*x0*x2
          enddo
        enddo

        do vb=1,vbmax
        if(biobc_type(vb)==2) then !*************>
         do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
         i=i0 !01-05-13
         do k=1,kmax
          do j=1,jmax
            bio_t(i,j,k,vb)=xy_t(i,j,0)*bio_t(i,j,k,vb)             &
                           +xy_t(i,j,1)*bio_relax_west(j,i0,k,vb,1) &
                           +xy_t(i,j,2)*bio_relax_west(j,i0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>
        enddo

       endif                                        !--------->

! Si frontière est (i=imax) ouverte:
       if(obcstatus(est  ) == 1) then !--------->

        do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
         i=imax+1-i0 !01-05-13
          do j=1,jmax
           xy_t(i,j,0)= 1.-sponge_t(i,j,1)*x0
           xy_t(i,j,1)=    sponge_t(i,j,1)*x0*x1
           xy_t(i,j,2)=    sponge_t(i,j,1)*x0*x2
          enddo
        enddo

        do vb=1,vbmax
        if(biobc_type(vb)==2) then !*************>
         do i0=1,bio_relax_size ! i0 indice de largeur de la zone de rappel
         i=imax+1-i0 !01-05-13
         do k=1,kmax
          do j=1,jmax
            bio_t(i,j,k,vb)=xy_t(i,j,0)*bio_t(i,j,k,vb)             &
                           +xy_t(i,j,1)*bio_relax_east(j,i0,k,vb,1) &
                           +xy_t(i,j,2)*bio_relax_east(j,i0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>
        enddo

       endif                                        !--------->

! Si frontière sud (j=1) ouverte:
       if(obcstatus(sud  ) == 1) then !--------->
        istr_=1 ; iend_=imax
        if(obcstatus(ouest)==1)istr_=bio_relax_size+1
        if(obcstatus(est  )==1)iend_=imax-bio_relax_size


        do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
        j=j0
!         do i=bio_relax_size+1,imax-bio_relax_size !i=1,imax!Ne pas appliquer 2 fois le rappel!01-05-13
          do i=istr_,iend_                                   !Ne pas appliquer 2 fois le rappel!01-05-13
           xy_t(i,j,0)= 1.-sponge_t(i,j,1)*x0
           xy_t(i,j,1)=    sponge_t(i,j,1)*x0*x1
           xy_t(i,j,2)=    sponge_t(i,j,1)*x0*x2
          enddo
        enddo

        do vb=1,vbmax

        if(biobc_type(vb)==2) then !*************>
         do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
         j=j0
         do k=1,kmax
!         do i=bio_relax_size+1,imax-bio_relax_size !i=1,imax!Ne pas appliquer 2 fois le rappel!01-05-13
          do i=istr_,iend_                                   !Ne pas appliquer 2 fois le rappel!01-05-13
            bio_t(i,j,k,vb)=xy_t(i,j,0)*bio_t(i,j,k,vb)              &
                           +xy_t(i,j,1)*bio_relax_south(i,j0,k,vb,1) &
                           +xy_t(i,j,2)*bio_relax_south(i,j0,k,vb,2)

          enddo
         enddo
         enddo

        endif                      !*************>
        enddo

       endif                                        !--------->

! Si frontière nord (j=jmax) ouverte:
       if(obcstatus(nord ) == 1) then !--------->
        istr_=1 ; iend_=imax
        if(obcstatus(ouest)==1)istr_=bio_relax_size+1
        if(obcstatus(est  )==1)iend_=imax-bio_relax_size

        do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
        j=jmax+1-j0
!         do i=bio_relax_size+1,imax-bio_relax_size !i=1,imax!Ne pas appliquer 2 fois le rappel!01-05-13
          do i=istr_,iend_                                   !Ne pas appliquer 2 fois le rappel!01-05-13
           xy_t(i,j,0)= 1.-sponge_t(i,j,1)*x0
           xy_t(i,j,1)=    sponge_t(i,j,1)*x0*x1
           xy_t(i,j,2)=    sponge_t(i,j,1)*x0*x2
          enddo
        enddo

        do vb=1,vbmax
        if(biobc_type(vb)==2) then !*************>
         do j0=1,bio_relax_size ! j0 indice de largeur de la zone de rappel
         j=jmax+1-j0
         do k=1,kmax
!         do i=bio_relax_size+1,imax-bio_relax_size !i=1,imax!Ne pas appliquer 2 fois le rappel!01-05-13
          do i=istr_,iend_                                   !Ne pas appliquer 2 fois le rappel!01-05-13
            bio_t(i,j,k,vb)=xy_t(i,j,0)*bio_t(i,j,k,vb)              &
                           +xy_t(i,j,1)*bio_relax_north(i,j0,k,vb,1) &
                           +xy_t(i,j,2)*bio_relax_north(i,j0,k,vb,2)
          enddo
         enddo
         enddo
        endif                      !*************>
        enddo

       endif                                        !--------->

      end subroutine biobc_nudging_area

!..........................................................

      subroutine biobc_nudging_area_full
      implicit none
      if(bio_relax_size==0)return

! Tableaux de melange de la solution avec la solution externe:
       relax_bio=dti_fw ! Rappel 100% sur la frontiere
       x0=min(dti_fw/max(real(relax_bio,kind=8),small1),1.d0)
       x1=(1.-rap_biobc)
       x2=    rap_biobc

        do j=1,jmax ; do i=1,imax
           xy_t(i,j,0)= 1.-sponge_t(i,j,1)*x0
           xy_t(i,j,1)=    sponge_t(i,j,1)*x0*x1
           xy_t(i,j,2)=    sponge_t(i,j,1)*x0*x2
        enddo       ; enddo

        do vb=1,vbmax
        if(biobc_type(vb)==2) then !*************>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
            bio_t(i,j,k,vb)=xy_t(i,j,0)*bio_t(i,j,k,vb)            &
                           +xy_t(i,j,1)*bio_relax_full(i,j,k,vb,1) &
                           +xy_t(i,j,2)*bio_relax_full(i,j,k,vb,2)
         enddo       ; enddo       ; enddo
        endif                      !*************>
        enddo

      end subroutine biobc_nudging_area_full

!..........................................................
      end module module_biobc
