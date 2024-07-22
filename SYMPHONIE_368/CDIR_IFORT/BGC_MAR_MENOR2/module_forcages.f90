










      module module_forcages
!______________________________________________________________________
! S model
! release S26 - last update: 12-02-13
!______________________________________________________________________
! Version date     Description des modifications
!         29/12/08 Mise en service
!         13-05-09 bienvenue à TIDEVAR real*4
! 2009.3  01-10-09 ajout MULTI2DVAR
!                  ajout VAR_LONMIN,VAR_LATMIN,VAR_LONMAX,VAR_LATMAX
!                  ,VAR_SCALEFACTOR,VAR_MISVAL
!         02-11-09 symp_var1d & symp_z1d deviennent ogcm_var1d et ogcm_z1d
! 2010.2  16-12-09 des variables pour l'interpolation en ligne de la maree
! 2010.6  03-02-10 ajout meteo_lon meteo_lat
!         04-02-10 ajout meteo_kmax,meteozoom_istr,meteozoom_iend
!                  ,meteozoom_jstr,meteozoom_jend
!         10-02-10 ajout meteo_resol
! 2010.9  06-06-10 ajout d'un cas ww3 netcdf
! 2010.10 15-06-10 ajout meteonemo meteo_resol_u v
! 2010.11 15-07-10 tout ce qui n'est pas allocatable desormais declare
!                  dans module_principal
! 2010.12 27-08-10 ajout d'une variable de decumul des cumuls d'ecmwf
!         14-09-10 ajout ww3_var2
! 2010.19 13-04-11 short2d et short3d pour champs PSY4V1R3
! 2010.25 23-02-12 transferts ww3 vers module_wave
!         31-05-12 ogcm_lon ogcm_lat declares en double precision
! S26     12-02-12 la routine allocate_forcage entre dans le module
!                  Ajout meteo_short
!______________________________________________________________________
      implicit none
      double precision,dimension(:,:),allocatable :: &
        symp_lat                                     &
       ,symp_lon                                     &
       ,ogcm_lat                                     & !31-05-12
       ,ogcm_lon                                     & !31-05-12
       ,multi2d_lon                                  &
       ,multi2d_lat

       real*4,dimension(:,:),allocatable ::                             &      !04-06-10
        meteonemo_lon                                                   &
       ,meteonemo_lat
       real*4,dimension(:,:,:),allocatable ::                           &
        meteonemo_data

      real*4,dimension(:,:,:),allocatable ::                            &
        ogcm_var3d                                                      &
       ,symp_var3d                                                      &
       ,symp_z                                                          &
       ,tidevar                                                         &
       ,multi2d_var                                                     &
       ,ogcm_z

      real*4,dimension(:,:),allocatable ::                              &
        meteo_var                                                       &
       ,meteo_cum                                                       &
       ,ogcm_var2d                                                      &
       ,symp_var2d                                                      &
       ,ssh_var

      real*4,dimension(:,:),allocatable :: &
        dust_var  &
       ,dust_cum

      real*4,dimension(:,:),allocatable :: &
        depnit_var  &
       ,depnit_cum

      real*4,dimension(:,:),allocatable :: &
        depammo_var  &
       ,depammo_cum

      real*4,dimension(:),allocatable ::                                &
        ssh_lat                                                         &
       ,ssh_lon                                                         &
       ,ogcm_var1d                                                      &
       ,ogcm_z1d

      double precision,dimension(:,:),allocatable ::                      &
        meteo_lon                                                       &
       ,meteo_lat

      double precision,dimension(:,:),allocatable ::   &
        dust_lon                                      &
       ,dust_lat

      double precision,dimension(:,:),allocatable ::   &
        depnit_lon                                      &
       ,depnit_lat

      double precision,dimension(:,:),allocatable ::   &
        depammo_lon                                      &
       ,depammo_lat

      integer,dimension(:,:),allocatable ::            &
        short2d                                        &
       ,meteo_short                                    &
       ,dust_short                                     &
       ,depnit_short                                   &
       ,depammo_short

      integer,dimension(:,:,:),allocatable ::          &
        short3d

contains

      subroutine allocate_forcages(ki1,ki2,dim1_,dim2_,dim3_)
      implicit none
      integer ki1,ki2,dim1_,dim2_,dim3_

! ----------------------------------------------------------
! Forcages par un model_ oceanique de grande echelle:
      if(ki2.eq.1) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(ogcm_var3d    (dim1_,dim2_,dim3_))
      allocate(short3d       (dim1_,dim2_,dim3_))
      allocate(ogcm_var2d    (dim1_,dim2_         ))
      allocate(ogcm_lat      (dim1_,dim2_         ))
      allocate(ogcm_lon      (dim1_,dim2_         ))
      allocate(short2d       (dim1_,dim2_         ))
      allocate(ogcm_z        (dim1_,dim2_,dim3_))
      allocate(ogcm_var1d    (                  dim3_))
      allocate(ogcm_z1d      (                  dim3_))

      endif                !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(ogcm_var3d)
      deallocate(ogcm_var2d)
      deallocate(ogcm_lat)
      deallocate(ogcm_lon)
      deallocate(ogcm_z  )
      deallocate(ogcm_var1d)
      deallocate(ogcm_z1d)
      deallocate(short2d)
      deallocate(short3d)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ oceanique de grande echelle:
      return
      endif
! ----------------------------------------------------------


! ----------------------------------------------------------
! Forcages par un model_ meteo
      if(ki2.eq.2) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(meteo_var     (dim1_,dim2_))
      allocate(meteo_cum     (dim1_,dim2_))
      allocate(meteo_lon     (dim1_,dim2_))
      allocate(meteo_lat     (dim1_,dim2_))
      allocate(meteo_short   (dim1_,dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(meteo_var)
      deallocate(meteo_cum)
      deallocate(meteo_lon)
      deallocate(meteo_lat)
      deallocate(meteo_short)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ meteo
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par un model_ dust
      if(ki2.eq.11) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(dust_var     (dim1_,dim2_))
      allocate(dust_cum     (dim1_,dim2_))
      allocate(dust_lon     (dim1_,dim2_))
      allocate(dust_lat     (dim1_,dim2_))
      allocate(dust_short   (dim1_,dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(dust_var)
      deallocate(dust_cum)
      deallocate(dust_lon)
      deallocate(dust_lat)
      deallocate(dust_short)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ dust
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par un model_ depnit
      if(ki2.eq.12) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(depnit_var     (dim1_,dim2_))
      allocate(depnit_cum     (dim1_,dim2_))
      allocate(depnit_lon     (dim1_,dim2_))
      allocate(depnit_lat     (dim1_,dim2_))
      allocate(depnit_short   (dim1_,dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(depnit_var)
      deallocate(depnit_cum)
      deallocate(depnit_lon)
      deallocate(depnit_lat)
      deallocate(depnit_short)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ depnit
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par un model_ depammo
      if(ki2.eq.13) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(depammo_var     (dim1_,dim2_))
      allocate(depammo_cum     (dim1_,dim2_))
      allocate(depammo_lon     (dim1_,dim2_))
      allocate(depammo_lat     (dim1_,dim2_))
      allocate(depammo_short   (dim1_,dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(depammo_var)
      deallocate(depammo_cum)
      deallocate(depammo_lon)
      deallocate(depammo_lat)
      deallocate(depammo_short)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ depammo
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par aviso
      if(ki2.eq.3) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(ssh_var     (dim2_,dim1_))
      allocate(ssh_lon     (dim1_))
      allocate(ssh_lat     (dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(ssh_var)
      deallocate(ssh_lon)
      deallocate(ssh_lat)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par aviso
      return
      endif
! ----------------------------------------------------------



! ----------------------------------------------------------
! Forcages par un model_ de maree
      if(ki2.eq.6) then
! ----------------------------------------------------------
      stop ' erreur 6 dans allocate_forcage'

!     if(ki1.eq.1) then !111111111111111>

!     allocate(tidevar   (0:dim1_-1,0:dim2_-1,1))                 !09-10-09

!     endif             !111111111111111>

!     if(ki1.eq.2) then !222222222222222>

!     deallocate(tidevar)

!     endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ ogcm
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par un model_ 2D multivariables                             !01-10-09
      if(ki2.eq.7) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(multi2d_var (1:dim1_,1:dim2_,1:dim3_))
      allocate(multi2d_lon (1:dim1_,1:dim2_))
      allocate(multi2d_lat (1:dim1_,1:dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>

      deallocate(multi2d_var)
      deallocate(multi2d_lon)
      deallocate(multi2d_lat)

      endif             !222222222222222>

! ----------------------------------------------------------
! Forcages par un model_ 2D multivariables
      return
      endif
! ----------------------------------------------------------

! ----------------------------------------------------------
! Forcages par meteo nemo                              !15-06-10
      if(ki2.eq.9) then
! ----------------------------------------------------------

      if(ki1.eq.1) then !111111111111111>

      allocate(meteonemo_data (1:dim1_,1:dim2_,1:dim3_))
      allocate(meteonemo_lon (1:dim1_,1:dim2_))
      allocate(meteonemo_lat (1:dim1_,1:dim2_))

      endif             !111111111111111>

      if(ki1.eq.2) then !222222222222222>
      deallocate(meteonemo_data)
      deallocate(meteonemo_lon)
      deallocate(meteonemo_lat)

      endif             !222222222222222>


! ----------------------------------------------------------
      return
      endif
! ----------------------------------------------------------


      end subroutine allocate_forcages



      end module module_forcages

