      subroutine river_bio_upd
!______________________________________________________________________
!
! S model
! release 2010.20  - last update: 18-04-11
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http:// 
!
!______________________________________________________________________

      use ModuleComBlock
      use UserDeclaration
      use ModuleDeclaration
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      real ox_A0,ox_A1,ox_A2,ox_A3,ox_A4,ox_A5,ox_B0,ox_B1,ox_B2,  &
           ox_B3,ox_A,ox_B,ox_C,ox_D,ox_C0,Ts,A                            
      integer runofftimecounter_
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
 
! modifs: 22/05/06: correction bug sur K1 (KR au lieu de KR)
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes
! 31-05-2013 Oxygene a saturation à partir de la temperature avec salinité nulle
! 07-07-2015 deplacement call river_upd hors boucle sur rivieres

!----------------------------------------------------------------------
!     BUT:
!     ----
!     Ce programme lit les concentrations des fleuves dan un fichier
!     et met a jour à chaque pas de temps
!
!
!     OUTPUT:
!     ------
!     RIVER_BIO
!
!
!     INPUT:
!     ------
!     ILD_TO_SD  ! liquid discharge to solid discharge
!
!
!     MODIFS:
!     -------
!
!-----------------------------------------------------------------------
! mise en service: 05/04/06
!    mises a jour: 06/04/06 argument ICHOIX inutile
!-----------------------------------------------------------------------


      if(nriver.eq.0)return ! Pas de riviere?  alors on quitte la routine.
      if(vbmax.eq.0) return ! Pas de biologie? alors on quitte la routine.


! LE CAS PARTICULIER DE L'ETAT INITIAL
! 1. Definir les parametres d'avancee dans le temps (date repere & echantillonage)
! Debut:
      if(iteration3d==kount0) then !.................!       !18-04-11

        do kr=1,nriver ! debut boucle sur les rivieres
         i1=bioriv_date(1,kr)
         i2=bioriv_date(2,kr)
         i3=bioriv_date(3,kr)
         i4=bioriv_date(4,kr)
         i5=bioriv_date(5,kr)
         i6=bioriv_date(6,kr)
!         write(6,*)'appel à date_to_kount dans river_bio_upd'
!        call date_to_kount
         call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!        bioriv_dt(kr,2)=xdtk_out                     ! date repere NC=1 en kount
!        bioriv_dt(kr,1)=bioriv_info(kr)*3600./dti_fw ! duree, en kount, entre 2 infos
         bioriv_dt(kr,2)=elapsedtime_out              !18-04-11
         bioriv_dt(kr,1)=bioriv_info(kr)*3600.        !18-04-11

        enddo          ! fin boucle sur les rivieres

      endif                    !.................!
! Fin.

! modif OXYGENE ici
      ox_A0= 2.00907
      ox_A1= 3.22014
      ox_A2= 4.05010
      ox_A3= 4.94457
      ox_A4= -0.256847
      ox_A5= 3.88767

      ox_B0= -6.24523e-3
      ox_B1= -7.37614e-3
      ox_B2= -1.03410e-2
      ox_B3= -8.17083e-3

      ox_C0= -4.88682e-7

! coeff pour le calcul du nombre de schmidt: from Keeling et al 1998
      ox_A=1638
      ox_B=81.83
      ox_C=1.483
      ox_D=0.008004


      call river_upd(2) ! pour le calcul de river_t  ! claude 07-07-15
      do 1000 kr=1,nriver ! boucle sur les rivieres

! modif OXYGENE ici
     
!     call river_upd(2) ! pour le calcul de river_t

      Ts = log((298.5-river_t(kr,1))/(273.15+river_t(kr,1)))

      A= (ox_A0 + ox_A1*Ts + ox_A2 * Ts**2 + ox_A3 *Ts**3     &
           + ox_A4 *Ts**4 + ox_A5 * Ts**5)

      river_bio(ioxygen,kr,1) = (1000./22.3916) * exp(A)


      if(ild_to_sd(kr).eq.1) then ! ILD ILD ILD ILD >                  !22/05/06


      ksecu=0 ! A priori on n'ouvre pas le fichier de "donnee riviere"

! Est on de part et d'autre d'une nouvelle echeance?
      x2=(elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)     !18-04-11
      x1=(elapsedtime_bef-bioriv_dt(kr,2))/bioriv_dt(kr,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
! Si la reponse est oui alors il faut lire le fichier:
      if(j2-j1.eq.1)         ksecu=1  ! indique que nous lirons le fichier
      if(iteration3d==kount0)ksecu=1  ! L'initialisation impose de lire le fichier !18-04-11


      if(ksecu.eq.1) then !*************************************>

      if(iteration3d==kount0)then !---------------->        !18-04-11
          k4=0 ! à l'initialisation on lit 2 champs
      else                        !---------------->                       !27/04/04
! sinon l'ancienne echeance est sauvée avant de lire la suivante:
          k4=2
          do vb=1,vbmax
           river_bio(vb,kr,0)=river_bio(vb,kr,2)
          enddo
      endif                       !---------------->



! Ouverture du fichier:
      open(unit=3,file=bioriv_file(kr)                                  &
                 ,access='direct'                                       &
                 ,recl=4*vbmax                                          &
!                ,recl=4*kmax                                          &
                 ,form='unformatted')

!        print*,'rec river',


! lecture du fichier:
      do k3=k4,2,2
         x2=(elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)      !18-04-11
         nc=int(1.+x2)+k3/2

!         print*,'rec river',nc,vbmax,elapsedtime_now,bioriv_dt(kr,2),bioriv_dt(kr,1)
!         print*,'file',bioriv_file(kr),kr,k3,4*vbmax,nc

         read(3,rec=nc)(river_bio(vb,kr,k3),vb=1,vbmax)
      enddo

      close(3)
      endif               !*************************************>



! Interpolation temporelle entre 2 echeances du fichier:
      rap=      (elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1)  & !18-04-11
          -int( (elapsedtime_now-bioriv_dt(kr,2))/bioriv_dt(kr,1) )

      do vb=1,vbmax
       river_bio(vb,kr,1)=(1.-rap)*river_bio(vb,kr,0)                  &
                         +    rap *river_bio(vb,kr,2)
      enddo


      endif                       ! ILD ILD ILD ILD >
 1000 continue ! fin de boucle sur les rivieres

      if(flag_nemoffline==1) then !wwwwwwwwwwwwwwwww>
! lecture des runoff
! suivant le mois, lecture des runoffs: couverture spatiale correspondant au proc
! cas initial ? lecture deux fichiers suivant jour de l'année?? en utilisant une formule comme ci dessus
!avec un espacement de 30 jours
! si oui interpolation des runoffs 
! puis peut etre ajouter une subroutine appelée dans module_biobalance_gateC comme oxygen_surface_flux
! pour l'oxygene ça s'ajoute au flux air mer?? qui est dans fluxbio
      ksecu=0 ! A priori on n'ouvre pas le fichier de "donnee riviere"
      if(iteration3d==kount0) then  !!!!
!definition de runoff_dt
      call datetokount(year_now-1,1,15,0,0,0)   
      runoff_dt(2)=elapsedtime_out          
      runoff_dt(1)=30.4375*86400.    ! un mois moyen
      endif                         !!!!

! Est on de part et d'autre d'une nouvelle echeance?
      x2=(elapsedtime_now-runoff_dt(2))/runoff_dt(1)    
      x1=(elapsedtime_bef-runoff_dt(2))/runoff_dt(1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
! Si la reponse est oui alors il faut lire le fichier:
      if(j2-j1.eq.1)         ksecu=1  ! indique que nous lirons le fichier
      if(iteration3d==kount0)ksecu=1  ! L'initialisation impose de lire le fichier

      if(ksecu.eq.1) then !*************************************>   il faut lire

      if(iteration3d==kount0)then !---------------->       
          k4=0 ! à l'initialisation on lit 2 champs
      else                        !---------------->     
! sinon l'ancienne echeance est sauvée avant de lire la suivante:
          k4=2
          do kr=1,nriver ! boucle sur les rivieres
           runoff_w(kr,0)=runoff_w(kr,2)
          enddo
      endif                       !---------------->

! lecture fichier : 
! condition dans le proc : rivertrc_inout(kr)=1

!     filename_= &
!      '/sirocco/mex/SIMU_SIMED/METEO/rnf_clim_ruiss_par_bassin_jo.nc'
      status=nfmpi_open(par%comm2d,filename_runoff,nf_nowrite, MPI_INFO_NULL,ncid1)
      status=nfmpi_inq_varid(ncid1,'sorunoff',var_id)
      do k3=k4,2,2

         x2=(elapsedtime_now-runoff_dt(2))/runoff_dt(1)    
         runofftimecounter_=int(1.+x2)+k3/2
         runofftimecounter_=mod(runofftimecounter_-1,12)+1

      start(3)=runofftimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                      ,anyvar2d(0:imax+1,0:jmax+1))

       do kr=1,nriver
        if(rivertrc_inout(kr)==1)then
         runoff_w(kr,k3)=anyvar2d(iriver(kr,1),jriver(kr,1))
!      if(par%rank==0.and.rivertrc_inout(kr)==1)write(20,*)year_now,month_now,day_now,k3,runofftimecounter_, &
!         kr,iriver(kr,1)+par%timax(1),jriver(kr,1)+par%tjmax(1),runoff_w(kr,k3)
        endif
       enddo

      enddo ! k3
      status=nfmpi_close(ncid1)

      endif         !*************************************>      ksecu  : fin lecture

! Interpolation temporelle entre 2 echeances du runoff des rivieres:
      rap=      (elapsedtime_now-runoff_dt(2))/runoff_dt(1)  & !18-04-11
          -int( (elapsedtime_now-runoff_dt(2))/runoff_dt(1) )

      do kr=1,nriver ! boucle sur les rivieres
       if(rivertrc_inout(kr)==1) &
       runoff_w(kr,1)=((1.-rap)*runoff_w(kr,0)                  &
                   +    rap *runoff_w(kr,2))                  &
              /1000.         ! passer de kg/m2/s en m3/m2/s
!      if(kr==1.and.rivertrc_inout(kr)==1)write(21,*)month_now,day_now,runoff_w(kr,0),runoff_w(kr,1),runoff_w(kr,2)
      enddo


      endif !wwwwwwwwwwwwwwwww>  flag_nemoffline


      return
      end
