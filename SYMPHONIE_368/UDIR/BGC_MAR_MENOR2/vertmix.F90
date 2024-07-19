      subroutine vertmix_tem(flagsolver_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 345 - last update: 06-04-22
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_cpl_oasis
      implicit none
      integer flagsolver_
#ifdef synopsis
       subroutinetitle='vertmix_tem'
       subroutinedescription=' Updates temperature field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Computes turbulent mixing and surface forcing terms for temperature

!...............................................................................
! Version date      Description des modifications
!         18/07/01: passage à la coordonnée sigma généralisée
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         11/12/01: pour que les coef de la matrice principale soient les
!                   mêmes pour toutes les variables _Z 1/2 niveaux on deplace
!                   l'eponge dans le partie explicite
!         20/06/02: le commentaire precedent est absurde: la matrice principale
!                   est la même pour T et S qu'il y ait ou non le terme
!                   implicite de relaxation. Il est donc remis en place.
!         28/02/03: amenagements pour coordonnee verticale hybride
!         28/12/05: Remise en service de l'albedo et d'un calcul de l'absorption
!                   sous reserve que IALBEDO=1
!         02/02/06: nouvelle penetration de la lumiere: combinaison de 2
!                   exponentielles
!         16/02/06: finalisation du point precedent avec des valeurs à lire
!                   dans un nouveau notebook: notebook_optical
!         30/03/06: appel (au choix) a TRIDIAGSYST_DBP
!         10/12/07: Condition limite surface. Les flux ne sont plus appliqués
!                   sur le niveau de surface mais répartis dur la couche de
!                   mélange de surface dont la profondeur est donnee par KSL_Z
!         31/01/08: Eponges transferees dans routine apply_sponge_layer
!         11/11/08: Division par facteur d'echelle verticale pris au temps t+1
! 2009.3  01-10-09  l'appel au calcul de l'albedo passe dans model_3d
!                   et model_offline
!         02-10-09  dti_lp remplace dti
!         03-10-09  utilisation des nouveaux facteurs d'echelle verticale
!                   on note que xy_z(i,j,3) est prealablement divisé par
!                   hz_z pour preserver l'homogeneité du nouveau calcul à
!                   suivre
!         08-10-09  albedo devient albedo_z en 2d
! 2010.3  12-01-10  bancs decouvrants
!         15-01-10  subroutine temperature renommée vertmix_tem
! 2010.7  14-02-10  Rationnalisation des calculs
! 2010.8  19-03-10  - Ajout heatrelax
!                   - kmax remplace kmax, k+1 remplace k+1, km1 remplace k-1
! 2010.9  06-06-10  tridiagsyst_dbp renommé tridiagonalsolver
! 2010.15 27-12-10  - Prise en compte de l'advection de la temperature
!                   en meme temps que le melange vertical. Modification
!                   du dernier argument de anyv3d pour harmonisation
!                   avec l'etape d'advection.
!                   - filtre temporel ordre élevé
! 2010.16 12-01-11 merge des modifs sss et filtre temporel élevé
! 2010.25 08-02-12 kz_w renomme kh_w
! S.26    26-06-14 modif arg 1 de routine tridiagonalsolver(0,....)
!         15-08-14 coef de la matrice de melange calcules dans routine vertmix_matrix
!         12-10-14 utiliser lv plutot qu'une valeur numerique en dur
!         12-12-14 c.l. sous le fond
!         09-04-15 Ajout d'un commentaire
!         11-07-15 - multiplier T,S par mask_t (pour eviter que certaines options
!                  de compilation ne bloquent inutilement sur des valeurs hors
!                  gamme du masque continental).
!                  - melange additionel echelle moleculaire
!         05-08-15 methode vst continue
!         14-09-15 ajout subroutine vertmix_nrj
!         05-01-16 advection pas de temps separEs
!         12-01-16 zones decouvertes: derniere valeur valide prise au temps 1
!         28-01-16 nouveau filtre 2dt-Ndt
!         09-02-16 filtre FB et LP: nouveaux noms pour les coef tfilterfb et tfilterlf
!         10-03-16 Seuillage de l'argument de la fonction exp
!         27-03-16 slhf_w*const1+precipi remplacE par son equivalent omega_w(kmax+1)
!         24-04-16 - nouveau calcul pour la specific heat capacity
!                  - rho remplacE par rho+rhp dans le calcul des flux de surface
!         10-08-16 - ajout des subroutines vertmix_coarser2_xxxx dont l'appel est pour
!                  le moment commentE
!                  - neutralisation par etiquette "bidon" de la C.L. au fond de la
!                  coordonnee VST
!         29-09-16 couplage
!         25-04-17 coordonee z-s
!         17-10-17 suite du point precedent
!         07-12-17 tableaux kundermin reduise le calcul de la moyenne aux points
!                  effectifs
!         24-01-18 debug lignes sous clef de couplage
! v246    21-02-19 Modif du dimensionnement du coef de melange dans la couche fusionnEe
! v247    25-02-19 - Ajout d'un commentaire important
!                  - wetdrying shema 6: Ne pas appliquer flux evap/precip si wetmask_t/=1
!         28-02-19 omega_w remplacee par 0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))
!         01-03-19 ajout vertmix_merged_levels_ts(time_dz_)
!                     et vertmix_merged_levels_uv(time_dz_)
!         02-03-19 schema de flux d'eau A la surface coherente avec
!                  l'advection verticale implicite
! v249    07-03-19 ajout vertmix_merged_levels_bio et vertmix_merged_levels_2t_any
! v251    10-04-19 call obc_river_groundwater_s sources sous marine
! v252    17-04-19 subroutine obc_river_groundwater_s renommEe obc_river_groundwater_ts
! v253    01-05-19 evitement de la division par dz
! v256    28-05-19 cas particulier du modele 1DV et check5,check6
!         29-05-19 deplacement de la C.L. sur tridia_in(:,:,k=1,1)
! v260    19-10-19 correction bug tridia_in(i,j,kmax,3)=0.
! v269    04-12-19 couches fusionnees des vitesses: c'est l'advection que l'on moyenne
!                  et non plus le champs total
! v285    03-06-20 ajout lien https
! v287    22-08-20 if(oasis_symsym_onoff==1)call cpl_oasis_sal_delta , cpl_oasis_tem_delta
! v296    19-02-21 if(flag_surfriver==0)call obc_river_surfriver_tem  !19-02-21
! v297    05-03-21 rivieres de surface et sous marine integrees au schema d'advection
! v303    03-08-21 modif couche fusionnee
! v310    02-11-21 Ne pas amplifier KH(kmerged) si melange dans couche merged supprimE
!         04-11-21 modif de la routine vertmix_merged_levels_uv avec un
!                  melange progressif proportionnel au chevauchement de la couche
!                  de fond virtuelle et de la couche kmerged
! v325    12-02-22 call my_outputs_3dflux_s
! v327    19-02-22 Application d'un facteur correctif de la lumiere reflechie sur le fond 
! v332    28-02-22 +wetmask_wi_t(i,j)*5. melanger couche wetdry 
! v333    02-03-22 +(1-wetmask_wi_t(i,j))*5. melanger couche wetdry 
! v334    03-03-22 version unifiee de vertmix_matrix
! v337    18-03-22 tridia_in(:,:,k=1,1) renforcE de 2 A kmin_w
! v338    19-03-22 suite du point precedent
! v345    06-04-22 je supprime wetmask_t(i,j)* qui me semble inutile 
!...............................................................................
!    _________                    .__                  .__                     !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................


!$ standard (3D) case: temperature has to be updated from the intermediate
!$ "after" solution just provided by the advection subroutine:
!     if(flag3d.eq.1)kts=after  ! standard 3D case. after=2
!$ special 1D vertical case: advection is not computed and thus temperature
!$ has to be updated from the value of the previous time step
!     if(flag3d.eq.0)kts=before ! special 1D vertical case. before=0


! PENSER A INITIALISER DES VALEURS POUR SSR (surface solar radiation)
! ET POUR SNSF (surface non solar radiation)
! ET POUR SLHF (surface latent heat flux)
!     Par convention un flux descendant est positif (T augmente)
!     un flux negatif refroidit la surface
! Calcul de l'angle zenithal et de l'albedo:
! par defaut ALBEDO=0.
!     albedo=0.
!     if(ialbedo.eq.1)call albedo_upd

!     const1=dti_lp/(rho*cp)
!     const4=const1*(1.-albedo)
!     const4=const1*(1.-albedo_z(imax/2,jmax/2))

! evaporation calculee a partir du flux de chaleur latente.
! convention: si eva est positif (SLHF negatif) il y a evaporation
! et SAL_Z augmente

! Depuis le 09-04-15 les flux a l'interface air/mer modifient le bilan de masse, autrement
! dit la ssh monte et descend en fonction de l'evaporation/condensation/precipitation,
! grace a omega_w(i,j,kmax+1) dans module_external_mode.F90. Pour autant, nous n'avons
! pas changer la façon de calculer la variation de T et S afin de minimiser les modifications
! dans le modèle (et risque d'erreur), en particulier parce que l'algo de la couche de surface
! convective (sur laquelle se repartissent les flux) n'est pas simple. Pour que omega_w(kmax+1)
! n'impacte pas le bilan (puisqu'il continue de se faire comme avant) on utilise une condition
! de gradient nul au sommet de la couche. L'advection de la couche kmax se fait de surcroit
! par methode implicite (voir advection_scal.F90).

!     write(6,*)'albdo',albedo_w
!     stop 'rourou'

      const1=dti_lp/(rho*cp)
      do j=1,jmax ; do i=1,imax

! Formule pour Cp expliquee dans: !03-06-20
! https://docs.google.com/document/d/1xJmiLSxhQIl4Z_1RH2h3jjf5AHUw88E00-Go4KOdNbA/edit

      xy_t(i,j,iDtOvRhCp)=dti_lp/( & !pmx>                   ! 2.Dt /     !24-04-16
                                   (rhp_t(i,j,kmax)+rho)    &! Surface density
                        *(4190.-5.7*sal_t(i,j,kmax,before)) &! Cp (heat capacity)
                                 )   !pmx>


      xy_t(i,j,1)=                                                    &
        light_rat1*exp(max(-600.,( depth_w(i,j,1)-ssh_int_w(i,j,1)) &
                         *light_kpar1))                                &
       +light_rat2*exp(max(-600.,( depth_w(i,j,1)-ssh_int_w(i,j,1)) &
                         *light_kpar2_w(i,j)))

      enddo       ; enddo

      do j=1,jmax ; do i=1,imax

       xy_t(i,j,3)=xy_t(i,j,1)/hz_w(i,j,1)                               !03-10-09

! Application d'un facteur correctif de la lumiere reflechie sur le fond expliquE dans: !19-02-22
! https://docs.google.com/document/d/1xJmiLSxhQIl4Z_1RH2h3jjf5AHUw88E00-Go4KOdNbA/edit
       ssr_w(i,j,1)=ssr_w(i,j,1)*(1-xy_t(i,j,1)**2)

      enddo       ; enddo

      do j=1,jmax ; do i=1,imax
       xy_t(i,j,4)=xy_t(i,j,iDtOvRhCp)*(1.-albedo_w(i,j))*ssr_w(i,j,1)   !24-04-16
      enddo       ; enddo

!$ Standard case: the coef of the tridiagonal system are now computed:
      !if(flag3d==1) then !ooo> !28-05-19
      ! thom_comment je pense mieux ainsi
      if(flag_1dv/=1) then !ooo> !31-05-19
       if(anyv3d(imax/2,jmax/2,kmax-1,id_tdiv)/=check5)stop 'check5'
      else               !ooo> !28-05-19
! Cas particulier du modele 1DV
          anyv3d(:,:,:,id_tdiv)=0.
      endif              !ooo> !28-05-19


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Note: le filtre temporel est de la forme A(t+1)-A(t) avec A(t+1)=tfilterfb*B(t+1)+tfilterlf*C(t+1)
! B(t+1)=dz(t+1)*(T(t+1)-T(t))+dz(t)*(T(t)-T(t-1)) c.a.d. la moyenne de la derivee sur 2 temps successifs
! C(t+1)=dz(t+1)*(T(t)-T(t-1))
! De sorte que A(t+1)-A(t)=tfilterfb*[ dz(t+1)*(T(t+1)-T(t))-dz(t-1)*(T(t-1)-T(t-2)) ]
!                         +tfilterlf*[ dz(t+1)*(T(t)-T(t-1))-dz(t  )*(T(t-1)-T(t-2)) ]

! coef4 --->
      xy_t(i,j,2)=                                                         &
        light_rat1*exp(max(-600.,( depth_w(i,j,k+1)-ssh_int_w(i,j,1))*light_kpar1)) &
       +light_rat2*exp(max(-600.,( depth_w(i,j,k+1)-ssh_int_w(i,j,1))*light_kpar2_w(i,j)))

          tridia_in(i,j,k,4)=(                               &

! time filter:
!!!!!!       dz_t(i,j,k,before)*tem_t(i,j,k,before)          &
!            dz_t(i,j,k,before)*tem_t(i,j,k,2)               & !05-01-16

!           +dz_t(i,j,k,after)*(tfc1*tem_t(i,j,k,now)        & ! t      !27-12-10
!                              +tfc2*tem_t(i,j,k,before)     & ! t-dt
!                              +tfc3*tem_t(i,j,k,before2))   & ! t-2dt
!           -dz_t(i,j,k,now)  *(tfc0*tem_t(i,j,k,now)        & ! t
!                              +tfc1*tem_t(i,j,k,before)     & ! t-dt
!                              +tfc2*tem_t(i,j,k,before2)    & ! t-2dt
!                              +tfc3*tem_t(i,j,k,before3))   & ! t-3dt

       dz_t(i,j,k,before)*( tem_t(i,j,k,2)                    & !05-01-16

                        +tfilterfb*( tem_t(i,j,k,before2)         & !28-01-16
                                    -tem_t(i,j,k,before ) ))      &

      +dz_t(i,j,k,now  )*tfilterlf*( tem_t(i,j,k,before2)         &
                                    -tem_t(i,j,k,before ))        &

      +dz_t(i,j,k,after)*( (tfilterlf-tfilterfb)*tem_t(i,j,k,now   )  &
                                     -tfilterlf *tem_t(i,j,k,before)) &

!      -wetmask_t(i,j)*anyv3d(i,j,k,id_tdiv) &
! je supprime wetmask_t(i,j)* qui me semble inutile le !06-04-22
                      -anyv3d(i,j,k,id_tdiv) & 

! Advection:
     ! +( anyv3d(i+1,j  ,k  ,1)                              & !05-01-16
     !   -anyv3d(i  ,j  ,k  ,1)                              &
     !   +anyv3d(i  ,j+1,k  ,2)                              &
     !   -anyv3d(i  ,j  ,k  ,2) )/dxdy_t(i,j)                &
     !   +anyv3d(i  ,j  ,k+1,3)                              & !13-11-16
     !   -anyv3d(i  ,j  ,k  ,3)                              &

! Solar flux:
              +xy_t(i,j,4)*    ( xy_t(i,j,2)                 & !08-10-09
                                -xy_t(i,j,1)                 &
                                +xy_t(i,j,3)*dz_t(i,j,k,1) ) &

                             )                                          !12-08-14

      xy_t(i,j,1)=xy_t(i,j,2)

      enddo
      enddo
      enddo

! Retroaction du modele enfant sur le modele parent:
      if(oasis_symsym_onoff==1)call cpl_oasis_tem_delta !22-08-20

! Additional mixing if non-zero heat molecular diffusivity: !11-07-15
! Note: compte tenu de la valeur faible de kmol_h le calcul
! est fait en explicite:
      if(kmol_h>0.) then !kmkmkmkm>
       do k=1,kmax
        kp1=min(k+1,kmax)
        km1=max(k-1,1)
! Note a propos de km1: le fait que T soit homogene entre
! k=1 et k=kmin_w(i,j) garantit la condition de flux nul a
! travers le fond si kmin/=1
        do j=1,jmax
        do i=1,imax
         tridia_in(i,j,k,4)=                &
         tridia_in(i,j,k,4)+kmol_h*dti_lp*( & !oooooo>

             (  tem_t(i,j,kp1,bef)-tem_t(i,j,k,bef)) &
            /(depth_t(i,j,k+1)  -depth_t(i,j,k))     &

            -(  tem_t(i,j,k,bef)-tem_t(i,j,km1,bef)) &
            /(depth_t(i,j,k)  -depth_t(i,j,k-1))     &

                                          )   !oooooo>
        enddo
        enddo
       enddo
      endif              !kmkmkmkm>


!......................................................                !28/02/03
!$ Surface and bottom boundary conditions for coef4, coef1, coef3
      do j=1,jmax                                                      !28/02/03
      do i=1,imax                                                      !28/02/03

      x1=xy_t(i,j,iDtOvRhCp)                       &                   !24-04-16
               *(snsf_w(i,j,1)                     &
                +slhf_w(i,j,1)                     &
                +sshf_w(i,j,1)                     &
           +heatrelax_w(i,j,1))                    & !19-03-10
             /( depth_w(i,j,kmax+1)- depth_w(i,j,ksl_t(i,j)))

       do k=kmax,ksl_t(i,j),-1
       tridia_in(i,j,k,4)=                                     &
       tridia_in(i,j,k,4)+x1*(depth_w(i,j,k+1)-depth_w(i,j,k))
       enddo

      enddo
      enddo
!......................................................

!     call tridiagonalsolver(1,1,imax,1,jmax,kmax)         !30/03/06
      if(flagsolver_/=1)stop 'Unexpected flagsolver_ in T equation'
      call tridiagonalsolver(flagsolver_,1,imax,1,jmax,kmax)

        do k=1,kmax
        do j=1,jmax
        do i=1,imax

                   tem_t(i,j,k,2)=     &          !10-10-09
              tridia_out(i,j,k)        &
                 *mask_t(i,j,kmax)     &          !11-07-15
              *wetmask_t(i,j)          &          !12-01-10
          +(1.-wetmask_t(i,j))         &
                  *tem_t(i,j,k,1) !12-01-16

        enddo
        enddo
        enddo

! note: le cas des sources sous marines est traite dans l'etape salinite
! par la subroutine obc_river_groundwater_ts !17-04-19


      end subroutine vertmix_tem

!......................................................................

      subroutine vertmix_sal(flagsolver_)
      use module_principal ; use module_parallele ; use module_cpl_oasis
      use module_my_outputs
      implicit none
      integer flagsolver_
#ifdef synopsis
       subroutinetitle='vertmix_sal'
       subroutinedescription='Updates salinity field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!$ Computes turbulent mixing and surface forcing terms for salinity

!...............................................................................
! Version date      Description des modifications
!         18/07/01: passage à la coordonnée sigma généralisée
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         11/12/01: pour que les coef de la matrice principale soient les
!                   mêmes pour toutes les variables _Z 1/2 niveaux on deplace
!                   l'eponge dans le partie explicite
!         20/06/02: le commentaire precedent est absurde: la matrice principale
!                   est la même pour T et S qu'il y ait ou non le terme
!                   implicite de relaxation. Il est donc remis en place.
!         24/02/03: ajout de commentaires sur le flux de surface
!         28/02/03: arrangements pour grille hybride
!         19/03/03: Introduction des precipitation, suppression de EVA_Z
!                   car ce stockage est inutile et on economise ainsi un tableau.
!         26/03/06: modifications des commentaires en tête de routine
!         30/03/06: appel (au choix) à TRIDIAGSYST_DBP
!         10/12/07: Condition limite surface: les flux ne sont plus appliqués
!                   sur le niveau de surface mais répartis sur la couche de
!                   mélange de surface dont la profondeur est donnee par KSL_Z
!         31/01/08: Eponges transferees dans routine apply_sponge_layer
!         11/11/08: Division par facteur d'echelle verticalle pris au temps t+1
! 2009.3  02-10-09: dti_lp remplace dti
! 2010.3  12-01-10  bancs decouvrants
!         15-01-10  subroutine salinity renommée vertmix_sal
! 2010.7  14-02-10  Rationnalisation des calculs
! 2010.9  06-06-10  tridiagonalsolver renommé tridiagonalsolver
! 2010.15 27-12-10  - Prise en compte de l'advection de la salinite
!                   en meme temps que le melange vertical. Modification
!                   du dernier argument de anyv3d pour harmonisation
!                   avec l'etape d'advection.
!                   - Filtre temporel ordre élevé.
! 2010.16 12-01-11 merge des modifs sss et filtre temporel élevé
! 2010.25 08-02-12 kz_w renomme kh_w
! S.26    26-06-14 modif arg 1 de routine tridiagonalsolver(0,....)
!...............................................................................


!..................................................................... !24/02/03
!  La constante A=lv=2.47E6 (J/Kg) est la chaleur latente de condensation.
!  Le flux de chaleur latente SLHF (en W/m2 avec W=J/s) divisé par A donne
!  l'eau évaporée (ou condensée) correspondante en kg par seconde et par m2.
!  En 1 pas de temps (multiplication par DTI) on a donc la quantité d'eau en
!  kg par m2 évaporée ou condensée. Si on divise par la densité de l'eau (RHO)
!  on convertit cette info en volume d'eau évaporée ou condensée par m2.
!  Si on multiplie par DX*DY, on a le volume évaporé ou condensé par maille
!  horizontale.
!
!  La salinité est une quantité de sel (Q) par unité de volume (V): S=Q/V.
!  Soit dV la variation de volume liée à l'évaporation ou à la condensation.
!  Pour la même quantité de sel Q, la salinité est maintenant Q/(V+dV).
!  La variation de salinité DS est Q/(V+dV)-Q/V=-QdV/V/(V+dV).
!  Si dV<<V DS=-S*dV/V                                                 !26/03/06
!
!  En un temps DTI on a dV=SLHF*DTI/RHO/A * DX*DY                      !26/03/06
!
!  Le volume V de la maille est DX*DY*DZ (DZ maille verticale du niveau
!  sous la surface).
!
!  Donc DS=-S*SLHF*DTI/RHO/A/DZ est la variation de salinité de la maille
!  de surface                                                          !26/03/06
!  (K=NR-1) en un pas de temps
!
!  Comme le model_ travaille avec les quantité volumique SHZ=S*DZ on remultiplie
!  la quantité DS par DZ pour retomber sur les unités du modèle.
!
!  Quid des precipitations? On suit un raisonnement analogue au precedent.
!  Questions unités on remarque que SLHF/A/RHO est homogene à des m/s tout
!  comme PRECIPI_Z. Convention de signes: des precipitaions (>0) diminuent
!  la salinité, d'où le signe - dans la boucle 101
!..............................................................................



! standard (3D) case: salinity has to be updated from the intermediate
! "after" solution just provided by the advection subroutine:
!     if(flag3d.eq.1)kts=after  ! standard 3D case. after=2
! special 1D vertical case: advection is not computed and thus salinity
! has to be updated from the value of the previous time step
!     if(flag3d.eq.0)kts=before ! special 1D vertical case. before=0


! Special case: matrix coef have been already computed for another scalar variable:

!     tfilterfb=0.01  ! Coef filtre FB  definis dans module_principal
!     tfilterlf=0.005 ! Coef filtre LP
!     stop 'ne pas oublier tfb0=1-tfilterfb'

!$ Standard case: the coef of the tridiagonal system are now computed:
      !if(flag3d==1) then !ooo> !28-05-19
      ! thom_comment je pense mieux ainsi
      if(flag_1dv/=1) then !ooo> !31-05-19
       if(anyv3d(imax/2,jmax/2,kmax-1,id_sdiv)/=check6)stop 'check6'
      else               !ooo> !28-05-19
! Cas particulier du modele 1DV
          anyv3d(:,:,:,id_sdiv)=0.
      endif              !ooo> !28-05-19


      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Note: le filtre temporel est de la forme A(t+1)-A(t) avec A(t+1)=tfilterfb*B(t+1)+tfilterlf*C(t+1)
! B(t+1)=dz(t+1)*(S(t+1)-S(t))+dz(t)*(S(t)-S(t-1)) c.a.d. la moyenne de la derivee sur 2 temps successifs
! C(t+1)=dz(t+1)*(S(t)-S(t-1))
! De sorte que A(t+1)-A(t)=tfilterfb*[ dz(t+1)*(S(t+1)-S(t))-dz(t-1)*(S(t-1)-S(t-2)) ]
!                         +tfilterlf*[ dz(t+1)*(S(t)-S(t-1))-dz(t  )*(S(t-1)-S(t-2)) ]

! coef4 of the tridiagonal system:
      tridia_in(i,j,k,4)=(                                    &

! time filter:
!!!!!!!!     dz_t(i,j,k,before)* sal_t(i,j,k,before)          &
       dz_t(i,j,k,before)*( sal_t(i,j,k,2)                    & !05-01-16

                        +tfilterfb*( sal_t(i,j,k,before2)         &
                                    -sal_t(i,j,k,before ) ))      &

      +dz_t(i,j,k,now  )*tfilterlf*( sal_t(i,j,k,before2)         &
                                    -sal_t(i,j,k,before ))        &

      +dz_t(i,j,k,after)*( (tfilterlf-tfilterfb)*sal_t(i,j,k,now   )  &
                                     -tfilterlf *sal_t(i,j,k,before)) &

!            -wetmask_t(i,j)*anyv3d(i,j,k,id_sdiv)   &
! je supprime wetmask_t(i,j)* qui me semble inutile le !06-04-22
                            -anyv3d(i,j,k,id_sdiv)   &

!           +dz_t(i,j,k,after)*(tfc1*sal_t(i,j,k,now)        & ! t      !27-12-10
!                              +tfc2*sal_t(i,j,k,before)     & ! t-dt
!                              +tfc3*sal_t(i,j,k,before2))   & ! t-2dt
!           -dz_t(i,j,k,now)  *(tfc0*sal_t(i,j,k,now)        & ! t
!                              +tfc1*sal_t(i,j,k,before)     & ! t-dt
!                              +tfc2*sal_t(i,j,k,before2)    & ! t-2dt
!                              +tfc3*sal_t(i,j,k,before3))   & ! t-3dt

! Advection:
     ! +( anyv3d(i+1,j  ,k  ,4)                              & !05-01-16
     !   -anyv3d(i  ,j  ,k  ,4)                              &
     !   +anyv3d(i  ,j+1,k  ,5)                              &
     !   -anyv3d(i  ,j  ,k  ,5) )/dxdy_t(i,j)                &
     !   +anyv3d(i  ,j  ,k+1,6)                              & !13-11-16
     !   -anyv3d(i  ,j  ,k  ,6)                              &

          )

      enddo
      enddo
      enddo

! Retroaction du modele enfant sur le modele parent:
      if(oasis_symsym_onoff==1)call cpl_oasis_sal_delta !22-08-20

! Additional mixing if non-zero salt molecular diffusivity: !11-07-15
! Note: compte tenu de la valeur faible de kmol_s le calcul
! est fait en explicite:
      if(kmol_s>0.) then !kmkmkmkm>
       do k=1,kmax
        kp1=min(k+1,kmax)
        km1=max(k-1,1)
! Note a propos de km1: le fait que S soit homogene entre
! k=1 et k=kmin_w(i,j) garantit la condition de flux nul a
! travers le fond si kmin/=1
        do j=1,jmax
        do i=1,imax
         tridia_in(i,j,k,4)=                &
         tridia_in(i,j,k,4)+kmol_s*dti_lp*( & !oooooo>

             (  sal_t(i,j,kp1,bef)-sal_t(i,j,k,bef)) &
            /(depth_t(i,j,k+1)  -depth_t(i,j,k))     &

            -(  sal_t(i,j,k,bef)-sal_t(i,j,km1,bef)) &
            /(depth_t(i,j,k)  -depth_t(i,j,k-1))     &

                                          )   !oooooo>
        enddo
        enddo
       enddo
      endif              !kmkmkmkm>

!     if(flag_merged_levels==1) then !m°v°m>
! Melange du RHS des equations dans la couche mergEe
!      do j=1,jmax ; do i=1,imax
!       x1=0.5*( tridia_in(i,j,   kmin_w(i,j),4) &
!               +tridia_in(i,j,kmerged_t(i,j),4) )

!                tridia_in(i,j,kmerged_t(i,j),4)=x1
!                tridia_in(i,j,   kmin_w(i,j),4)=x1

!      enddo       ; enddo
!     endif                          !m°v°m>

!     call tridiagonalsolver(0,1,imax,1,jmax,kmax)    !26-06-14
      if(flagsolver_/=0)stop 'Unexpected flagsolver_ in S equation'
      call tridiagonalsolver(flagsolver_,1,imax,1,jmax,kmax)

!       do k=1,kmax
        do k=kmax,1,-1
        do j=1,jmax
        do i=1,imax
                   sal_t(i,j,k,2)=     &          !10-10-09
              tridia_out(i,j,k)        &
                 *mask_t(i,j,kmax)     &          !11-07-15
              *wetmask_t(i,j)          &          !12-01-10
          +(1.-wetmask_t(i,j))         &
                  *sal_t(i,j,k,1) !12-01-16
        enddo
        enddo
        enddo

#ifdef bilan_s_3d
      call my_outputs_3dflux_s('implicite',0) !12-02-22
      call my_outputs_3dflux_s('check',0)     !12-02-22
#endif

#ifdef bidon
!......................................................
!$ update salinity with surface water flux                             !02-03-19
      const1=1./(rho*lv)                                               !12-10-14
      do j=1,jmax                                                      !10/12/07
      do i=1,imax

! Pour garantir une parfaite coherence avec la maniere dont omega_w
! est reliee A slhf_w et precipi_w dans module_airseaflux.F90
! (je note que j'avais oubliE de reporter wetmask par exemple...)
! je commente les lignes suivantes et les reprend en remplaCant
! slhf_w(i,j,1)*const1+precipi_w(i,j,1) par -omega (attention signe oppposE)
!     x1=-dti_lp*(slhf_w(i,j,1)*const1+precipi_w(i,j,1))               &
!                 *sal_t(i,j,kmax,1)                                   &
!             /( depth_w(i,j,kmax+1)- depth_w(i,j,ksl_t(i,j)))

!     x1= dti_lp*omega_w(i,j,kmax+1,1) & ! omega_w(kmax+1)=-(slhf_w*const1+precipi_w) module_airseaflux) !27-03-16
!                 *sal_t(i,j,kmax,1)                           &
!             /( depth_w(i,j,kmax+1)- depth_w(i,j,ksl_t(i,j))) &
!     *int(wetmask_u(i,j)*wetmask_u(i+1,j)*wetmask_v(i,j)*wetmask_v(i,j+1)) !25-02-19

! Puisqu'avec les bancs decouvrants omega_w(kmax+1) peut representer autre chose que
! la vitesse verticale due au bilan E-P, on utilise A la place 0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))
! Noter qu'on utilise sal_t au temps after pour etre coherent avec le fait
! que l'etape d'advection verticale est implicite et donc basee sur sal_t au temps after...
! Details dans https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
#ifdef key_oasis_symp_surfex
       x1=dti_lp*(precipi_w(i,j,1)/rho)                         &
                  *sal_t(i,j,kmax,2)                            &
               /( depth_w(i,j,kmax+1)- depth_w(i,j,ksl_t(i,j)))   !24-01-18
#else
      x1= dti_lp*0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1)) & !28-02-19
                  *sal_t(i,j,kmax,2)                                 &
              /( depth_w(i,j,kmax+1)- depth_w(i,j,ksl_t(i,j)))
#endif
       do k=kmax,ksl_t(i,j),-1
       tridia_in(i,j,k,4)=                                             &
       tridia_in(i,j,k,4)+x1*( depth_w(i,j,k+1)- depth_w(i,j,k)) 
       enddo
      enddo
      enddo
!......................................................                !28/02/03

!     call tridiagonalsolver(0,1,imax,1,jmax,kmax)    !26-06-14
      if(flagsolver_/=0)stop 'Unexpected flagsolver_ in S equation'
      call tridiagonalsolver(flagsolver_,1,imax,1,jmax,kmax) 

! Under the virtual bottom:
!      if(vststep==1) then !sssssssss>

! Methode: simple persistence depuis le niveau kmin_w
!       do j=1,jmax
!       do i=1,imax


!        do k=kmin_w(i,j),kmax
!                  sal_t(i,j,k,2)=     &          !10-10-09
!             tridia_out(i,j,k)        &
!                *mask_t(i,j,kmax)     &          !11-07-15
!             *wetmask_t(i,j)          &          !12-01-10
!         +(1.-wetmask_t(i,j))         &
!                 *sal_t(i,j,k,1) !12-01-16
!        enddo

!        do k=1,kmin_w(i,j)-1
!                sal_t(i,j,k,2)=           &
!                sal_t(i,j,kmin_w(i,j),2)
!        enddo

!       enddo
!       enddo

!      else                !sssssssss>

        do k=1,kmax
        do j=1,jmax
        do i=1,imax
                   sal_t(i,j,k,2)=     &          !10-10-09
              tridia_out(i,j,k)        &
                 *mask_t(i,j,kmax)     &          !11-07-15
              *wetmask_t(i,j)          &          !12-01-10
          +(1.-wetmask_t(i,j))         &
                  *sal_t(i,j,k,1) !12-01-16
        enddo
        enddo
        enddo

!      endif               !sssssssss>

! Ligne commentee car melange augmente en amont dans la matrice implicite
!     if(flag_merged_levels==1)call vertmix_merged_levels_sal(2)
#endif

! Meme chose mais pour les sources sous marine
! La temperature de fond est traitEe en mEme temps
! commentees le 05-03-21 car rivieres de surface et sous marine integrees au schema d'advection
!     if(flag_groundwater==1)call obc_river_groundwater_ts !17-04-19
!       if(flag_surfriver==1)call obc_river_surfriver_tem  !19-02-21

      end subroutine vertmix_sal

!..........................................................................

      subroutine vertmix_matrix
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='vertmix_matrix'
       subroutinedescription= &
          'Computes the coefficients of the vertical mixing' &
       //' matrix'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!......................................................
! Coefficients for the following linear equations system:
! tridia_in(k,1)*X(k-1)+tridia_in(k,2)*X(k)+tridia_in(k,3)*X(k+1)=tridia_in(k,4)
! dz(after)*X(after)=X(before)*dz(before)-dti_lp*D[ Kh.DX(after)/Dz ]/Dk
!......................................................

! Melange intense dans la couche fusionnee:
       do j=1,jmax ; do i=1,imax !-i,j->
        do k=2,kmin_w(i,j) !18-03-22
! La formulation suivante revient A considerer que le coef de melange est K=10*(dz**2)/dt
! soit une valeur forte compte tenu de la resolution et du pas de temps, mais sans Etre exagErEe.
! Notamment on considere dz(k-1) plutot que dz(k) pour eviter une valeur trop grande
! entre kmin_w et kmin_w-1. Il faut aussi garder en memoire que les traceurs des sources 
! soumarines transitent par les couches ecrasEes depuis k=1
         tridia_in(i,j,k,1)=-10.*dz_t(i,j,k-1,2)  !19-03-22
        enddo
       enddo       ; enddo       !-i,j->

! Melange normal pour les autres niveaux:
       do j=1,jmax ; do i=1,imax !-i,j->
        do k=kmin_w(i,j)+1,kmax !02-11-21!18-03-22
                tridia_in(i,j,k,1)=-dti_lp*mask_t(i,j,kmax)*( & !pmx>
                     kh_w(i,j,k  )/(depth_t(i,j,k  )-depth_t(i,j,k-1)) &
         +(1-wetmask_wi_t(i,j))*5. & ! melanger couche wetdry !02-03-22
                                                         )   !pmx>
        enddo
       enddo       ; enddo       !-i,j->

! Deduire coef3 de coef1
       do k=1,kmax-1 ; do j=1,jmax ; do i=1,imax
         tridia_in(i,j,k  ,3) &
        =tridia_in(i,j,k+1,1)
       enddo       ; enddo       ; enddo

!CL fond surface
       do j=1,jmax ; do i=1,imax
         tridia_in(i,j,kmax,3)=0. !19-10-19
         tridia_in(i,j,1   ,1)=0. !02-03-22
       enddo       ; enddo       

! Deduire coef2 de coef3 coef1
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
          tridia_in(i,j,k,2)=tfb0*dz_t(i,j,k,after) &
         -tridia_in(i,j,k,1)                        &
         -tridia_in(i,j,k,3)
       enddo       ; enddo       ; enddo


!     stop 'vertmix_matrix'
      end subroutine vertmix_matrix

!......................................................................
#ifdef bidon
      subroutine vertmix_nrj
      use module_principal
      implicit none
      if(nbdom/=1)return


! Ici on calcule le bilan de grav*z*d(Dz.rhp )/dt
! Conseil: pour verifier l'equilibrage de l'advection de la
! densite avec le terme de pression, annuler toutes les autres contribution
! a la variation d T et S. Annuler Kh le melange vertical. Annuler le filtre
! d'asselin. Utiliser une equation d'etat lineraire.

      id_rhpa=2 ! identifiant densite potentielle au temps after  t+1
      id_rhpb=0 ! identifiant densite potentielle au temps before t-1

! Densite potentielle au temps t-1 (before)
      id_tem=3 ; id_sal=4 ; id_rhp=id_rhpb
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_tem)=tem_t(i,j,k,before)
       anyv3d(i,j,k,id_sal)=sal_t(i,j,k,before)
      enddo       ; enddo         ; enddo
      call equation_of_state_pot_anyv(2,imax-1,2,jmax-1)

! Densite potentielle au temps t+1 (after)
      id_tem=3 ; id_sal=4 ; id_rhp=id_rhpa
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_tem)=tem_t(i,j,k,after)
       anyv3d(i,j,k,id_sal)=sal_t(i,j,k,after)
      enddo       ; enddo         ; enddo
      call equation_of_state_pot_anyv(2,imax-1,2,jmax-1)
     
      sum1=0.
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax-1

       sum1=sum1+dxdy_t(i,j)*(grav/dti_lp)*depth_t(i,j,k)*(  &
        dz_t(i,j,k,2)*anyv3d(i,j,k,id_rhpa)                  &
       -dz_t(i,j,k,0)*anyv3d(i,j,k,id_rhpb)                )

      enddo       ; enddo         ; enddo

!     write(6,*)'grav*z*d(Dz.rhp )/dt',sum1
      write(6,*)'pot',sum1


      sum0=0.
      do k=kmax,1,-1
      do j=2,jmax-1
      do i=2,imax
       sum0=sum0+rho*veldydz_u(i,j,k,1)*presgrad_u(i,j,k,1)
      enddo
      enddo
      enddo
      do k=kmax,1,-1
      do j=2,jmax
      do i=2,imax-1
       sum0=sum0+rho*veldxdz_v(i,j,k,1)*presgrad_v(i,j,k,1)
      enddo
      enddo
      enddo

! ce terme provient du remplacement de du/dx par dw/dk+d(Dz)/dt, plus precisement il
! s'agit de p.d(Dz/dt). En remplacant Dz par d(z)/dk, on obtient p.d(dz/dk)/dt puis dz/dt*dp/dk
! soit sigma*dssh/dt*grav*rhp*Dz
! ATTENTION COMME A CE STADE rhp_t n'est surement plus rhp_t au moment du calcul du gradient
! de pression. On recalcule (dans un tableau anyv3d par securite) l'EOS potentielle
      id_rhp=2 ; id_tem=1 ; id_sal=0
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax-1
       anyv3d(i,j,k,id_tem)=tem_t(i,j,k,now)
       anyv3d(i,j,k,id_sal)=sal_t(i,j,k,now)
      enddo       ; enddo         ; enddo
      call equation_of_state_pot_anyv(2,imax-1,2,jmax-1)

      sum6=0.
      do k=1,kmax+1
      do j=2,jmax-1
      do i=2,imax-1  ! 1,imax
       sum6=sum6+dxdy_t(i,j)*mask_t(i,j,kmax)*(      &
       (depth_w(i,j,k)+h_w(i,j))/hz_w(i,j,1)*(hz_w(i,j,2)-hz_w(i,j,0))/dti_lp &   ! sigma/dti
      *(dz_t(i,j,k  ,1)*anyv3d(i,j,k  ,id_rhp)             &
       +dz_t(i,j,k-1,1)*anyv3d(i,j,k-1,id_rhp)))*grav*0.5    ! rhp.dz.grav
      enddo
      enddo
      enddo
      write(6,*)'prs',sum0-sum6

      end subroutine vertmix_nrj
#endif
!........................................................................

! ATTENTION NE PAS UTILISER CETTE ROUTINE DANS UNE ZONE OU anyv3d(:,:,:,0:1)
! EST MONOPOLISE
!...............................................................................

      subroutine vertmix_merged_levels_u(time_dz_,id_var_,id_varbef_)
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_,id_varbef_


      if(flag_merged_levels==1) then !m°v°m>
! v269 04-12-19: on ne moyenne plus le champ total mais l'ecart entre
!                le champs et vel_u(0) ce qui revient approximativement
!                A moyenner l'advection seulement

       do j=2,jmax-1 ; do i=2,imax
         sum1=0. ; sum2=0.
!        do k=kundermin_u(i,j),kmerged_u(i,j)
         do k=1,kmerged_u(i,j)
          sum1=sum1+dz_u(i,j,k,time_dz_)
!         sum2=sum2+                     anyv3d(i,j,k,id_var_) !01-05-19
! Moyenner la tendance advective !04-12-19
          sum2=sum2+anyv3d(i,j,k,id_var_)-anyv3d(i,j,k,id_varbef_)*dz_u(i,j,k,time_dz_)
         enddo
         x1=sum2/max(sum1,small1)*mask_u(i,j,kmax)
         do k=1,kmerged_u(i,j)
!         anyv3d(i,j,k,id_var_)=x1
          anyv3d(i,j,k,id_var_)=x1+anyv3d(i,j,k,id_varbef_) !04-12-19
         enddo
         do k=kmerged_u(i,j)+1,kmax
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_u(i,j,k,time_dz_)
         enddo
       enddo         ; enddo

      else                           !m°v°m>

        do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_u(i,j,k,time_dz_)
        enddo       ; enddo         ; enddo

      endif                          !m°v°m>

      end subroutine vertmix_merged_levels_u

!...............................................................................

      subroutine vertmix_merged_levels_v(time_dz_,id_var_,id_varbef_)
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_,id_varbef_

      if(flag_merged_levels==1) then !m°v°m>
! v269 04-12-19: on ne moyenne plus le champ total mais l'ecart entre
!                le champs et vel_u(0) ce qui revient approximativement
!                A moyenner l'advection seulement

        do j=2,jmax ; do i=2,imax-1
         sum1=0. ; sum2=0.
!        do k=kundermin_v(i,j),kmerged_v(i,j)
         do k=1,kmerged_v(i,j)
          sum1=sum1+dz_v(i,j,k,time_dz_)
!         sum2=sum2+                     anyv3d(i,j,k,id_var_) !01-05-19
! Moyenner la tendance advective !04-12-19
          sum2=sum2+anyv3d(i,j,k,id_var_)-anyv3d(i,j,k,id_varbef_)*dz_v(i,j,k,time_dz_)
         enddo
         x1=sum2/max(sum1,small1)*mask_v(i,j,kmax)
         do k=1,kmerged_v(i,j)
!         anyv3d(i,j,k,id_var_)=x1
          anyv3d(i,j,k,id_var_)=x1+anyv3d(i,j,k,id_varbef_) !04-12-19
         enddo
         do k=kmerged_v(i,j)+1,kmax
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_v(i,j,k,time_dz_)
         enddo
        enddo         ; enddo

      else                           !m°v°m>

        do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_v(i,j,k,time_dz_)
        enddo       ; enddo         ; enddo

      endif                          !m°v°m>

      end subroutine vertmix_merged_levels_v

!...............................................................................

      subroutine vertmix_merged_levels_any(time_dz_,id_var_)
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_

        if(flag_merged_levels==1) then !m°v°m>

        do j=1,jmax ; do i=1,imax

! Melange sur la couche dynamiquement concernEe (kundermin_t : kmerged_t)
         sum1=0. ; sum2=0.
!        do k=kundermin_t(i,j),kmerged_t(i,j)
         do k=1,kmerged_t(i,j)
          sum1=sum1+dz_t(i,j,k,time_dz_)
!         sum2=sum2+dz_t(i,j,k,time_dz_)*anyv3d(i,j,k,id_var_)
          sum2=sum2+                     anyv3d(i,j,k,id_var_) !01-05-19
         enddo

         x1=sum2/max(sum1,small1)*mask_t(i,j,kmax)

! Attribution de la moyenne, y compris pour les niveaux sous kundermin_t
         do k=1,kmerged_t(i,j)
          anyv3d(i,j,k,id_var_)=x1
         enddo
         do k=kmerged_t(i,j)+1,kmax
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_t(i,j,k,time_dz_)
         enddo

        enddo ; enddo ! i,j

        else                           !m°v°m>

         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_var_)=anyv3d(i,j,k,id_var_)  &
                                 /dz_t(i,j,k,time_dz_)
         enddo ; enddo ; enddo
        
        endif                          !m°v°m>

      end subroutine vertmix_merged_levels_any

!...............................................................................

      subroutine vertmix_merged_levels_ts(time_dz_) !01-03-19
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_


        do j=1,jmax ; do i=1,imax

! Melange sur la couche dynamiquement concernEe (kundermin_t : kmerged_t)
         sum1=0. ; sum2=0. ; sum3=0.
!        do k=kundermin_t(i,j),kmerged_t(i,j)
         do k=1,kmerged_t(i,j)
          sum1=sum1+dz_t(i,j,k,time_dz_)
          sum2=sum2+dz_t(i,j,k,time_dz_)*tem_t(i,j,k,2)
          sum3=sum3+dz_t(i,j,k,time_dz_)*sal_t(i,j,k,2)
         enddo

         x1=sum2/max(sum1,small1)*mask_t(i,j,kmax)
         x2=sum3/max(sum1,small1)*mask_t(i,j,kmax)

! Attribution de la moyenne, y compris pour les niveaux sous kundermin_t
         do k=1,kmerged_t(i,j)
          tem_t(i,j,k,2)=x1
          sal_t(i,j,k,2)=x2
         enddo

        enddo
        enddo

      end subroutine vertmix_merged_levels_ts

!...............................................................................

      subroutine vertmix_merged_levels_uv(time_dz_) !01-03-19!04-11-21
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_

! Melange kmerged-1 et kmerged proportionnel au chevauchement de la
! couche de fond virtuelle et la couche kmerged !04-11-21

      do j=2,jmax-1 ; do i=2,imax
         sum1=0. ; sum2=0.
         do k=1,kmerged_u(i,j)-1
          sum1=sum1+dz_u(i,j,k,time_dz_)
          sum2=sum2+dz_u(i,j,k,time_dz_)*vel_u(i,j,k,2)
         enddo
! Le melange avec la couche kmerged est partiel et proportionnel au
! chevauchement de la couche virtuelle de fond et la couche kmerged.    
! L'epaisseur du chevauchement est x2=dz(kmerged)-sum1, 
! avec un seuil pour eviter epaisseur negative.
! Si x2=0 (c.a.d. si dz(kmerged-1)>dz(kmerged) alors pas de melange.
         k=kmerged_u(i,j)
         x2=max(dz_u(i,j,k,time_dz_)-sum1,0.)
         x1=(sum2+x2*vel_u(i,j,k,2)) &
           /(sum1+x2)
         do k=1,kmerged_u(i,j)-1
          vel_u(i,j,k,2)=x1
         enddo
         k=kmerged_u(i,j)
         vel_u(i,j,k,2)=(x2*x1+(dz_u(i,j,k,time_dz_)-x2)*vel_u(i,j,k,2)) &
                               /dz_u(i,j,k,time_dz_)
      enddo         ; enddo

      do j=2,jmax ; do i=2,imax-1
         sum1=0. ; sum2=0.
         do k=1,kmerged_v(i,j)-1
          sum1=sum1+dz_v(i,j,k,time_dz_)
          sum2=sum2+dz_v(i,j,k,time_dz_)*vel_v(i,j,k,2)
         enddo
! Le melange avec la couche kmerged est partiel et proportionnel au
! chevauchement de la couche virtuelle de fond et la couche kmerged.    
! L'epaisseur du chevauchement est x2=dz(kmerged)-sum1, 
! avec un seuil pour eviter epaisseur negative.
! Si x2=0 (c.a.d. si dz(kmerged-1)>dz(kmerged) alors pas de melange.
         k=kmerged_v(i,j)
         x2=max(dz_v(i,j,k,time_dz_)-sum1,0.)
         x1=(sum2+x2*vel_v(i,j,k,2)) &
           /(sum1+x2)
         do k=1,kmerged_v(i,j)-1
          vel_v(i,j,k,2)=x1
         enddo
         k=kmerged_v(i,j)
         vel_v(i,j,k,2)=(x2*x1+(dz_v(i,j,k,time_dz_)-x2)*vel_v(i,j,k,2)) &
                               /dz_v(i,j,k,time_dz_)
      enddo         ; enddo


#ifdef bidon
      do j=2,jmax-1 ; do i=2,imax
         sum1=0. ; sum2=0.
         do k=kundermin_u(i,j),kmerged_u(i,j)
          sum1=sum1+dz_u(i,j,k,time_dz_)
          sum2=sum2+dz_u(i,j,k,time_dz_)*vel_u(i,j,k,2)
         enddo
         x1=sum2/max(sum1,small1)*mask_u(i,j,kmax)
         do k=1,kmerged_u(i,j)
          vel_u(i,j,k,2)=x1
         enddo
      enddo         ; enddo

      do j=2,jmax ; do i=2,imax-1
         sum1=0. ; sum2=0.
         do k=kundermin_v(i,j),kmerged_v(i,j)
          sum1=sum1+dz_v(i,j,k,time_dz_)
          sum2=sum2+dz_v(i,j,k,time_dz_)*vel_v(i,j,k,2)
         enddo
         x1=sum2/max(sum1,small1)*mask_v(i,j,kmax)
         do k=1,kmerged_v(i,j)
          vel_v(i,j,k,2)=x1
         enddo
      enddo         ; enddo
#endif

      end subroutine vertmix_merged_levels_uv

!...............................................................................


      subroutine vertmix_merged_levels_bio !07-03-19
      use module_principal ; use module_parallele
      implicit none


        do j=1,jmax ; do i=1,imax

! Melange sur la couche dynamiquement concernEe (kundermin_t : kmerged_t)
         sum1=0.
!        do k=kundermin_t(i,j),kmerged_t(i,j)
         do k=1,kmerged_t(i,j) !03-08-21
          sum1=sum1+dz_t(i,j,k,1)+dz_t(i,j,k,2)
         enddo
         sum1=max(sum1,small1)


         do vb=1,vbmax

          sum2=0.
!        do k=kundermin_t(i,j),kmerged_t(i,j)
         do k=1,kmerged_t(i,j) !03-08-21
           sum2=sum2+(dz_t(i,j,k,1)+dz_t(i,j,k,2))*bio_t(i,j,k,vb)
          enddo
          x1=(sum2/sum1)*mask_t(i,j,kmax)
! Attribution de la moyenne, y compris pour les niveaux sous kundermin_t
          do k=1,kmerged_t(i,j)
           bio_t(i,j,k,vb)=x1
          enddo

         enddo ! vb loop

        enddo ; enddo

      end subroutine vertmix_merged_levels_bio

!...............................................................................


      subroutine vertmix_merged_levels_2t_any(time_dz_,id_var_) !07-03-19
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_

      if(flag_merged_levels==1) then !m°v°m>

        do j=1,jmax ; do i=1,imax

! Melange sur la couche dynamiquement concernEe (kundermin_t : kmerged_t)
         sum1=0. ; sum2=0.
!        do k=kundermin_t(i,j),kmerged_t(i,j)
         do k=1,kmerged_t(i,j)
!         sum1=sum1+(dz_t(i,j,k,time_dz_)+dz_t(i,j,k,time_dz_-1))
!         sum2=sum2+(dz_t(i,j,k,time_dz_)+dz_t(i,j,k,time_dz_-1))*anyv3d(i,j,k,id_var_)
          sum1=sum1+0.5*(dz_t(i,j,k,time_dz_)+dz_t(i,j,k,time_dz_-1))
          sum2=sum2+anyv3d(i,j,k,id_var_) !01-05-19
         enddo

         x1=sum2/max(sum1,small1)*mask_t(i,j,kmax)

! Attribution de la moyenne, y compris pour les niveaux sous kundermin_t
         do k=1,kmerged_t(i,j)
          anyv3d(i,j,k,id_var_)=x1
!         if(i+par%timax(1)==280.and.j+par%tjmax(1)==172)write(10+par%rank,*)k,anyv3d(i,j,k,id_var_)
         enddo
         do k=kmerged_t(i,j)+1,kmax
          anyv3d(i,j,k,id_var_)=    &
          anyv3d(i,j,k,id_var_)/(0.5*(dz_t(i,j,k,time_dz_)+dz_t(i,j,k,time_dz_-1)))
         enddo

        enddo ; enddo

      else                           !m°v°m>
       do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,id_var_)=    &
          anyv3d(i,j,k,id_var_)/(0.5*(dz_t(i,j,k,time_dz_)+dz_t(i,j,k,time_dz_-1)))
       enddo       ; enddo       ; enddo
      endif                          !m°v°m>

      end subroutine vertmix_merged_levels_2t_any

!...............................................................................

      subroutine vertmix_merged_levels_divu(time_dz_,id_var_)
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_


! Note: la correction du terme de divergence est de la forme vel=vel-divterm/dz
! La moyenne  de divterm/dz est donc sum([divterm/dz]*dz)/sum(dz) soit sum(divterm)/sum(dz)
! le nouveau divterm, apres moyennage, est donc [sum(divterm)/sum(dz)]*dz

       do j=2,jmax-1 ; do i=2,imax
         sum1=0. ; sum2=0.
         do k=1,kmerged_u(i,j)
          sum1=sum1+dz_u(i,j,k,time_dz_)
          sum2=sum2+anyv3d(i,j,k,id_var_)
         enddo
         x1=sum2/max(sum1,small1)*mask_u(i,j,kmax)
         do k=1,kmerged_u(i,j)
          anyv3d(i,j,k,id_var_)=x1*dz_u(i,j,k,time_dz_)
         enddo
       enddo         ; enddo

      end subroutine vertmix_merged_levels_divu

!...............................................................................

      subroutine vertmix_merged_levels_divv(time_dz_,id_var_)
      use module_principal ; use module_parallele
      implicit none
      integer time_dz_,id_var_


! Note: la correction du terme de divergence est de la forme vel=vel-divterm/dz
! La moyenne  de divterm/dz est donc sum([divterm/dz]*dz)/sum(dz) soit sum(divterm)/sum(dz)
! le nouveau divterm, apres moyennage, est donc [sum(divterm)/sum(dz)]*dz

       do j=2,jmax   ; do i=2,imax-1
         sum1=0. ; sum2=0.
         do k=1,kmerged_v(i,j)
          sum1=sum1+dz_v(i,j,k,time_dz_)
          sum2=sum2+anyv3d(i,j,k,id_var_)
         enddo
         x1=sum2/max(sum1,small1)*mask_v(i,j,kmax)
         do k=1,kmerged_v(i,j)
          anyv3d(i,j,k,id_var_)=x1*dz_v(i,j,k,time_dz_)
         enddo
       enddo         ; enddo


      end subroutine vertmix_merged_levels_divv

!...............................................................................

      subroutine vertmix_merged_levels_obc(t_) !26-10-21
      use module_principal ; use module_parallele
      implicit none
      integer t_


        do j=1,jmax ; do i=1,imax

         sum1=0. ; sum2=0. ; sum3=0.
         do k=1,kmerged_t(i,j)
          sum1=sum1+dsig_t(i,j,k)
          sum2=sum2+dsig_t(i,j,k)*temobc_t(i,j,k,t_)
          sum3=sum3+dsig_t(i,j,k)*salobc_t(i,j,k,t_)
         enddo

         x1=sum2/max(sum1,small1)*mask_t(i,j,kmax)
         x2=sum3/max(sum1,small1)*mask_t(i,j,kmax)

         do k=1,kmerged_t(i,j)
          temobc_t(i,j,k,t_)=x1
          salobc_t(i,j,k,t_)=x2
         enddo

        enddo
        enddo

      end subroutine vertmix_merged_levels_obc

!...............................................................................
