!______________________________________________________________________
! SYMPHONIE ocean model
! release 310  - last update: 02-11-21
!______________________________________________________________________
!...............................................................................
! Version date      Description des modifications
!         18/04/01: pour éviter que le niveau de stress ne s'effondre
!                   en cas d'emergence de bruit sur la vitesse
!         19/07/01: passage à la coordonnée sigma généralisée
!         27/07/01: passage à la coordonnée sigma généralisée
!         28/08/01: suppression advection tke supposée négligeable par rapport
!                   aux termes de production / consommation et finalement
!                   consommateur de cpu & eventuellement producteur de bruit 2DX
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         03/03/03: passage à coordonnees hybride
!         09/07/03: suite et debugage
!                   ajout de la condition limite surface Craig et Banner
!         03/11/03: bilan d'energie sur terme de production de vitesse
!                   + calcul conservatif sur ce dernier
!         04/11/03: bilan d'energie sur terme de flottabilité
!                   + calcul conservatif sur ce dernier
!         17/06/04: modif sur diffusion horizontale: on multiplie
!                  les flux direction x par   mask_c(i,j,k)*  mask_c(i-1,j,k)
!                   les flux direction y par   mask_c(i,j,k)*  mask_c(i,j-1,k)
!                   pour ne pas diffuser à travers les murs.
!         14/11/04: Rappeller la density(choix=1) n'a d'interet que si on utilise
!                   l'equation d'etat lineaire sans effet de pression. Son appel
!                   dépend donc maintenant des variables EQS_STATE
!         30/03/06: appel (au choix) a TRIDIAGSYST_DBP
!                   Dissipation de TKE tend vers zero quand profil rho devient
!                   instable
!         10/05/06: fonctions compatibles avec double precision
!         02/07/06: equation re-ecrite en schema forward et advection upstream
!         21/04/07: passage à la coordonnée curviligne
!         14/03/08: Correction flux advection vertical du premier niveau
!         11/11/08: Division par facteur d'echelle verticale pris au temps t+1
!         09/03/09: Regroupement des C.L. dans une routine obc_tken pour clarifier
!                   la parallelisation
! 2009.3  02-10-09  dti_fw remplace dti
!         03-10-09  - changement de nom de routine: tkenergy.F devient
!                   turbulent_kinetic_energy.F
!                   - utilisation des nouveaux facteurs d'echelle verticale
! 2010.7  16-02-10  vitesse forward = veldxdz(2)
! 2010.8  03-05-10  Avec nouveau schema forward, plus besoin de omega(2) et veldxdz(2)
! 2010.9  06-06-10  tridiagonalsolver renommé tridiagonalsolver
! 2010.11 23-07-10  suppression de l'appel à density
! 2010.12 28-09-10  remise à nievau de l'ancienne condition aux limites de Craig & Banner
! 2010.19 14-04-11  Ajout de commentaires et d'une possibilite de supprimer "proprement"
!                   l'advection (voir explication detaillee en tête de routine)
! 2010.22 01-05-11  Forward time stepping scheme conservatif
!         06-05-11  Ajout C.L. flux d'energie à la surface calculé par ww3 et au passage
!                   correction d'un bug dans condition C.B.
! S.26    18-09-14  advection verticale implicite
!         02-11-14  retour a dti_fw....
!         23-05-15  call turbulence_k_eps_sbc(0) + distribution z de l'input de tke en surface
!         11-07-15  ajout kmol_m à km_w
!         19-08-15  Possibilite de retrouver la longueur de melange de NEMO avec routine
!                   turbulence_gaspar_nemo_l
!         05-01-16  Advection pas de temps sEparEs
!         03-02-16 - Possibilite d'une friction lineaire
!                  - call turbulence_k_eps_tkebbc(const2)  
!         26-04-16 ajout turbulence_constant_kz
!                        turbulence_driver
!         01-03-17  Production Vagues ajoutE dans anyv3d(:,:,:,id_prod)
!         02-03-17  0.64hsw devient 1.6hsw
!         21-04-17  coordonnee s-z
!         31-10-17  amenagement cas kz=constant
!         01-11-17  schema de fermeture 0 equation
!         22-03-18  correction schema 0 equation
!         26-08-18  Prise en compte de ce que certaines options font demarrer 
!                   la turbulence en iteration3d==1
!         01-09-18  ajout constant_km et constant_kh
!         03-09-18 la turbulence est deplacEe apres le calcul des moments et par
!                  consequent le terme de production turbulente est maintenant calculE
!                  avec vel(2) au lieu de vel(1)
! v245    03-02-19 correction shema 0 equation
! v246    03-02-19 ajout "tke NEMO" et "tke NEMO" neutral
!         16-02-19 call turbulence_otherproduction
!                  subroutine turbulence_mangrove 
! v255    13-05-19 debug indice vertical
! v269    10-12-19 call turbulence_rhp_linear_profil !10-12-19
!                  A l'interieur de la couche fusionnEe recontruire un profil 
!                  lineaire pour la densite
! v270    13-12-19 Bornes et CL Surface et fond et ajout km=kh=0 pour k<kmin
! v293    02-12-20 iturbulence==6
! v295    18-01-21 condition craig-banner special vagues deterministes
! v296    03-02-21 subroutine turbulence_tkeproduction !03-02-21
!         27-02-21 *mask
! v302    20-05-21 adaptation a la norme gfortran
! v310    02-11-21 reconstruction d'un profil lineaire pour la densite 
!                  supprimE le 02-11-21
!...............................................................................
!    _________                    .__                  .__                     !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      subroutine turbulence_driver !26-04-16
      use module_principal ; use module_parallele
      implicit none

! A l'interieur de la couche fusionnEe recontruire un profil 
! lineaire pour la densite: (supprimE le 02-11-21)
!     if(flag_merged_levels==1)call turbulence_rhp_linear_profil !10-12-19!02-11-21

        if(iturbulence==1)call turbulence_k_eps  

        if(    iturbulence==0 &
           .or.iturbulence==4 &
           .or.iturbulence==5 &
           .or.iturbulence==6)call turbulence_gaspar 

        if(iturbulence==2) then
           if(iteration3d<=1)call turbulence_constant_kz !31-10-17!26-08-18
        endif 

        if(iturbulence==3)call turbulence_zeroequation 


! Bornes !13-12-19
      km_w=min(km_w,50.)
      kh_w=min(kh_w,50.)

! Surface et fond: !13-12-19
      do j=1,jmax ; do i=1,imax
       km_w(i,j,kmax+1)=0.
       km_w(i,j,1:kmin_w(i,j))=0. !13-12-19
       kh_w(i,j,kmax+1)=0.
       kh_w(i,j,1:kmin_w(i,j))=0. !13-12-19
      enddo ; enddo

      end subroutine turbulence_driver

!...............................................................................

      subroutine turbulence_gaspar
      use module_principal
      implicit none

! Turbulent kinetic energy:
      call turbulence_gaspar_tke

! OBC for tke at time 'n'
      call obc_turbulence_tken

! turbulent length scales:
      if(iturbulence==0)call turbulence_gaspar_l             !09-02-19
      if(iturbulence==4)call turbulence_gaspar_nemo_l        !19-08-15!09-02-19
      if(iturbulence==5)call turbulence_gaspar_nostrat       !09-02-19
      if(iturbulence==6)call turbulence_craig_banner_nostrat !02-12-20

! vertical mixing coefficient:
      call turbulence_gaspar_kz

      end subroutine turbulence_gaspar

!...........................................................

      subroutine turbulence_gaspar_tke
      use module_principal
      use module_parallele
      implicit none
      double precision dz1_,dz2_,cst1_,cst2_
#ifdef synopsis
       subroutinetitle='turbulence_gaspar_tke'
       subroutinedescription= &
       'Computes TKE according to the Gaspar scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Advection step for tke field:
      call turbulence_adv_tken !05-01-16

! TKE production
      call turbulence_tkeproduction !03-02-21

      cst1_=dti_fw*grav/rho                                      !04/11/03
      cst2_=0.5*0.5*dti_fw ! pour demi somme     km_w...!30-09-09!02-11-14
      do 200 k=2,kmax ! Attention boucle obligatoirement croissante

      do 300 j=2,jmax-1
      do 300 i=2,imax-1


          dz2_=0.25*( dz_t(i,j,k-1,2)+dz_t(i,j,k,2) &  ! dz_w(k,t+1/2) !01-05-11
                     +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))
          dz1_=0.25*( dz_t(i,j,k-1,0)+dz_t(i,j,k,0) &  ! dz_w(k,t-1/2) !01-05-11
                     +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))
!         dz2_=0.25*( dz_t(i,j,k-1, 1)+dz_t(i,j,k, 1) &  ! dz_w(k,t+1/2) !01-05-11
!                    +dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0))
!         dz1_=0.25*( dz_t(i,j,k-1, 0)+dz_t(i,j,k, 0) &  ! dz_w(k,t-1/2) !01-05-11
!                    +dz_t(i,j,k-1,-1)+dz_t(i,j,k,-1))

      tridia_in(i,j,k,1)=-dti_fw*0.5*(km_w(i,j,k-1)+km_w(i,j,k))/dz_t(i,j,k-1,2)*mask_t(i,j,k-1)
      tridia_in(i,j,k,3)=-dti_fw*0.5*(km_w(i,j,k+1)+km_w(i,j,k))/dz_t(i,j,k  ,2)*mask_t(i,j,k-1)

      tridia_in(i,j,k,2)=-tridia_in(i,j,k,1)-tridia_in(i,j,k,3)        &
        +dz2_*(1.+dti_fw*ctke2*sqrt(tken_w(i,j,k))                  &
                                  /(tkle_w(i,j,k)+small2))

! flux advection vertical superieur:
!     x1=cst2_*(omega_w(i,j,k+1,1)+omega_w(i,j,k,1))     !  dt*omega/2
!     xy_t(i,j,2)=-(x1+abs(x1))*tken_w(i,j,k)         &   ! -dt*(omega.TKE)
!                 -(x1-abs(x1))*tken_w(i,j,k+1)

! legende de TRIDIA_IN(I,J,4):

       tridia_in(i,j,k,4)=                                         &
       (                                                           &

          tken_w(i,j,k)*dz1_                                    & ! tken_w*dz_w au temps precedent

! Advection:
#ifdef zeroadvection
         +tken_w(i,j,k)*(dz2_-dz1_)                          &
#else
!        +( xy_u(i+1,j  ,1)                                        & !05-01-16
!          -xy_u(i  ,j  ,1)                                        &
!          +xy_v(i  ,j+1,1)                                        &
!          -xy_v(i  ,j  ,1) )/dxdy_t(i,j)                          &
!          +xy_t(i  ,j  ,2)                                        &
!          -xy_t(i  ,j  ,1)                                        &
#endif

! Production par cisaillement de vitesse:
         +anyv3d(i,j,k,id_prod)                                    & !01-03-17

! note que la division par dz n'apparait pas car l'equation
! est re-multipliée par dz
         +cst1_*km_w(i,j,k)*(rhp_t(i,j,k)-rhp_t(i,j,k-1))         & !04/11/03

       )*mask_t(i,j,k-1)                                           !09/07/03

!     if(iteration3d==2) then
!      if(i==imax/2.and.j==jmax/2)write(6,*)'GROSSE BIDOUILLE TKEgaspar'
!      tridia_in(i,j,k,1)=0.
!      tridia_in(i,j,k,3)=0.
!      tridia_in(i,j,k,2)=dz2_
!      tridia_in(i,j,k,4)=tken_w(i,j,k)*dz1_  &
!        +( xy_u(i+1,j  ,1)                   &
!          -xy_u(i  ,j  ,1)                   &
!          +xy_v(i  ,j+1,1)                   &
!          -xy_v(i  ,j  ,1) )/dxdy_t(i,j) 
!     endif

!      if(iteration3d==10) then
! Lignes utiles pour verifier les proprietEs de conservation du
! schema d'advection
!       if(i==imax/2.and.j==jmax/2.and.k==kmax/2)write(6,*)'BIDOUE TKE'
!       tridia_in(i,j,k,2)=dz2_
!       tridia_in(i,j,k,1)=0.
!       tridia_in(i,j,k,3)=0.
!       tridia_in(i,j,k,4)=tken_w(i,j,k)*dz1_*mask_t(i,j,k-1)
!      endif



!     xy_t(i,j,1)=xy_t(i,j,2) !update flux advectif vertical facette inferieure
#ifndef zeroadvection
! Implicit vertical advection
      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+cst2_*(          &
               (     omega_w(i,j,k  ,1)+omega_w(i,j,k+1,1)   &
                +abs(omega_w(i,j,k  ,1)+omega_w(i,j,k+1,1))) &
              +(    -omega_w(i,j,k-1,1)-omega_w(i,j,k  ,1)   &
                +abs(omega_w(i,j,k-1,1)+omega_w(i,j,k  ,1))))

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+cst2_*           &
               (   -omega_w(i,j,k-1,1)-omega_w(i,j,k  ,1)    &
               -abs(omega_w(i,j,k-1,1)+omega_w(i,j,k  ,1)))

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+cst2_*           &
               (     omega_w(i,j,k  ,1)+omega_w(i,j,k+1,1)   &
                -abs(omega_w(i,j,k  ,1)+omega_w(i,j,k+1,1)))
#endif

  300 continue
  200 continue

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! DEBUT:
!..............................................................................
! conditions aux limites sur coefs coef1 coef3 coef4
! Changement par rapport au schéma précédent, les équations
! triviale en K=NR et K=KLOW_Z sont maintenant résolue comme
! des équations normales (sauf que les coefs sont tout simple).
! Equation à la surface
      if(tke_surf.eq.0) then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>        !09/07/03
      const3=1./(rho*sqrt(ctke1*ctke2))
      do j=2,jmax-1
      do i=2,imax-1
        tridia_in(i,j,kmax+1,1)=0.
        tridia_in(i,j,kmax+1,2)=1.
        tridia_in(i,j,kmax+1,3)=0.
        tridia_in(i,j,kmax+1,4)=const3*wstress_w(i,j)
      enddo
      enddo
      endif                  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>        !09/07/03

      if(tke_surf.eq.1) then !-------------------------------->        !09/07/03

       do j=2,jmax-1
       do i=2,imax-1
! Dans le cas de la C.L. de Craig et banner le niveau kmax+1 ne sert
! qu'à s'assurer que le flux kz*dtken/dz entre kmax et kmax+1 ne joue
! pas. On a tken(kmax+1)=tken(kmax) que nous codons avec les lignes
! suivantes dites d'equation de surface. Le flux d'energie que Craig
! et banner estiment à 100*.ustar**3 soit 100.*(wstress/rho)**1.5
! est appliqué au niveau kmax que nous codons avec les lignes suivantes
! dite de flux de C&B.
! Equation de surface tken(kmax+1)=tken(kmax):                !28-09-10
        tridia_in(i,j,kmax+1,1)=-1.
        tridia_in(i,j,kmax+1,2)= 1.
        tridia_in(i,j,kmax+1,3)=0.
        tridia_in(i,j,kmax+1,4)=0.

       enddo
       enddo

! Surface Boundary Condition (tke injection from waves field & wind):
! CommentE car deplace dans le terme anyv3d(:,:,:,id_prod) !01-03-17
!      call turbulence_k_eps_sbc(0)  !23-05-15

      endif                  !-------------------------------->        !09/07/03

! Equation au fond:
      const2=0.5/sqrt(ctke1*ctke2)                                !18/04/01
      call turbulence_k_eps_tkebbc(const2) !Equilibre Production diffusion

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! FIN.
!..............................................................................

      call tridiagonalsolver(1,2,imax-1,2,jmax-1,kmax+1)                     !30/03/06

      do k=1,kmax+1
      do j=2,jmax-1
      do i=2,imax-1

         tken_w(i,j,k)=max(tridia_out(i,j,k)*mask_t(i,j,kmaxp1),emin) !27-02-21

      enddo
      enddo
      enddo

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif

      end subroutine turbulence_gaspar_tke

!................................................................................

      subroutine turbulence_gaspar_l
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='turbulence_gaspar_l'
       subroutinedescription= &
       'Computes the length scales of the Gaspar turbulent scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!...............................................................................
! Version Date      Description des modifications
!         01/11/01: inversion de l ordre des boucles: J passe avant I
!         05/03/03: passage à la coordonnee verticale hybride
!         03/06/04: correction d un commentaire
!         18/04/06: fonctions compatibles avec double precision
!         10/05/06: fonctions compatibles avec double precision
!         09/03/09: regroupement des C.L. dans obc_tkll pour clarifier la parallelisation
!         05-06-09: Parallelisation: suppression des C.L. en calculant les longueurs
!                   de i=1:imax & j=1:jmax
! 2010.3  16-01-10  mixinglength renommée turbulent_length_scale
! 2010.12 18-09-10  z0s donné par karman*1.6*hs_wave_t (voir these Nicolas Rascle)
!         20-09-10  Possibilité de calcul en simple precision
! 2010.22 06-05-11  C.L. surface
!...............................................................................

! CONST1: raccourci calculs
      const1=0.5*grav/rho

      do 100 j=1,jmax   ! 2,NECO-1                                     !05-06-09
      do 100 i=1,imax   ! 2,MECO-1

#ifdef stokes
!     z0s=1.6*hs_wave_t(i,j,1) ! these Nicolas Rascle page 43 & Eq 2.12  !18-09-10
      z0s=1.6*hsw_wave_t(i,j,1) ! these Nicolas Rascle page 43 & Eq 2.12  !06-05-11!02-03-17
#endif

      k=kmin_w(i,j)-1                                                  !05/03/03
         rhp_t(i,j,k  )=                                              & !05/03/03
         rhp_t(i,j,k+2)-(rhp_t(i,j,k+2)-   rhp_t(i,j,k+1))*           & !05/03/03
                     ( depth_t(i,j,k+2)- depth_t(i,j,k  ))            & !05/03/03
                    /( depth_t(i,j,k+2)- depth_t(i,j,k+1))             !05/03/03

      rhp_t(i,j,kmax+1)=rhp_t(i,j,kmax-1)+(rhp_t(i,j,kmax)-rhp_t(i,j,kmax-1))*  &
                              ( depth_t(i,j,kmax+1)- depth_t(i,j,kmax-1))     &
                             /( depth_t(i,j,kmax)- depth_t(i,j,kmax-1))


!     DO 111 K=2,NR-1
      do 111 k=kmin_w(i,j)+1,kmax                                      !05/03/03
      itest=0
      lup=0.
      sumb=0.
! definition de la densite (T et S) au niveau d energie considere. De part la
! distribution "etagee", un niveau K d energie est situe entre les niveaux
! K et K-1 de densite:
      rap=( depth_w(i,j,k)- depth_t(i,j,k-1))                           &
         /( depth_t(i,j,k)- depth_t(i,j,k-1))
      rbase=(1.-rap)*(rhp_t(i,j,k-1))+rap*rhp_t(i,j,k)

!     CALCUL DE L HAUT: LUP
! calcul de l'integrale:
      do 44 k2=k,kmax+1

      if(k2.eq.k)then
! cas particulier pour K2=k, puisque la repartition en demi-niveaux
! fait debuter l'integration sur le DEMI-intervale superieur,
        suma=sumb
        sumb=sumb-const1*(rhp_t(i,j,k2)-rbase)                          &
                        *( depth_t(i,j,k2)- depth_w(i,j,k))
      else
! par la suite l'integrale s'effectue sur des intervales entiers.
        suma=sumb
        sumb=sumb-const1*(rhp_t(i,j,k2)+rhp_t(i,j,k2-1)-2.*rbase)       &
                        *( depth_t(i,j,k2)- depth_t(i,j,k2-1))
      endif

      if(sumb.ge.tken_w(i,j,k))then
! test: si l'integrale devient plus grande que l'energie,
! alors la longueur est comprise entre les niveaux K2 et K2-1.
! on sort de la boucle. L'avant derniere valeur de l'integrale
! est archivee dans SUMA.
      itest=1
      goto 55
      endif

   44 continue

! rdv itest.eq.1:
   55  if(itest.eq.1)then

! longueur (lup ou ldown) sera solution de
! (a/2)*longueur**2 + b * longueur + c = 0
! COA=a                                                                !03/06/04
! COB=b
! COC=-(TKEN_Z(I,J,K)-SUMA)/G
! le calcul, classique, passe par le calcul du discriminant:
! DISCRI = b**2-4(a/2)c=COB**2-2*COA*COC
! cas general deux solutions (sol1,sol2)=( (-b-sqrt(DISCRI))/a
!                                        , (-b+sqrt(DISCRI))/a )
! cas a=0 ou considéré comme tel longueur=-c/b

!__________________________________________________________!
! coef c du polynome (a/2)*longueur**2+b*longueur+c :      !
          coc=-(tken_w(i,j,k)-suma)/grav                   !
!                                                          !
! coef b du polynome (a/2)*longueur**2+b*longueur+c :      !
          cob=-(rhp_t(i,j,k2-1)-rbase)/rho                 !
             if(k2.eq.k) then                              !
               cob=0.                                      !
               profm1= depth_w(i,j,k)                      !
             else                                          !
               profm1= depth_t(i,j,k2-1)                   !
             endif                                         !
!                                                          !
! coef a du polynome (a/2)*longueur**2+b*longueur+c :      !
          coa=-(rhp_t(i,j,k2)-rhp_t(i,j,k2-1))/            & !
               (( depth_t(i,j,k2)- depth_t(i,j,k2-1))*rho) !
!                                                          !
! calcul du discriminant afin d en trouver les racines     !
          discri=max(cob**2-2*coa*coc,zero)                !
!__________________________________________________________!


          if(abs(coa/(cob**2+1.e-10)).lt.1.e-10)then
! cas particulier coef a est nul:
             cob=sign(max(abs(cob),small1),cob)
             sol1=min( -coc/cob , ( depth_t(i,j,k2)-profm1+z0s) )
             lup=profm1- depth_w(i,j,k)+sol1
          else
! cas general (la plus grande des racines est retenue pour solution):
             coa=sign(max(abs(coa),small1),coa)
             sol1=min((-cob+sqrt(discri))/coa,                          &
                        ( depth_t(i,j,k2)-profm1+z0s))
             sol2=min((-cob-sqrt(discri))/coa,                          &
                        ( depth_t(i,j,k2)-profm1+z0s))
             lup=profm1- depth_w(i,j,k)+                                &
                  max(max(sol1,sol2),zero)
          endif

! rdv itest.eq.1:
       else

         lup= depth_w(i,j,kmax+1)- depth_w(i,j,k)+z0s

! rdv itest.eq.1:
       endif

!     CALCUL DE L BAS  LDOWN
      ldown=0.
      itest=0
      sumb=0.

!     DO 66 K2=K,1,-1
      do 66 k2=k,kmin_w(i,j),-1                                        !05/03/03

      if(k2.eq.k) then
      suma=sumb
      sumb=sumb+const1*(rhp_t(i,j,k2-1)-rbase)                          &
                      *( depth_w(i,j,k)- depth_t(i,j,k2-1))
      else
      suma=sumb
      sumb=sumb+const1*(rhp_t(i,j,k2)+rhp_t(i,j,k2-1)-2.*rbase)         &
                      *( depth_t(i,j,k2)- depth_t(i,j,k2-1))
      end if

      if(sumb.ge.tken_w(i,j,k))then
      itest=1
      goto 77
      endif

   66 continue

! rdv itest.eq.1:
   77  if(itest.eq.1)then

! longueur (lup ou ldown) sera solution de
! (a/2)*longueur**2 + b * longueur + c = 0
! COA=a                                                                !03/06/04
! COB=b
! COC=-(TKEN_Z(I,J,K)-SUMA)/G
! le calcul, classique, passe par le calcul du discriminant:
! DISCRI = b**2-4(a/2)c=COB**2-2*COA*COC
! cas general deux solutions (sol1,sol2)=( (-b-sqrt(DISCRI))/a
!                                        , (-b+sqrt(DISCRI))/a )
! cas a=0 ou considéré comme tel longueur=-c/b

!__________________________________________________________!
! coef c du polynome (a/2)*longueur**2+b*longueur+c :      !
          coc=-(tken_w(i,j,k)-suma)/grav                   !
!                                                          !
! coef b du polynome (a/2)*longueur**2+b*longueur+c :      !
          cob=(rhp_t(i,j,k2)-rbase)/rho                    !
            if(k2.eq.k) then                               !
              cob=0.                                       !
              profp1= depth_w(i,j,k)                       !
            else                                           !
              profp1= depth_t(i,j,k2)                      !
            endif                                          !
!                                                          !
! coef a du polynome (a/2)*longueur**2+b*longueur+c :      !
          coa=-(rhp_t(i,j,k2)-rhp_t(i,j,k2-1))/            & !
               (( depth_t(i,j,k2)- depth_t(i,j,k2-1))*rho) !
!                                                          !
! calcul du discriminant afin d en trouver les racines     !
          discri=max( cob**2-2*coa*coc , zero )            !
!__________________________________________________________!

          if(abs(coa/(cob**2+1.e-10)).lt.1.e-10)then
! cas particulier coef a est nul:
             cob=sign(max(abs(cob),small1),cob)
             sol1=min( -coc/cob, profp1- depth_t(i,j,k2-1) )
             ldown= depth_w(i,j,k)-profp1+sol1
          else
! cas general (la plus grande des racines est retenue pour solution):
             coa=sign(max(abs(coa),small1),coa)
             sol1=min((-cob+sqrt(discri))/coa,                          &
             (profp1- depth_t(i,j,k2-1)))

             sol2=min((-cob-sqrt(discri))/coa,                          &
             (profp1- depth_t(i,j,k2-1)))

             ldown= depth_w(i,j,k)-profp1+max(max(sol2,sol1),zero)
          endif

! rdv itest.eq.1:
       else

         ldown= depth_w(i,j,k)- depth_w(i,j,kmin_w(i,j))               !05/03/03

! rdv itest.eq.1:
       endif

      tkll_w(i,j,k)=min(ldown,lup)
!cc  1             *0.5+0.5*TKLL_Z(I,J,K)     ! modif 24/04/01 filtre sur temps
      tkle_w(i,j,k)=((abs(ldown*lup))**0.5)
!cc  1             *0.5+0.5*TKLE_Z(I,J,K)     ! modif 24/04/01 filtre sur temps

!     tkll_w(i,j,k)=1./(1./ldown+1./lup)
!     tkll_w(i,j,k)=tkle_w(i,j,k)

!.......................................................................
! Pour retrouver schema NEMO (voir Eq9 Reffray et al, Ocean Science 2015)
!      x2=-(grav/rho)*(rhp_t(i,j,k)  -rhp_t(i,j,k-1)) &  ! N**2
!                  /(depth_t(i,j,k)-depth_t(i,j,k-1))
!      x3=sqrt(4.*tken_w(i,j,k)/ (small1+x2+abs(x2) )) ! l=sqrt(2k)/N
!.......................................................................

  111 continue

      tkll_w(i,j,kmin_w(i,j))=z0_w(i,j)                                !05/03/03
      tkll_w(i,j,kmax+1     )=z0s
      tkle_w(i,j,kmin_w(i,j))=z0_w(i,j)                                !05/03/03
      tkle_w(i,j,kmax+1     )=z0s

  100 continue

      end subroutine turbulence_gaspar_l

!...........................................................................

      subroutine turbulence_gaspar_nemo_l  !19-08-15
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='turbulence_gaspar_nemo_l'
       subroutinedescription= &
       'Computes the turbulent lenght scale as in TKE NEMO'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Reference:
! Equation (9) Reffray et al 2015: Geosci. Model Dev., 8, 69–86, 2015
!                                  www.geosci-model-dev.net/8/69/2015/
!                                  doi:10.5194/gmd-8-69-2015


      x0=1./sqrt(grav/rho)
      do k=1,kmax+1
      do j=1,jmax  
      do i=1,imax   

      x1=x0*sqrt(-4.*tken_w(i,j,k)*   (depth_t(i,j,k)-depth_t(i,j,k-1))  &
                                  /(     rhp_t(i,j,k)  -rhp_t(i,j,k-1)   &
                                    -abs(rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
                                    -small1 ) ) 

! Lk=min(lup,ldown):
! ldown=min(z-zbot ,x1)
! lup  =min(zsurf-z,x1)
               tkll_w(i,j,k)=                               &
      min(min(depth_w(i,j,k)     -depth_w(i,j,1)+z0_w(i,j), & 
              depth_w(i,j,kmax+1)-depth_w(i,j,k)+z0s)     , & 

                                 x1)

! Leps=sqrt(lup*ldown):
               tkle_w(i,j,k)=   & !09-02-19
       sqrt( & !ooo>
             min(depth_w(i,j,k)     -depth_w(i,j,1)+z0_w(i,j),x1) &
            *min(depth_w(i,j,kmax+1)-depth_w(i,j,k)+z0s      ,x1) &
           )   !ooo>

      enddo
      enddo
      enddo

      end subroutine turbulence_gaspar_nemo_l

!...........................................................................

      subroutine turbulence_gaspar_nostrat  !09-02-19
      use module_principal
      implicit none

! Reference:
! No stratification case of equation (9) Reffray et al 2015: 
! Geosci. Model Dev., 8, 69–86, 2015 www.geosci-model-dev.net/8/69/2015/ 
! doi:10.5194/gmd-8-69-2015


      do k=1,kmax+1
      do j=1,jmax  
      do i=1,imax   

! Lk=min(lup,ldown):
! ldown=z-zbot
! lup  =zsurf-z
               tkll_w(i,j,k)=                               &
          min(depth_w(i,j,k)     -depth_w(i,j,1)+z0_w(i,j), &
              depth_w(i,j,kmax+1)-depth_w(i,j,k)+z0s)        

! Leps=sqrt(lup*ldown):
               tkle_w(i,j,k)=                                         & !09-02-19
       sqrt( (depth_w(i,j,k)     -depth_w(i,j,1)+z0_w(i,j)) &
            *(depth_w(i,j,kmax+1)-depth_w(i,j,k)+z0s)      )

      enddo
      enddo
      enddo

      end subroutine turbulence_gaspar_nostrat

!...........................................................................

      subroutine turbulence_craig_banner_nostrat  !02-12-20
      use module_principal
      implicit none

! Reference:
! Craig Banner 1994, JPO, vol 24, pp 2546-2559

      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax   

! Lk=min(lup,ldown):
! ldown=z-zbot
! lup  =zsurf-z
               tkll_w(i,j,k)=                               &
          min(depth_w(i,j,k)     -depth_w(i,j,1)+z0_w(i,j), &
              depth_w(i,j,kmax+1)-depth_w(i,j,k)+z0s)        

! Craig Banner 1994, JPO, vol 24, pp 2546-2559
               tkle_w(i,j,k)=tkll_w(i,j,k)

      enddo ; enddo ; enddo

      end subroutine turbulence_craig_banner_nostrat

!...........................................................................

      subroutine turbulence_gaspar_kz
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='turbulence_gaspar_kz'
       subroutinedescription= &
       'Computes the vertical coef mixing according to the Gaspar' &
       //' turbulent scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!..............................................................................
!Version  date      Description des modifs
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         03/11/03: masquage par multiplication directe par mask
!         28/11/03: application d'un seuil
!         30/03/06: diminution de ce seuil à 1E-7. (avant seuil=VIS)
!         10/05/06: fonctions compatibles avec double precision
!         04/07/06: reorganisation des boucles, possibilie de melange
!                   "accelere" en situation drho/dz instable, suppression
!                   du seuil (on seuille tke), supression du masque sur
!                       kz_z puisque deja applique sur tke.
!         04/08/06: L'augmentation artificielle du coef de melange vertical
!                   en situation convective est bornée par un critere de
!                   stabilite: Kmax=delta_z**2/DTI/2
! 2009.3  30-09-09: difver_x et difver_y supprimés et remplacés par des 1/2 sommes
!                   de     kz_z dans uv_upd.F
! 2010.3  16-01-10: routine difver_upd renommée vertmix_coef
! 2010.25 08-02-12  kz_w renomme km_w et ajout kh_w
!..............................................................................


! note sur masquage. Telle qu'est construite la grille hybride le masque
! des quantités u,v,t,s s'applique à tke et cie avec le même numéro
! d'indice vertical.

      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax

!         km_w(i,j,k)=ctke1*tkll_w(i,j,k)*sqrt(tken_w(i,j,k))
!         kh_w(i,j,k)=km_w(i,j,k)
          kh_w(i,j,k)=ctke1*tkll_w(i,j,k)*sqrt(tken_w(i,j,k)) !11-07-15
          km_w(i,j,k)=kh_w(i,j,k)+kmol_m                      !11-07-15

      enddo
      enddo
      enddo

      if(convect_yn.eq.2) then !------------------>
! Augmentation artificielle de K quand drho/dz>0:
       const2=1.e8
       do j=1,jmax
       do i=1,imax
        if(  mask_t(i,j,kmax+1).eq.1) then !!!!!!!!!!!!!!>
         do k=kmin_w(i,j)+1,kmax
              kh_w(i,j,k)=min(    kh_w(i,j,k)                           &
               +const2*max(zero,(rhp_t(i,j,k)   -rhp_t(i,j,k-1))        &
                            /( depth_t(i,j,k)- depth_t(i,j,k-1)))       &
         ,0.5*( depth_t(i,j,k)- depth_t(i,j,k-1))**2/dti_lp)
         enddo
        endif                          !!!!!!!!!!!!!!>
       enddo
       enddo
      endif                    !------------------>

      end subroutine turbulence_gaspar_kz

!............................................................................................

      subroutine turbulence_gaspar_allocate_init
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='turbulence_gaspar_allocate_init'
       subroutinedescription= &
       'Arrays allocation required by the Gaspar turbulent scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(initial_main_status==1)then !::::::::::>
       write(6,*)'If Restart files have been loaded it'
       write(6,*)'is too late for dynamic allocation in routine'
       write(6,*)'turbulence_gaspar_allocate_init'
       stop 'stop in turbulence_gaspar_allocate_init'
      endif                          !::::::::::>

      allocate(tken_w(1:imax,1:jmax,1:kmax+1)) ; tken_w=0 
      allocate(tkle_w(1:imax,1:jmax,1:kmax+1)) ; tkle_w=0
      allocate(tkll_w(1:imax,1:jmax,1:kmax+1)) ; tkll_w=0

       do k=1,kmax+1
        do j=1,jmax
        do i=1,imax
           km_w(i,j,k)=vis
           kh_w(i,j,k)=vis
         tken_w(i,j,k)=emin
         tkll_w(i,j,k)=km_w(i,j,k)/(ctke1*sqrt(tken_w(i,j,k)))
         tkle_w(i,j,k)=tkll_w(i,j,k)
        enddo
        enddo
       enddo

      end subroutine turbulence_gaspar_allocate_init

!...........................................................................

      subroutine turbulence_constant_kz !26-04-16
      use module_principal
      implicit none

! Constant Kz:

          kh_w=constant_kh !01-09-18
          km_w=constant_km

      end subroutine turbulence_constant_kz

!...........................................................................

      subroutine turbulence_zeroequation
      use module_principal ; use module_q  ; use module_parallele
      implicit none

! Km=(karman*DistanceFromSurfaceBottomBoundary)**2 * sqrt( (du/dz)**2+(dv/dz)**2 )

       do k=2,kmax ; do j=1,jmax ; do i=1,imax

         km_w(i,j,k)=                                          &

! (karman*min( (Zsurf-Z) , (Z-Zbot) ))**2
       (0.4*min( depth_w(i,j,kmax+1)-depth_w(i,j,k)             &
                ,depth_w(i,j,k     )-depth_w(i,j,1)))**2        & !22-03-18 !03-02-19
! sqrt( (du/dz)**2+(dv/dz)**2 ):
         *sqrt(0.5*( & ! OOO>
               ( (  vel_u(i  ,j,k,2) -vel_u(i  ,j,k-1,2) )      &
                /(depth_u(i  ,j,k)- depth_u(i  ,j,k-1)   ))**2  &
              +( (  vel_u(i+1,j,k,2) -vel_u(i+1,j,k-1,2) )      &
                /(depth_u(i+1,j,k)- depth_u(i+1,j,k-1)   ))**2  &

              +( (  vel_v(i,j  ,k,2) -vel_v(i,j  ,k-1,2) )      &
                /(depth_v(i,j  ,k)- depth_v(i,j  ,k-1)   ))**2  &
              +( (  vel_v(i,j+1,k,2) -vel_v(i,j+1,k-1,2) )      &
                /(depth_v(i,j+1,k)- depth_v(i,j+1,k-1)   ))**2  &
               )   )   ! OOO>

       enddo   ;  enddo   ;  enddo

       kh_w=km_w

      end subroutine turbulence_zeroequation

!...........................................................................

      subroutine turbulence_otherproduction !16-02-19
      use module_principal ; use module_parallele ; use module_q
      implicit none

! additional production terms from wave field: !01-03-17
! Deux cas appellent turbulence_k_eps_sbc: !02-12-20
! 1: le cas hydrostatique avec tke_surf==1
! 2: le cas non-hydrostatique avec tke_surf==1 avec time splitting flag_nh3d_timesplit_tsuv
      if(tke_surf==1) then !pmx>
       if(flag_nh3d==flag_nh3d_none.or.       & !02-12-20!18-01-21
          flag_nh3d==flag_nh3d_timesplit_tsuv)call turbulence_k_eps_sbc(0)
      endif                !pmx>

! Modelisation deterministe des vagues: deferlement
!    if(flag_nh3d==flag_nh3d_nosplit_uv.or. & !02-12-20
!       flag_nh3d==flag_nh3d_nosplit_tsuv)call q_tken_breaking

! additional production terms from mangrove friction:
      if(coef_diss_mangrove>0.)call turbulence_mangrove !16-02-19

      end subroutine turbulence_otherproduction

!...........................................................................

      subroutine turbulence_mangrove !16-02-19
      use module_principal ; use module_parallele
      implicit none

! Ajout d'un terme de production d'energie cinetique turbulente equilibrant la perte 
! d'energie cinetique par la friction "de la mangrove". 
! Details dans:
! https://docs.google.com/document/d/1V9cnPdDir-x3Ic3DjFHDpJ8osOkKobqs4pt5z4iHTPU/edit

! Terme de friction dans l'equation des vitesses: -r*u et -r*v
! dont on deduit un terme de production de turbulence +r*(u**2+v**2)

! Le tableau anyv3d(i,j,k,id_prod) contient DZ*r*(u**2+v**2) donc on
! multiplie par le facteur d'echelle verticale 0.5*(depth_w(k+1)-depth_w(k-1))

! x0 contient les moyennes:
! 1/2   puisque 0.5*(depth_w(k+1)-depth_w(k-1))
! 1/2   puisque 0.5*(r_mangrove(i,j,k  )+r_mangrove(i,j,k+1))
! 1/16  puisque (  (u(i+1,k)+u(i+1,k-1)+u(i,k)+u(i,k-1))/4 )**2
! soit 1/64 au total

! Note: les valeurs absolues pour eviter effets compensatoires du mode numerique (c.a.d. +1-1=0)

        x0=1./64.
        do k=2,kmax
        do j=2,jmax-1
        do i=2,imax-1

        anyv3d(i,j,k,id_prod)=             &
        anyv3d(i,j,k,id_prod)+             &

         x0*( depth_w(i,j,k+1)             & ! DZ
             -depth_w(i,j,k-1))            &

        *( r_mangrove(i,j,k  )             & ! r
          +r_mangrove(i,j,k-1))            & !13-05-19

        *( & !ooo>
          ( abs(vel_u(i+1,j  ,k  ,2))      & ! u**2
           +abs(vel_u(i+1,j  ,k-1,2))      &
           +abs(vel_u(i  ,j  ,k  ,2))      &
           +abs(vel_u(i  ,j  ,k-1,2)) )**2 &

         +( abs(vel_v(i  ,j+1,k  ,2))      & ! v**2
           +abs(vel_v(i  ,j+1,k-1,2))      &
           +abs(vel_v(i  ,j  ,k  ,2))      &
           +abs(vel_v(i  ,j  ,k-1,2)) )**2 &
         )   !ooo>

        enddo
        enddo
        enddo

      end subroutine turbulence_mangrove

!...........................................................................

      subroutine turbulence_rhp_linear_profil !10-12-19
      use module_principal ; use module_parallele
      implicit none

! A l'interieur de la couche fusionnEe recontruire un profil 
! lineaire pour la densite

! Profil lineaire de densite dans la couche fusionnee
         do j=1,jmax ; do i=1,imax
! x1 est le gradient entre la couche moyenne composee de kmin et
! kmerged et la couche au dessus de kmerged
          x1=(  rhp_t(i,j,kmerged_t(i,j)+1)                            &
              -( dz_t(i,j,   kmin_w(i,j),2)*rhp_t(i,j,   kmin_w(i,j))  &
                +dz_t(i,j,kmerged_t(i,j),2)*rhp_t(i,j,kmerged_t(i,j))) &
              /( dz_t(i,j,   kmin_w(i,j),2)                            &
                +dz_t(i,j,kmerged_t(i,j),2)))                          &
            /(      depth_t(i,j,kmerged_t(i,j)+1)                      &
              -0.5*(depth_w(i,j,kmin_w(i,j))+depth_w(i,j,kmerged_t(i,j)+1)))

! Ensuite ce gradient moyen est utilise pour reconstruire une
! stratification  dans la couche kmin & kmerged
          rhp_t(i,j,kmerged_t(i,j))=rhp_t(i,j,kmerged_t(i,j)+1) &
         -x1*(depth_t(i,j,kmerged_t(i,j)+1)-depth_t(i,j,kmerged_t(i,j)))

          rhp_t(i,j,kmin_w(i,j))=rhp_t(i,j,kmerged_t(i,j)+1) &
         -x1*(depth_t(i,j,kmerged_t(i,j)+1)-depth_t(i,j,kmin_w(i,j)))

         enddo       ; enddo

      end subroutine turbulence_rhp_linear_profil

!...........................................................................

      subroutine turbulence_tkeproduction
      use module_principal
      use module_parallele
      implicit none
      double precision cst2_


      cst2_=0.5*0.5*dti_fw ! pour demi somme     km_w...!30-09-09!02-11-14

!...................
! terme de production
      do k=2,kmax 
        do j=2,jmax
        do i=2,imax
! production par cisaillement de vitesse:
! Note on divize par dz et non dz**2 car l'equation
! de la tke est multipliee par dz. Pas le même dz, certes,
! mais justement c'est plus consistant avec le bilan d'energie.
! On note egalement qu'on demi/somme le kz au point de vitesse
! afin d'être coherent avec le melange vertical de l'equation des vitesses
          xy_u(i,j,2)=cst2_*(   km_w(i,j,k)+    km_w(i-1,j,k))*      &  ! dt/2*Kz*Dz*(du/dz)**2
                            (  vel_u(i,j,k,2) -vel_u(i,j,k-1,2) )*   &
                            (  vel_u(i,j,k,2) -vel_u(i,j,k-1,2) )/   &
                              (depth_u(i,j,k)- depth_u(i,j,k-1)   )

          xy_v(i,j,2)=cst2_*(   km_w(i,j,k)+    km_w(i,j-1,k))*      &  ! dt/2*Kz*Dz*(dv/dz)**2
                            (  vel_v(i,j,k,2) -vel_v(i,j,k-1,2) )*   &
                            (  vel_v(i,j,k,2) -vel_v(i,j,k-1,2) )/   &
                            (depth_v(i,j,k)- depth_v(i,j,k-1)   )
        enddo 
        enddo 
        do j=2,jmax-1
        do i=2,imax-1
          anyv3d(i,j,k,id_prod)    & !01-03-17
         =xy_u(i  ,j  ,2)          &
         +xy_u(i+1,j  ,2)          &
         +xy_v(i  ,j  ,2)          &
         +xy_v(i  ,j+1,2)
        enddo 
        enddo 
      enddo ! fin de boucle sur k
! additional production terms from wave fields, mangroves,...       !16-02-19
      call turbulence_otherproduction                               !16-02-19
!...................

     end subroutine turbulence_tkeproduction
