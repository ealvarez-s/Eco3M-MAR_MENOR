










      subroutine turbulence_k_eps_allocate_init
!______________________________________________________________________
! SYMPHONIE ocean model
! release 310 - last update: 03-11-21
!______________________________________________________________________
! Version date     Description des modifications
! S.26    14-02-14 tableaux tkeb et epsb pour besoins du subcycling
!         18-09-14 advection verticale implicite
!         02-11-14 wetmask_t et retour a dti_fw
!         23-05-15 call turbulence_k_eps_sbc(case_) Surface Boundary Condition 
!         11-07-15 ajout kmol_m a km_w
!         22-07-15 z0s determine la profondeur d'application de la condition de 
!                  surface de la TKE
!         05-08-15 Suite aux pb de precisions de tridiagonalsolver lorsque Kz
!                  atteind des valeurs enormes hors gamme, application d'un seuil
!         11-08-15 choix possible entre les fonctions de stabilite de Galperin 1988
!                  Kantha et Claysson 1994 et Canuto 2001
!         19-08-15 Choix de 3 conditions aux limites (fond et surface) pour epsilon
!         05-01-16 Advection separee 
!         03-02-16 - Possibilite d'une friction lineaire
!                  - call turbulence_k_eps_tkebbc(const2)    
!         14-02-16 Pas de calcul de l'advection si modele 1DV
!         10-03-16 Seuillage de l'argument de la fonction exp
!         06-08-16 subroutines d'advection transferrees dans turbulence_adv.F90
!         16-09-16 C.L. sup et inf sur epsilon 
!         16-12-16 Ajout d'un choix possible pour une C.L. type CB94 pour epsilon
!                  reference Reffray 2015
!         19-12-16 - on s'oriente vers une distribution lineaire de la QDM (attention le
!                  cas des vagues reste encore A traiter).
!                  - d'autre part modif sur la C.L. sur epsilon
!         23-12-16  momentum_input_depth est la profondeur de penetration
!                   des inputs de QDM, defini dans notebook_visco
!         24-02-17 Cas du forcage par WW3, profondeur de distribution du FOC indexEe
!                  sur 0.64*hs_wave_t 
!         01-03-17 Production Vagues ajoutE dans anyv3d(:,:,:,id_prod)
!         02-03-17 0.64*hsw devient 1.6*hsw si produit par 0.4 vient ensuite
!         10-03-17 constante de Gaspar dans C.L sur epsilon
!         21-04-17 coordonnee s-z
!         03-09-18 la turbulence est deplacEe apres le calcul des moments et par
!                  consequent le terme de production turbulente est maintenant calculE
!                  avec vel(2) au lieu de vel(1)
!         12-01-19 conditions aux limites de fond pour espilon revues pour le cas
!                  des couches merged et au passage le choix des constantes est plus
!                  coherent avec le schema k-epsilon
! v245    01-02-19 possibilite d'ajouter une limite sur epsilon dans le cas non-stratifie
! v246    16-02-19 call turbulence_otherproduction
! v269    06-12-19 Dans la condition de fond prodi la vitesse de fond est la vitesse 
!                  moyennee entre z=-h et z=-h+dz(kmerged)
!         11-12-19 suite point precedent
! v270    13-12-19 CL Surface et fond deplacees dans driver gaspar !13-12-19
! v275    03-03-20 Parce que les couches fusionnEes n'annulent pas mask_t, la multiplication par
!                  mask_t des coef de la matrice ne suffit pas A produire une solution correcte
!                  dans les couches ecrasees. La condition de fond est donc programmee de k=1 A
!                  kmin_w !03-03-20
! v296    03-02-21 call turbulence_tkeproduction
!         27-02-21 *mask_t
! v310    03-11-21 adaptation A la possibilitE que dz_t(i,j,kmin_w(i,j),1) > dz_t(i,j,kmerged_t(i,j),1))
!...............................................................................
!    _________                    .__                  .__                     !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      use module_principal
      implicit none

      if(initial_main_status==1)then !::::::::::>
       write(6,*)'If Restart files have been loaded it'
       write(6,*)'is too late for dynamic allocation in routine'
       write(6,*)'turbulence_k_eps_allocate_init'
       stop 'stop in turbulence_k_eps_allocate_init'
      endif                          !::::::::::>

      allocate(tkeb_w(1:imax,1:jmax,1:kmax+1)) ; tkeb_w=0.
      allocate(tken_w(1:imax,1:jmax,1:kmax+1)) ; tken_w=0.
      allocate(tkea_w(1:imax,1:jmax,1:kmax+1)) ; tkea_w=0.
      allocate(epsb_w(1:imax,1:jmax,1:kmax+1)) ; epsb_w=0
      allocate(epsn_w(1:imax,1:jmax,1:kmax+1)) ; epsn_w=0
      allocate(epsa_w(1:imax,1:jmax,1:kmax+1)) ; epsa_w=0

       do k=1,kmax+1
        do j=1,jmax
        do i=1,imax

           km_w(i,j,k)=vis
           kh_w(i,j,k)=vis
         tken_w(i,j,k)=emin
         epsn_w(i,j,k)=0.1*tken_w(i,j,k)**2/kh_w(i,j,k)

         tkea_w(i,j,k)=tken_w(i,j,k)
         epsa_w(i,j,k)=epsn_w(i,j,k)
         tkeb_w(i,j,k)=tken_w(i,j,k)
         epsb_w(i,j,k)=epsn_w(i,j,k)

        enddo
        enddo
       enddo

      end subroutine turbulence_k_eps_allocate_init

!............................................................................................

      subroutine turbulence_k_eps
      use module_principal
      use module_parallele
      implicit none
!     double precision dz1_,dz2_,cst1_,cst2_,cst3_
      double precision dz1_,     cst1_,cst2_,cst3_
      real*4,parameter ::         &
!     real*4           ::         &
                g0_=0.6674699  &
               ,g1_=15.62      &
               ,g2_=0.3932723  &
               ,g3_=3.085786   &
               ,g4_=34.67640   &
               ,g5_=6.1272     &
               ,g6_=0.4939277  &
              ,kc1_=0.4939277  &
              ,kc2_=-30.19200  &
              ,kc3_=-6.127200  &
              ,kc4_=0.3932723  &
              ,kc5_=17.07336   &
            ,cst_k_=1.3        &
         ,cst_prod_=1.44       &
          ,cst_eps_=1.92       &
        ,cst_buo_p_=1.         &
        ,cst_buo_m_=-0.52        ! -0.404
!       ,cst_c0cub=0.1704       ! 0.5544**3
!     integer loop_,spinup_


! TKE Equation:

! Advection step for tke field:
      call turbulence_adv_tkea !05-01-16

! TKE production
      call turbulence_tkeproduction !03-02-21
      check1=anyv3d(imax/2,jmax/2,kmax/2,id_prod)
       
      gravoverrho=grav/rho
      cst1_=dti_fw*grav/rho
      cst2_=0.5*0.5*dti_fw        !02-11-14
      cst3_=0.5*dti_fw

      do 200 k=2,kmax ! An increasing k loop is required

      do 300 j=2,jmax-1
      do 300 i=2,imax-1

! Buoyancy term:
       anyv3d(i,j,k,id_buoy)=cst1_*kh_w(i,j,k)*(rhp_t(i,j,k)-rhp_t(i,j,k-1))


!         dz2_=0.25*( dz_t(i,j,k-1,2)+dz_t(i,j,k,2) &  ! dz_w(k,t+1/2) !01-05-11
!                    +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))
!         dz1_=0.25*( dz_t(i,j,k-1,0)+dz_t(i,j,k,0) &  ! dz_w(k,t-1/2) !01-05-11
!                    +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))

! On exprime les facteurs d'echelle en passant par le tableau depth_w qui est "ajustE" dans la couche fusionnee
! en k=kmerged pour eviter la singularitE constateE. On renonce ce faisant A la definition rigoureuse de dz1 qui
! garantissait strictement la conservation de l'advection, mais ce n'est pas grave.
      dz1_=0.5*(depth_w(i,j,k+1)-depth_w(i,j,k-1)) !06-12-19

!     tridia_in(i,j,k,1)=-dti_fw*0.5*(km_w(i,j,k-1)+km_w(i,j,k))/dz_t(i,j,k-1,2)*mask_t(i,j,k-1)
      tridia_in(i,j,k,1)=-dti_fw*0.5*(km_w(i,j,k-1)+km_w(i,j,k))/(depth_w(i,j,k)-depth_w(i,j,k-1))*mask_t(i,j,k-1) !06-12-19

!     tridia_in(i,j,k,3)=-dti_fw*0.5*(km_w(i,j,k+1)+km_w(i,j,k))/dz_t(i,j,k  ,2)*mask_t(i,j,k-1)
      tridia_in(i,j,k,3)=-dti_fw*0.5*(km_w(i,j,k+1)+km_w(i,j,k))/(depth_w(i,j,k+1)-depth_w(i,j,k))*mask_t(i,j,k-1) !06-12-19

      tridia_in(i,j,k,2)=dz1_-tridia_in(i,j,k,1)-tridia_in(i,j,k,3) &
!       +dti_fw*dz2_*(epsn_w(i,j,k)+epsb_w(i,j,k))                   & ! epsilon implicite
!                      /(tken_w(i,j,k)+tkeb_w(i,j,k))
!       +dti_fw*dz2_*(epsn_w(i,j,k)              )                   & ! epsilon implicite
!                     /(tken_w(i,j,k)              )
        +dti_fw*dz1_*0.5*(epsn_w(i,j,k)/tken_w(i,j,k)              & !mm2
                         +epsb_w(i,j,k)/tkeb_w(i,j,k))


! Explicit right hand side terms of the tke equation:
       tridia_in(i,j,k,4)=                                         &
       (                                                           &

!         tken_w(i,j,k)*dz1_                                    & ! tken_w*dz_w au temps precedent
          tkea_w(i,j,k)*dz1_  & ! Prendre tkea si etape advection preliminaire !05-01-16


! Advection:

! Production par cisaillement de vitesse:
         +anyv3d(i,j,k,id_prod)                             &

! note que la division par dz n'apparait pas car l'equation
! est re-multipli�e par dz
         +anyv3d(i,j,k,id_buoy)                                          &

! Epsilon explicite:
!       -dti_fw*dz2_*epsn_w(i,j,k)                                &

       )*mask_t(i,j,k-1)                                           !09/07/03

!      if(iteration3d==10) then
! Lignes utiles pour verifier les proprietEs de conservation du
! schema d'advection
!       if(i==imax/2.and.j==jmax/2.and.k==kmax/2)write(6,*)'BIDOUE TKE'
!       tridia_in(i,j,k,2)=dz2_
!       tridia_in(i,j,k,1)=0.
!       tridia_in(i,j,k,3)=0.
!       tridia_in(i,j,k,4)=tkea_w(i,j,k)*dz1_*mask_t(i,j,k-1)
!      endif

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

  300 continue
  200 continue
      check2=anyv3d(imax/2,jmax/2,kmax/2,id_buoy)

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! DEBUT:
!..............................................................................
! conditions aux limites sur coefs coef1 coef3 coef4
! Changement par rapport au sch�ma pr�c�dent, les �quations
! triviale en K=NR et K=KLOW_Z sont maintenant r�solue comme
! des �quations normales (sauf que les coefs sont tout simple).
! Equation � la surface

      if(tke_surf.eq.0) then !------------------------------->        !09/07/03
      const3=1./(rho*sqrt( sqrt(2.)*cst_c0cub*g2_ )) ! Equilibre Production Dissipation
      do j=2,jmax-1
      do i=2,imax-1
        tridia_in(i,j,kmax+1,1)=0.
        tridia_in(i,j,kmax+1,2)=1.
        tridia_in(i,j,kmax+1,3)=0.
        tridia_in(i,j,kmax+1,4)=const3*wstress_w(i,j)
      enddo
      enddo
      endif                  !-------------------------------->        !09/07/03

      if(tke_surf.eq.1) then !-------------------------------->        !09/07/03

       do j=2,jmax-1
       do i=2,imax-1
! Dans le cas de la C.L. de Craig et banner le niveau kmax+1 ne sert
! qu'� s'assurer que le flux kz*dtken/dz entre kmax et kmax+1 ne joue
! pas. On a tken(kmax+1)=tken(kmax) que nous codons avec les lignes
! suivantes dites d'equation de surface. Le flux d'energie que Craig
! et banner estiment � 100*.ustar**3 soit 100.*(wstress/rho)**1.5
! est appliqu� au niveau kmax que nous codons avec les lignes suivantes
! dite de flux de C&B.
! Equation de surface tken(kmax+1)=tken(kmax):                !28-09-10

        tridia_in(i,j,kmax+1,1)=-1.
        tridia_in(i,j,kmax+1,2)= 1.
        tridia_in(i,j,kmax+1,3)=0.
        tridia_in(i,j,kmax+1,4)=0.

       enddo
       enddo

! Surface Boundary Condition (tke injection from waves field & wind):
! DeplacE dans anyv3d(:,:,:,id_prod) le !01-03-17
!      call turbulence_k_eps_sbc(0)  !23-05-15

      endif                  !-------------------------------->        !09/07/03

! Equation au fond:
!     const2=0.5/sqrt(ctke1*ctke2)                                !18/04/01
      const2=0.5/sqrt( sqrt(2.)*cst_c0cub*g2_ ) ! Equilibre Production Dissipation
      call turbulence_k_eps_tkebbc(const2)      !03-02-16

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! FIN.
!..............................................................................

      call tridiagonalsolver(1,2,imax-1,2,jmax-1,kmax+1)                     !30/03/06

      do k=1,kmax+1
      do j=2,jmax-1
      do i=2,imax-1

         tkea_w(i,j,k)=max(tridia_out(i,j,k)*mask_t(i,j,kmaxp1),emin)*wetmask_t(i,j)   & !27-02-21
                                                            +emin*(1.-wetmask_t(i,j))    !02-11-14

      enddo
      enddo
      enddo


!............................................!
! Lateral Open Boundary Conditions
      call obc_turbulence_tkea
!............................................!

! eps Equation:

! Advection step for epsilon field:
      call turbulence_adv_epsa

! On retablit une boucle de 2 A kmax en sachant que pour certains choix de C.L.
! le calcul est en fait inutile en k=2 et k=kmax (on recalcule les coef tridia_in
! plus loin sur ces niveaux) mais que pour certaines C.L. (Reffray 2015 pour epsilon
! notamment) la C.L. porte vraiment sur kmax+1 !16-12-16
      if(check1/=anyv3d(imax/2,jmax/2,kmax/2,id_prod)) &
      stop  ' Err anyv3d(:,:,:,id_prod) corrupted'
      if(check2/=anyv3d(imax/2,jmax/2,kmax/2,id_buoy)) &
      stop  ' Err anyv3d(:,:,:,id_buoy) corrupted'
      do k=2,kmax   ! Attention boucle obligatoirement croissante

      do j=2,jmax-1
      do i=2,imax-1

!         dz2_=0.25*( dz_t(i,j,k-1,2)+dz_t(i,j,k,2) &  ! dz_w(k,t+1/2) !01-05-11
!                    +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))
!         dz1_=0.25*( dz_t(i,j,k-1,0)+dz_t(i,j,k,0) &  ! dz_w(k,t-1/2) !01-05-11
!                    +dz_t(i,j,k-1,1)+dz_t(i,j,k,1))

! On exprime les facteurs d'echelle en passant par le tableau depth_w qui est "ajustE" dans la couche fusionnee
! en k=kmerged pour eviter la singularitE constateE. On renonce ce faisant A la definition rigoureuse de dz1 qui
! garantissait strictement la conservation de l'advection, mais ce n'est pas grave.
      dz1_=0.5*(depth_w(i,j,k+1)-depth_w(i,j,k-1)) !06-12-19

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)/cst_k_
!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)/cst_k_

!     tridia_in(i,j,k,1)=-dti_fw*0.5*(km_w(i,j,k-1)+km_w(i,j,k))/dz_t(i,j,k-1,2)*mask_t(i,j,k-1)/cst_k_
      tridia_in(i,j,k,1)=-dti_fw*0.5*(km_w(i,j,k-1)+km_w(i,j,k))/(depth_w(i,j,k)-depth_w(i,j,k-1))*mask_t(i,j,k-1)/cst_k_ !06-12-19

!     tridia_in(i,j,k,3)=-dti_fw*0.5*(km_w(i,j,k+1)+km_w(i,j,k))/dz_t(i,j,k  ,2)*mask_t(i,j,k-1)/cst_k_
      tridia_in(i,j,k,3)=-dti_fw*0.5*(km_w(i,j,k+1)+km_w(i,j,k))/(depth_w(i,j,k+1)-depth_w(i,j,k))*mask_t(i,j,k-1)/cst_k_ !06-12-19

      tridia_in(i,j,k,2)=dz1_-tridia_in(i,j,k,1)-tridia_in(i,j,k,3) &
          +dti_fw*dz1_*cst_eps_                                  &
!       *(epsn_w(i,j,k) +epsb_w(i,j,k))/(tken_w(i,j,k)+tkea_w(i,j,k))   ! s1
!       *(               epsb_w(i,j,k))/(              tken_w(i,j,k))   ! s3
        *(epsn_w(i,j,k)/tken_w(i,j,k)+epsb_w(i,j,k)/tkea_w(i,j,k))*0.5  ! s2

! flux advection vertical superieur:
!     x1=cst2_*(omega_w(i,j,k+1,1)+omega_w(i,j,k,1))     !  dt*omega/2
!     xy_t(i,j,2)=-(x1+abs(x1))*epsn_w(i,j,k)         &   ! -dt*(omega.TKE)
!                 -(x1-abs(x1))*epsn_w(i,j,k+1)

!     x0=0.5*cst1_*kh_w(i,j,k)*(rhp_t(i,j,k)-rhp_t(i,j,k-1))
      x0=0.5*anyv3d(i,j,k,id_buoy)

       tridia_in(i,j,k,4)=                                         &
       (                                                           &

!         epsn_w(i,j,k)*dz1_                                     & ! epsn_w*dz_w au temps precedent
          epsa_w(i,j,k)*dz1_   & ! epsa depuis advection preliminaire !05-01-16

! Advection:
!        +( xy_u(i+1,j  ,1)                                        &
!          -xy_u(i  ,j  ,1)                                        &
!          +xy_v(i  ,j+1,1)                                        &
!          -xy_v(i  ,j  ,1) )/dxdy_t(i,j)                          &
!          +xy_t(i  ,j  ,2)                                        &
!          -xy_t(i  ,j  ,1)                                        &

          +( anyv3d(i,j,k,id_prod)*cst_prod_                        &
!         +( anyv3d(i,j,k,0)*cst_prod_                              &
            +(x0+abs(x0))*cst_buo_p_+(x0-abs(x0))*cst_buo_m_)    &
!       *(epsn_w(i,j,k) +epsb_w(i,j,k))/(tken_w(i,j,k)+tkea_w(i,j,k))  &! m1
        *(epsn_w(i,j,k)/tken_w(i,j,k)+epsb_w(i,j,k)/tkea_w(i,j,k))*0.5 &! m2
!       *(              epsb_w(i,j,k))/(              tken_w(i,j,k))   &! m3

       )*mask_t(i,j,k-1)                                           !09/07/03

!      if(iteration3d==10) then
! Lignes utiles pour verifier les proprietEs de conservation du
! schema d'advection
!       if(i==imax/2.and.j==jmax/2.and.k==kmax/2)write(6,*)'BIDOUE EPS'
!       tridia_in(i,j,k,2)=dz2_
!       tridia_in(i,j,k,1)=0.
!       tridia_in(i,j,k,3)=0.
!       tridia_in(i,j,k,4)=epsa_w(i,j,k)*dz1_*mask_t(i,j,k-1)
!      endif

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

      enddo
      enddo

      enddo      ! fin de boucle k

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! DEBUT:
!..............................................................................

! eps=c0**3*tke**1.5/l ! l=0.4*z0

! Condition � la surface
! La loi de paroie n'est pas valide sur la paroie. La condition s'applique donc
! au premier niveau au dessous de la surface, c.a.d. k=kmax. Par convenance "graphique" la
! valeur en k=kmax est egalement appliquee en k=kmax+1
      do j=2,jmax-1
      do i=2,imax-1

! C.L. convenance graphique:
        tridia_in(i,j,kmax+1,1)=-1.
        tridia_in(i,j,kmax+1,2)=1.
        tridia_in(i,j,kmax+1,3)=0.
        tridia_in(i,j,kmax+1,4)=0.

! Vraie C.L.:
        tridia_in(i,j,kmax,1)=0.
        tridia_in(i,j,kmax,2)=1.
        tridia_in(i,j,kmax,3)=0.
!       tridia_in(i,j,k,4) computed by turbulence_k_eps_bc_xxx, next.

      enddo
      enddo
!     call turbulence_k_eps_bc_e87('srf')     ! computes tridia_in(i,j,kmax,4)
!     call turbulence_k_eps_bc_neutral('srf') ! computes tridia_in(i,j,kmax,4)
      call turbulence_k_eps_bc_tkebuo('srf')  ! computes tridia_in(i,j,kmax,4)

! Attention cette C.L. (reffray 2015) s'applique en k=kmax+1 et suppose par consequent
! de commenter les CL sur tridian_in(:,:,:,1:3) precedentes
!     call turbulence_k_eps_bc_epsclp('srf')  !16-12-16

! La loi de paroie n'est pas valide sur la paroie. La condition au fond s'applique donc
! au premier niveau au dessus du fond, c.a.d. k=2 (kmin+1). Par convenance "graphique" la
! valeur en k=2 est egalement appliquee en k=1
      do j=2,jmax-1
      do i=2,imax-1

! Parce que les couches fusionnEes n'annulent pas mask_t, la multiplication par
! mask_t des coef de la matrice ne suffit pas A produire une solution correcte
! dans les couches ecrasees. La condition de fond est donc programmee de k=1 A kmin_w !03-03-20
       do k=1,kmin_w(i,j) !03-03-20
        tridia_in(i,j,k,1)=0.
        tridia_in(i,j,k,2)=1.
        tridia_in(i,j,k,3)=-1.
        tridia_in(i,j,k,4)=0.
       enddo !03-03-20

! Vraie C.L.:
        k=kmin_w(i,j)+1
        tridia_in(i,j,k,1)=0.
        tridia_in(i,j,k,2)=1.
        tridia_in(i,j,k,3)=0.

      enddo
      enddo
!     call turbulence_k_eps_bc_e87('bot')     ! computes tridia_in(i,j,kmin_w,4)
!     call turbulence_k_eps_bc_neutral('bot') ! computes tridia_in(i,j,kmin_w,4)
      call turbulence_k_eps_bc_tkebuo('bot')  ! computes tridia_in(i,j,kmin_w,4)

!..............................................................................
! CONDITIONS LIMITES SURFACE ET FOND:
! FIN.
!..............................................................................

      call tridiagonalsolver(1,2,imax-1,2,jmax-1,kmax+1)                     !30/03/06

! Details des condition aux limites (Galperin 88 et autres...) dans:
! https://docs.google.com/document/d/1D2dCVKGz45gSCdU16bRxYb6Qjylm8pnmubIJCa84l98/edit
      x0=0.25*cst_c0cub/0.53/sqrt(2.)
      x1=epsmin ! limite par defaut !01-02-19
      do k=1,kmax+1
       k1=1
       if(k==1)k1=0      !16-09-16
       if(k==kmax+1)k1=0 !16-09-16
       do j=2,jmax-1
       do i=2,imax-1

       nn=-k1*gravoverrho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1)) & !16-09-16
                       /(depth_t(i,j,k)-depth_t(i,j,k-1))

! future limite A tester !01-02-19
! https://docs.google.com/document/d/1D2dCVKGz45gSCdU16bRxYb6Qjylm8pnmubIJCa84l98/edit
!        x1=cst_c0cub*tkea_w(i,j,k)*sqrt(tkea_w(i,j,k))  & ! c0**3 * k**(3/2) / 0.4 / Llim
!          /(0.4*min( depth_w(i,j,kmax+1)-depth_w(i,j,k) & ! Llim=min(ssh-z,z+h)
!                    ,depth_w(i,j,k)+h_w(i,j))) 


         epsa_w(i,j,k)=max(max(tridia_out(i,j,k)*mask_t(i,j,kmaxp1) & !27-02-21
                               ,x1)                                 & ! Limit !01-02-19
                           ,x0*sqrt(nn+abs(nn))*tkea_w(i,j,k))      & ! Limit Galperin 88 on lenght scale
                       *wetmask_t(i,j)                              &
                   +(1.-wetmask_t(i,j))*epsmin                     !02-11-14
!        epsa_w(i,j,k)=tridia_out(i,j,k)               

      enddo
      enddo
      enddo


!      if(iteration3d==10)stop 'yeap!'
!............................................!
! Lateral Open Boundary Conditions
      call obc_turbulence_epsa
!............................................!

! Mixing coefficients:
      call turbulence_k_eps_kz

! Bornes deplacees dans driver gaspar !13-12-19
!     km_w=min(km_w,50.)
!     kh_w=min(kh_w,50.)

! CL Surface et fond deplacees dans driver gaspar !13-12-19
!     do j=1,jmax
!     do i=1,imax
!      km_w(i,j,kmax+1)=0.
!      km_w(i,j,kmin_w(i,j))=0.
!      kh_w(i,j,kmax+1)=0.
!      kh_w(i,j,kmin_w(i,j))=0.
!     enddo
!     enddo

!     write(20,*)'----------------------------------'
!     i=30 ; j=3
!     do k=kmax+1,kmax,-1
!      write(20,*)k,real(km_w(i,j,k)),real(tkea_w(i,j,k)),real(epsa_w(i,j,k))
!     enddo
      

! Move forward
      do k=1,kmax+1
      do j=1,jmax
      do i=1,imax

       tkeb_w(i,j,k)=tken_w(i,j,k) !14-02-14
       tken_w(i,j,k)=tkea_w(i,j,k)

       epsb_w(i,j,k)=epsn_w(i,j,k) !14-02-14
       epsn_w(i,j,k)=epsa_w(i,j,k)

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




!     enddo                    !spispispispispispispispispi>

      end subroutine turbulence_k_eps

!.............................................................

      subroutine turbulence_k_eps_sbc(case_) ! Surface Boundary Condition !23-05-15
      use module_principal
      use module_parallele
      implicit none
      integer case_



       if(iwve==1) then !foc.foc.foc.foc.foc.foc>
! Case where the Wave to Ocean TKE flux is provided by the ww3 model (iwve=1):
! See discussion in Uchiyama et al 2010 (sections 3.3 & 3.4) doi:10.1016/j.ocemod.2010.04.002
        x2=dti_fw/rho
        do j=2,jmax-1
        do i=2,imax-1

!..... inverse distance verticale de distribution ..... !24-02-17 
           x0=1./(0.64*hsw_wave_t(i,j,1)) ! Basee sur hauteur des vagues deferlante

!.....
! Profil "gradient d'une fonction lineaire" (avec normalisation de l'integrale)
! Debut:
           anyv1d(kmax,1)=1.
           do k=1,kmax-1
            anyv1d(k,1)=max(0.,1.+(depth_t(i,j,k)-depth_t(i,j,kmax))*x0)
           enddo
           x3=x2*foc_wave_t(i,j,1)/(1.-anyv1d(1,1)) ! normalisation
           do k=2,kmax   ! Les extremes tke(kmax+1) et tke(kmin) ne sont pas concernes
!           tridia_in(i,j,k,4)=tridia_in(i,j,k,4)+x3*(anyv1d(k,1)-anyv1d(k-1,1))
            anyv3d(i,j,k,id_prod)= & !01-03-17
            anyv3d(i,j,k,id_prod)+x3*(anyv1d(k,1)-anyv1d(k-1,1))
           enddo
! Fin.
!.....
! Les 3 profils de Uchiyama et al 2010, Eq 53
! Debut:
!          sum1=small2 ! pour la normalisation de l'integrale
!          do k=2,kmax  
!           anyv1d(k,1)= &
!            (1.-tanh(x0*(depth_w(i,j,kmax)-depth_w(i,j,k)))**4) ! Type I
!            (1.-tanh(x0*(depth_w(i,j,kmax)-depth_w(i,j,k)))**2) ! Type II
!            cosh(x0*(depth_t(i,j,k)+h_w(i,j)))                  ! Type III 
!            sum1=sum1+anyv1d(k,1) ! pour la normalisation de l'integrale
!          enddo
!          x3=x2*foc_wave_t(i,j,1)/sum1 ! normalisation
!          do k=2,kmax   ! Les extremes tke(kmax+1) et tke(kmin) ne sont pas concernes
!           anyv3d(i,j,k,id_prod)=anyv3d(i,j,k,id_prod)+x3*anyv1d(k,1)
!          enddo
! Fin.
!.....


        enddo
        enddo
       return
       endif            !foc.foc.foc.foc.foc.foc>



       if(iwve/=1) then !cbcbcbcbcbcbcbcbcbc>

!#ifdef bidon
! If the ww3 model "foc" is not available the Wave to Ocean TKE flux is estimated
! using the Craig & Banner 1994 scheme:
        x1=dti_fw*100./(rho**1.5)    ! Coef C.B. = 100. !06-05-11
        x0=1./momentum_input_depth   ! 23-02-17 inverse distance verticale de distribution
        do j=2,jmax-1
        do i=2,imax-1

! Ces deux extremes (0 et 1) garantissent que l'integrale verticale = 1
           anyv1d(kmax,1)=1.
           anyv1d(kmin_w(i,j) ,1)=0.
           do k=kmin_w(i,j)+1,kmax-1 ! kmax-1
! Profil d'application de la condition de CB semblable a celle de la QDM !22-07-15
             anyv1d(k,1)=max(0.,1.+(depth_t(i,j,k)-depth_t(i,j,kmax))*x0)
           enddo
           do k=kmin_w(i,j)+1,kmax   ! Les extremes tke(kmax+1) et tke(kmin) ne sont pas concernes
            anyv3d(i,j,k,id_prod)=                                   & !01-03-17
            anyv3d(i,j,k,id_prod)                                    & !01-03-17
                 +(anyv1d(k,1)-anyv1d(k-1,1))*x1*wstress_w(i,j)**1.5
           enddo

!.....
! Les 3 profils de Uchiyama et al 2010, Eq 53
! Debut:
!          sum1=small2 ! pour la normalisation de l'integrale
!          do k=2,kmax  
!           anyv1d(k,1)= &
!            (1.-tanh(x0*(depth_w(i,j,kmax)-depth_w(i,j,k)))**4) ! Type I
!            (1.-tanh(x0*(depth_w(i,j,kmax)-depth_w(i,j,k)))**2) ! Type II
!            cosh(x0*(depth_t(i,j,k)+h_w(i,j)))                  ! Type III 
!            sum1=sum1+anyv1d(k,1) ! pour la normalisation de l'integrale
!          enddo
!          x3=x1*(wstress_w(i,j)**1.5)/sum1 ! normalisation
!          do k=2,kmax   ! Les extremes tke(kmax+1) et tke(kmin) ne sont pas concernes
!           anyv3d(i,j,k,id_prod)=                 & !01-03-17
!           anyv3d(i,j,k,id_prod)+x3*anyv1d(k,1)
!          enddo

        enddo
        enddo

!#endif


!      stop 'coco'
       return
       endif            !cbcbcbcbcbcbcbcbcbc>

      end subroutine turbulence_k_eps_sbc        ! Surface Boundary Condition

!...............................................................................

      subroutine turbulence_k_eps_kz
      use module_principal
      implicit none
! Voir aussi mon doc perso de travail dans drive:
! https://docs.google.com/document/d/1d-uGXJ3e154duQmkMstjU1h2L4BDJImJmw8VUZMmLDo/edit?usp=sharing

! km and kh depending on the selected turbulence stability functions:
        if(flag_tke_stab_func=='ca01') then 
          call turbulence_k_eps_ca01 ; return ! Canuto 2001
        endif
        if(flag_tke_stab_func=='kc94') then 
          call turbulence_k_eps_kc94 ; return ! kantha Clayson 1994
        endif
        if(flag_tke_stab_func=='gp88') then 
          call turbulence_k_eps_gp88 ; return ! galperin 1988
        endif

      stop 'flag_tke_stab_func not defined'
      end subroutine turbulence_k_eps_kz

!...............................................................................

      subroutine turbulence_k_eps_ca01
      use module_principal
      implicit none
      real*4,parameter ::                                              &
                t1_=-23.84,t2_=2.68,t3_=75.574,t4_=-45.48,t5_=-0.2937  &
               ,s0_=0.5168,s1_=-7.848,s2_=-0.0545,s4_=0.5412,s5_=-2.04 &
               ,s6_=0.3964

! km and kh using stability fonction Canuto "A"

       gravoverrho=grav/rho
       const3=0.5*cst_c0cub**2
       const4=sqrt(2.)*cst_c0cub
!      i=1 ; j=1 ; k=2
       do k=2,kmax ; do j=1,jmax ; do i=1,imax

        tke2overeps=0.5*( tken_w(i,j,k)**2/epsn_w(i,j,k)     &
                         +tkea_w(i,j,k)**2/epsa_w(i,j,k) )

         tkeovereps=0.5*( tken_w(i,j,k)/epsn_w(i,j,k)        &
                         +tkea_w(i,j,k)/epsa_w(i,j,k) )

       nn=-gravoverrho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1))     &
                    /(depth_t(i,j,k)-depth_t(i,j,k-1))

! gh_lim: Burchard Deleersnijder Ocean Modelling 3 (2001) Eq23
       gh=min(-const3*nn*tkeovereps**2 & ! Gh unlimited
             ,0.0233d0)                  ! Gh limit: recommended by Galperin 1988
!            ,0.0673d0)                  ! Gh limit: recommended by Burchard & Bolding 2001
!            ,0.0497d0)                  ! Gh limit: my personnal suggestion, apparently more
                                         ! consistent with the following Gm_lim (practically
                                         ! the gh value corresponding to dotted line of Fig4
                                         ! of Burchard & Bolding 2001 crossing the "0 axis")

! gm_lim: Burchard Deleersnijder Ocean Modelling 3 (2001) Eq27
       gm=min( &

          const3*(tkeovereps**2)*                               & ! Gm unlimited   
             0.5*(                                              & 
                   ((  vel_u(i  ,j,k,2)-vel_u(i  ,j,k-1,2))     &
                   /(depth_u(i  ,j,k)-depth_u(i  ,j,k-1)))**2   &
                  +((  vel_u(i+1,j,k,2)-vel_u(i+1,j,k-1,2))     &
                   /(depth_u(i+1,j,k)-depth_u(i+1,j,k-1)))**2   &
                  +((  vel_v(i,j  ,k,2)-vel_v(i,j  ,k-1,2))     &
                   /(depth_v(i,j  ,k)-depth_v(i,j  ,k-1)))**2   &
                  +((  vel_v(i,j+1,k,2)-vel_v(i,j+1,k-1,2))     & 
                   /(depth_v(i,j+1,k)-depth_v(i,j+1,k-1)))**2 ) &

              ,(1.+t1_*gh+t3_*(gh**2))/(t2_+t4_*gh )     )        ! Gm limit 

!      do gh=-0.28,0.0233,0.001 ; do x1=0,10.,0.1 ; gm=min(x1,(1.+t1_*gh+t3_*(gh**2))/(t2_+t4_*gh ))
        
        x0=const4*tke2overeps/(1.+t1_*gh+t2_*gm+t3_*(gh**2)+t4_*gh*gm+t5_*(gm**2))

        km_w(i,j,k)=x0*(s0_+s1_*gh+s2_*gm) & ! SM*sqrt(2)*(c0**3)*k**2/epsilon 
                    +kmol_m                  ! +kz molecular

        kh_w(i,j,k)=x0*(s4_+s5_*gh+s6_*gm)   ! SH*sqrt(2)*(c0**3)*k**2/epsilon

! note: Km=sqrt(2)*(c0**3)*k**2/epsilon*SM=sqrt(2*k)*l*SM=q.l*SM 
!       Kh=sqrt(2)*(c0**3)*k**2/epsilon*SH=sqrt(2*k)*l*SH=q.l*SH as in 
! Burchard Deleersnijder Ocean Modelling 3 (2001) Eq 4
 
!      write(66,*)real(gh),real(x1) &
!                   ,(s4_+s5_*gh+s6_*gm)/(1.+t1_*gh+t2_*gm+t3_*(gh**2)+t4_*gh*gm+t5_*(gm**2)) &
!                   ,(s0_+s1_*gh+s2_*gm)/(1.+t1_*gh+t2_*gm+t3_*(gh**2)+t4_*gh*gm+t5_*(gm**2))
!      enddo ; enddo

       enddo ; enddo ; enddo

!     stop 'turbulence_k_eps_ca01'
      end subroutine turbulence_k_eps_ca01

!...............................................................................

      subroutine turbulence_k_eps_gp88
      use module_principal
      implicit none
      real*4,parameter ::      &
                g0_=0.6674699  &
               ,g1_=15.62      &
               ,g2_=0.3932723  &
               ,g3_=3.085786   &
               ,g4_=34.67640   &
               ,g5_=6.1272     &
               ,g6_=0.4939277   

! km and kh using stability fonction from Galperin 88

! Voir aussi mon doc perso de travail dans drive:
! https://docs.google.com/document/d/1d-uGXJ3e154duQmkMstjU1h2L4BDJImJmw8VUZMmLDo/edit?usp=sharing

      x0=0.02   ! gh_crit for smoothing function while approaching gh_max
      x2=0.0233 ! gh_max Galperin 88 Eq29
      x3=-0.28  ! gh_min Galperin 88 Eq30
      const3=0.5*cst_c0cub**2
      const4=sqrt(2.)*cst_c0cub
      gravoverrho=grav/rho
!     i=1 ; j=1 ; k=2
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

       tke2overeps=0.5*( tken_w(i,j,k)**2/epsn_w(i,j,k)     &
                        +tkea_w(i,j,k)**2/epsa_w(i,j,k) )

        tkeovereps=0.5*( tken_w(i,j,k)/epsn_w(i,j,k)        &
                        +tkea_w(i,j,k)/epsa_w(i,j,k) )

       nn=-gravoverrho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1))     &
                    /(depth_t(i,j,k)-depth_t(i,j,k-1))


       gh=max(-const3*nn*tkeovereps**2,x3)

!      do x1=-0.48,0.5233,0.0001  ; gh=max(x1,x3) ! Pour test rapide....

       gh=gh+((-max(gh,x0)+x0)**2)/(-gh-x2+2.*x0)   ! Smoothing function while approaching gh_max

       sh= g6_/(1. -g4_*gh)                                 ! SH Galperin88 (Eq25)
       kh_w(i,j,k)=const4*tke2overeps*sh
       km_w(i,j,k)=const4*tke2overeps                    &
                 *(g2_-g3_*gh)/((1.-g4_*gh)*(1.-g5_*gh)) &  ! SM Galperin88 (Eq24)
                   +kmol_m                                  !11-07-15

!     write(66,*)x1,g6_/(1. -g4_*gh),(g2_-g3_*gh)/((1.-g4_*gh)*(1.-g5_*gh))
!     enddo

      enddo ; enddo ; enddo

!     stop 'turbulence_k_eps_gp88'
      end subroutine turbulence_k_eps_gp88

!...............................................................................

      subroutine turbulence_k_eps_kc94
      use module_principal
      implicit none
      real*4,parameter ::      &
               kc1_=0.4939277  &
              ,kc2_=-30.19200  &
              ,kc3_=-6.127200  &
              ,kc4_=0.3932723  &
              ,kc5_=17.07336    

! Voir aussi mon doc perso de travail dans drive:
! https://docs.google.com/document/d/1d-uGXJ3e154duQmkMstjU1h2L4BDJImJmw8VUZMmLDo/edit?usp=sharing

! Mixing coefficients:
! Km=sqrt(2tke).L.Sm=sqrt(2)sqrt(tke)(c0**3*tke**1.5/eps).Sm=sqrt(2)*c0**3*tke**2/eps.Sm

!     la longueur de melange L � partir de tke et eps: L=(c**3)*(tke**1.5)/eps
!     Gh=L**2/(2Tke)*(grav/rho0)*drho/dz=0.5*(c**3)**2*(g/rho0)*(tke/eps)**2*drho/dz

      x0=0.02   ! gh_crit for smoothing function while approaching gh_max
      x2=0.029  ! gh_max KC94 Eq44
      x3=-0.28  ! gh_min KC94 Eq45
      const3=0.5*cst_c0cub**2
      const4=sqrt(2.)*cst_c0cub
      gravoverrho=grav/rho

!     i=1 ; j=1 ; k=2
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

       tke2overeps=0.5*( tken_w(i,j,k)**2/epsn_w(i,j,k)     &
                        +tkea_w(i,j,k)**2/epsa_w(i,j,k) )

        tkeovereps=0.5*( tken_w(i,j,k)/epsn_w(i,j,k)        &
                        +tkea_w(i,j,k)/epsa_w(i,j,k) )

       nn=-gravoverrho*(rhp_t(i,j,k)  -rhp_t(i,j,k-1))     &
                    /(depth_t(i,j,k)-depth_t(i,j,k-1))


       gh=max(-const3*nn*tkeovereps**2,x3)

!      do x1=-0.48,0.5233,0.0001  ; gh=max(x1,x3) ! pour test rapide

       gh=gh+((-max(gh,x0)+x0)**2)/(-gh-x2+2.*x0)   ! Smoothing function while approaching gh_max

       sh=kc1_/(1.+kc2_*gh)                                ! SH KC94 (Eq28)
       kh_w(i,j,k)=const4*tke2overeps*sh
       km_w(i,j,k)=const4*tke2overeps                    &
                         *(kc4_+kc5_*sh*gh)/(1.+kc3_*gh) & ! SM KC94 (Eq29)
                   +kmol_m                                 !11-07-15

!     write(66,*)x1,kc1_/(1.+kc2_*gh),(kc4_+kc5_*sh*gh)/(1.+kc3_*gh)
!     enddo


!      if(i==85.and.k==2.and.j<=2)write(6,*)tken_w(i,j,k)

      enddo ; enddo ; enddo

!     stop 'turbulence_k_eps_kc94'
      end subroutine turbulence_k_eps_kc94

!.....................................................................

      subroutine turbulence_k_eps_bc_e87(txt_) !19-08-15
      use module_principal
      implicit none
      character*3 txt_


      if(txt_=='srf') then !sssssssss>
       do j=2,jmax-1
       do i=2,imax-1

        x0=-anyv3d(i,j,kmax,2)/(small2+anyv3d(i,j,kmax,1))  ! Ri*Kh/Km

        if(x0<0.16) then !est87est87est87est87>

         tridia_in(i,j,kmax,4)=cst_c0cub*tkea_w(i,j,kmax)**1.5    &
         /(0.4*(dz_t(i,j,kmax,1)+z0s)                                 & ! Estournel et al 1987 (Eq9)
              *(1.-5.*0.5*(x0+abs(x0))) )                               ! http://dx.doi.org/10.1007/BF00121874

        else             !est87est87est87est87>

         tridia_in(i,j,kmax,4)=cst_c0cub*tkea_w(i,j,kmax)**1.5    &
         /(0.4*(dz_t(i,j,kmax,1)+z0s)                                 & !Estournel et al 1987 (Eq13)
              *(1./(1.+41.*x0))**0.8 )

        endif            !est87est87est87est87>

        enddo
        enddo
      return
      endif                !sssssssss>

      if(txt_=='bot') then !bbbbbbbbb>
       do j=2,jmax-1
       do i=2,imax-1

        k=kmin_w(i,j)+1

        x0=-anyv3d(i,j,k,2)/(small2+anyv3d(i,j,k,1))          ! Ri*Kh/Km

        if(x0<0.16) then !est87est87est87est87>

         tridia_in(i,j,k,4)=cst_c0cub*tkea_w(i,j,k)**1.5     &
        /(0.4*( dz_t(i,j,k-1,1)+z0_w(i,j) )                  &
              *(1.-5.*0.5*(x0+abs(x0))) )                               ! http://dx.doi.org/10.1007/BF00121874

        else             !est87est87est87est87>

         tridia_in(i,j,k,4)=cst_c0cub*tkea_w(i,j,k)**1.5     &
        /(0.4*( dz_t(i,j,k-1,1)+z0_w(i,j) )                  &
              *(1./(1.+41.*x0))**0.8 )

        endif            !est87est87est87est87>

        enddo
        enddo
      return
      endif                !bbbbbbbbb>

      stop 'Err 1072 txt_ not recognised in turbulence_k_eps_bc_e87'
      end subroutine turbulence_k_eps_bc_e87

!.....................................................................

      subroutine turbulence_k_eps_bc_neutral(txt_) !19-08-15
      use module_principal
      implicit none
      character*3 txt_


      if(txt_=='srf') then !sssssssss>
       do j=2,jmax-1
       do i=2,imax-1


         tridia_in(i,j,kmax,4)=cst_c0cub*tkea_w(i,j,kmax)**1.5    &
                                    /(0.4*(dz_t(i,j,kmax,1)+z0s)  )

        enddo
        enddo
      return
      endif                !sssssssss>

      if(txt_=='bot') then !bbbbbbbbb>
       do j=2,jmax-1
       do i=2,imax-1

        k=kmin_w(i,j)+1

         tridia_in(i,j,k,4)=cst_c0cub*tkea_w(i,j,k)**1.5     &
        /(0.4*( dz_t(i,j,k-1,1)+z0_w(i,j) )  )

        enddo
        enddo
      return
      endif                !bbbbbbbbb>

      stop 'Err 1072 txt_ not recognised in turbulence_k_eps_bc_neutral'
      end subroutine turbulence_k_eps_bc_neutral

!.....................................................................

      subroutine turbulence_k_eps_bc_tkebuo(txt_) !19-08-15
      use module_principal
      implicit none
      character*3 txt_

! https://docs.google.com/document/d/1D2dCVKGz45gSCdU16bRxYb6Qjylm8pnmubIJCa84l98/edit
! x1=length scale Eq(9) Reffray et al 2015
! Reference:
! Reffray et al 2015: Geosci. Model Dev., 8, 69–86, 2015
! www.geosci-model-dev.net/8/69/2015/
! doi:10.5194/gmd-8-69-2015

      x0=1./sqrt(grav/rho)

      if(txt_=='srf') then !sssssssss>
      k=kmax
      do j=2,jmax-1
      do i=2,imax-1

!#ifdef stokes
! Note: 0.4*1.6=0.64
!      z0s=1.6*hsw_wave_t(i,j,1) ! these Nicolas Rascle page 43 & Eq 2.12  !06-05-11!02-03-17
!#endif

      x1=                                                             &

! methode 1
!     min(min(depth_w(i,j,k)     -depth_w(i,j,kmin_w(i,j))+z0_w(i,j), &
!             depth_w(i,j,kmax+1)-depth_w(i,j,k)          +z0s)     , &

! La methode 2 fait l'hypothese que 0.5*(depth_t(i,j,k)+depth_t(i,j,k-1)) est plus
! representative que depth_w(i,j,k) compte tenu du facteur d'echelle dz des variables
! turbulente, cad depth_t(:,:,k)-depth_t(:,:,k-1). Un petit test (qui demanderait A Etre
! confirmE) donne l'impression que la methode 2 donne de meilleurs resultats que methodes 1 et 3
! methode 2 !19-12-16
!     min(min(0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))  -depth_w(i,j,kmin_w(i,j))+z0_w(i,j), &
!             depth_w(i,j,kmax+1)-0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))+z0s)     , &

! Sans prise en compte de z0s (meilleur comportement constatE) !10-03-17
      min(min(0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))  -depth_w(i,j,kmin_w(i,j)), &
                    depth_w(i,j,kmax+1)-0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))), &

! methode 3
!     min(min(0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))  -depth_t(i,j,kmin_w(i,j))+z0_w(i,j), &
!             depth_t(i,j,kmax  )-0.5*(depth_t(i,j,k)+depth_t(i,j,k-1))+z0s)     , &

      x0*sqrt(-2.*(tken_w(i,j,k)+tkea_w(i,j,k))                       &
                               *   (depth_t(i,j,k)-depth_t(i,j,k-1))  &
                               /(     rhp_t(i,j,k)  -rhp_t(i,j,k-1)   &
                                 -abs(rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
                                 -small1 )))


      tridia_in(i,j,kmax,4)=cst_c0cub*(tkea_w(i,j,kmax)**1.5)/(0.41*x1) !coherent avec k-epsilon
!     tridia_in(i,j,kmax,4)=ctke2*(tkea_w(i,j,kmax)**1.5)/x1            !coherent avec Gaspar


      enddo
      enddo

      return
      endif                !sssssssss>

      if(txt_=='bot') then !bbbbbbbbb>
      do j=2,jmax-1
      do i=2,imax-1

      k=kmin_w(i,j)+1
!     do k=kmin_w(i,j)+1,kmerged_t(i,j)+1 !12-01-19
!     do k=kmin_w(i,j)+1,kmerged_t(i,j)+5 !12-01-19
!     do k=kmin_w(i,j)+1,kmax-1

      x1=                                                             &
      min(min(depth_w(i,j,k)     -depth_w(i,j,kmin_w(i,j))          , &
              depth_w(i,j,kmax+1)-depth_w(i,j,k)              )     , &
!     min(min(depth_w(i,j,k)     -depth_w(i,j,kmin_w(i,j))+z0_w(i,j), &
!             depth_w(i,j,kmax+1)-depth_w(i,j,k)          +z0s)     , &

      x0*sqrt(-2.*(tken_w(i,j,k)+tkea_w(i,j,k))                       &
                               *   (depth_t(i,j,k)-depth_t(i,j,k-1))  &
                               /(     rhp_t(i,j,k)  -rhp_t(i,j,k-1)   &
                                 -abs(rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
                                 -small1 )))

      tridia_in(i,j,k,4)=cst_c0cub*(tkea_w(i,j,k)**1.5)/(0.4*x1) !12-01-19 coherent avec k-eps
!     tridia_in(i,j,k,4)=    ctke2*(tkea_w(i,j,k)**1.5)/x1       !         coherent avec gaspar
!     enddo

      enddo
      enddo


!#ifdef bidon
      if(flag_merged_levels==1) then !m°v°m> !11-12-19
! Si le niveau zw(kmerged) est collE au niveau kmin alors il faut appliquer le schema L
! au niveau kmerged+1. La variable rap introduit la progressivite necessaire A la continuite 
! des champs. Note: si kmerged=kmin rap=1 donc coef inchangEs
       do j=2,jmax-1 ; do i=2,imax-1

       k=kmerged_t(i,j)+1
       rap=min((dz_t(i,j,kmin_w(i,j),1)/dz_t(i,j,kmerged_t(i,j),1))**2,1.) !03-11-21

      x1=min( & !ooo>
          min(depth_w(i,j,k)     -depth_w(i,j,kmin_w(i,j)),depth_w(i,j,kmax+1)-depth_w(i,j,k)) &
          ,x0*sqrt(-2.*(tken_w(i,j,k)+tkea_w(i,j,k))                  &
                               *   (depth_t(i,j,k)-depth_t(i,j,k-1))  &
                               /(     rhp_t(i,j,k)  -rhp_t(i,j,k-1)   &
                                 -abs(rhp_t(i,j,k)  -rhp_t(i,j,k-1))  &
                                 -small1 ))                           &
            )   !ooo>

      tridia_in(i,j,k,1)=rap *tridia_in(i,j,k,1)
      tridia_in(i,j,k,3)=rap *tridia_in(i,j,k,3)
      tridia_in(i,j,k,2)=rap *tridia_in(i,j,k,2)+(1.-rap)*0.5*(depth_w(i,j,k+1)-depth_w(i,j,k-1))
      tridia_in(i,j,k,4)=rap *tridia_in(i,j,k,4)  &
                    +(1.-rap)*(cst_c0cub*(tkea_w(i,j,k)**1.5)/(0.4*x1))*0.5*(depth_w(i,j,k+1)-depth_w(i,j,k-1))

       enddo ; enddo

      endif                                             !m°v°m> !06-12-19
!#endif


      return
      endif                !bbbbbbbbb>


      stop 'Err 1072 txt_ not recognised in turbulence_k_eps_bc_tkebuo'
      end subroutine turbulence_k_eps_bc_tkebuo

!.................................................................

! subroutine transferee dans turbulence_adv.F90 !06-08-16
!     subroutine turbulence_adv_tkea  ! 05-01-16
!     end subroutine turbulence_adv_tkea

!.................................................................

! subroutine transferee dans turbulence_adv.F90 !06-08-16
!     subroutine turbulence_adv_epsa  ! 05-01-16
!     end subroutine turbulence_adv_epsa

!.................................................................

! subroutine transferee dans turbulence_adv.F90 !06-08-16
!     subroutine turbulence_adv_tken  ! 05-01-16
!     end subroutine turbulence_adv_tken

!.........................................................................................

      subroutine turbulence_k_eps_tkebbc(coefprodi_) ! tke bottom boundary condition
      use module_principal
      implicit none
      double precision coefprodi_

! Parce que les couches fusionnEes n'annulent pas mask_t, la multiplication par
! mask_t des coef de la matrice ne suffit pas A produire une solution correcte
! dans les couches ecrasees. La condition de fond est donc programmee de k=1 A
! kmin_w !03-03-20

! tke bottom boundary condition assuming production/dissipation balance
      do j=2,jmax-1 ; do i=2,imax-1 

!     do k=kmin_w(i,j),kmerged_t(i,j) !12-01-19
!        k=kmin_w(i,j)                !06-12-19
! Parce que les couches fusionnEes n'annulent pas mask_t, la multiplication par
! mask_t des coef de la matrice ne suffit pas A produire une solution correcte
! dans les couches ecrasees. La condition de fond est donc programmee de k=1 A
! kmin_w !03-03-20
      do k=1,kmin_w(i,j) !03-03-20

        tridia_in(i,j,k,1)=0.
        tridia_in(i,j,k,2)=1.
        tridia_in(i,j,k,3)=0.

      enddo  !03-03-20

      enddo ; enddo

      if(flag_linearfric==1) then !fricfricfric>

       do j=2,jmax-1 ; do i=2,imax-1

          k=kmin_w(i,j)                !06-12-19

        tridia_in(i,j,k,4)=                             & 
       coefprodi_*coef_linearfric*sqrt(                 & !03-02-16
                   velbot_u(i+1,j  )**2  & 
                  +velbot_u(i  ,j  )**2  & 
                  +velbot_v(i  ,j+1)**2  & 
                  +velbot_v(i  ,j  )**2)


       enddo ; enddo

      else                        !fricfricfric>

      do j=2,jmax-1 ; do i=2,imax-1

         k=kmin_w(i,j)                !06-12-19

        tridia_in(i,j,k,4)=                                          &
       coefprodi_*cdb_t(i,j)*(                                       & !09/07/03
                   velbot_u(i+1,j  )**2               & !09/07/03
                  +velbot_u(i  ,j  )**2               & !09/07/03
                  +velbot_v(i  ,j+1)**2               & !09/07/03
                  +velbot_v(i  ,j  )**2               & !09/07/03
                         )

      enddo ; enddo

      endif                       !fricfricfric>

      if(flag_merged_levels==1) then !m°v°m> !06-12-19

       do j=2,jmax-1 ; do i=2,imax-1
          k2=kmerged_t(i,j)
          k1=   kmin_w(i,j)
          rap=min((dz_t(i,j,k1,1)/dz_t(i,j,k2,1))**2,1.) !03-11-21
         tridia_in(i,j,k2,1)=rap *tridia_in(i,j,k2,1)  &
                         +(1-rap)*tridia_in(i,j,k1,1)*0.5*(depth_w(i,j,k2+1)-depth_w(i,j,k2-1))
         tridia_in(i,j,k2,2)=rap *tridia_in(i,j,k2,2)  &
                         +(1-rap)*tridia_in(i,j,k1,2)*0.5*(depth_w(i,j,k2+1)-depth_w(i,j,k2-1))
         tridia_in(i,j,k2,3)=rap *tridia_in(i,j,k2,3)  &
                         +(1-rap)*tridia_in(i,j,k1,3)*0.5*(depth_w(i,j,k2+1)-depth_w(i,j,k2-1))
         tridia_in(i,j,k2,4)=rap *tridia_in(i,j,k2,4)  &
                         +(1-rap)*tridia_in(i,j,k1,4)*0.5*(depth_w(i,j,k2+1)-depth_w(i,j,k2-1))
       enddo ; enddo

! Dans les couches ecrasees on reporte la condition prodi du niveau kmin_w !03-03-20
! sachant que les coef 1,2,3 ont ete plus haut determines
       do j=2,jmax-1 ; do i=2,imax-1
        do k=1,kmin_w(i,j)-1 !03-03-20
         tridia_in(i,j,k,4)=tridia_in(i,j,kmin_w(i,j),4)
        enddo
       enddo ; enddo

      endif                                             !m°v°m> !06-12-19


      end subroutine turbulence_k_eps_tkebbc

!...............................................................................

      subroutine turbulence_k_eps_bc_epsclp(txt_) !16-12-16
      use module_principal
      implicit none
      character*3 txt_

! Condition de surface pour epsilon
! Reference:
! Reffray et al 2015: Geosci. Model Dev., 8, 69–86, 2015
! www.geosci-model-dev.net/8/69/2015/
! doi:10.5194/gmd-8-69-2015
! La base consiste A considerer la formule 24c de l'article de Reffray et al 2015
! On suppose un equilibre entre epsilon et le melange turbulent qu'on approxime
! a la prise de l'input de TKE en surface (condition de CB) repartie sur une
! couche epaisse de z0s. Concretement epsilon=(Fluxsurface-Fluxsousjacent)/z0s
! avec Fluxsousjacent negligeable par rapport A Fluxsurface=flux de CB94
! Soit epsilon=(FluxCB94)/z0s
! Mais attention aux regimes "sans vagues" qui conduirait epsilon A etre nul ou disons
! trop petit (entrainant TKE immense...). On envisage donc un regime hybride
! basE sur equation 8 de Reffray 2015 epsilon=c0**3*tke**1.5/L avec L=karman*z0s (Eq18b)
! A noter qu'on utilise tkemin et non pas tke(i,j,k) pour eviter d'entrainer des oscillations
! liee A la grille "fusion des couche" 


      if(txt_=='srf') then !sssssssss>

      do j=2,jmax-1
      do i=2,imax-1


       tridia_in(i,j,kmax+1,1)=0.
       tridia_in(i,j,kmax+1,2)=1.
       tridia_in(i,j,kmax+1,3)=0.

       tridia_in(i,j,kmax+1,4)=max( & !pmxpmx> !18-09-16
                                  (100./z0s)*(wstress_w(i,j)/rho)**1.5 &!Eq 24c
                                 ,cst_c0cub*(emin**1.5)/(0.41*z0s)     &!Eq8
                                   )  !pmxpmx> !18-09-16

!      tridia_in(i,j,kmax+1,4)=(100./z0s)*(wstress_w(i,j)/rho)**1.5  !Eq 24c

      enddo
      enddo

      return
      endif                !sssssssss>

      if(txt_=='bot') then !bbbbbbbbb>
        stop 'turbulence_k_eps_bc_epsclp'
      return
      endif                !bbbbbbbbb>


      stop 'Err 1072 txt_ not recognised in turbulence_k_eps_bc_epsclp'
      end subroutine turbulence_k_eps_bc_epsclp

!.................................................................
