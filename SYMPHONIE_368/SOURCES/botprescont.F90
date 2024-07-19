      subroutine botprescont(txt_loc)
!______________________________________________________________________
! S model
! release S.26 - last update: 10-09-16 
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      implicit none
      double precision mean_ssh_loc,lagrange_coef_loc,ddz_loc    &
       ,rhp_avr_loc,restore_loc
      character*10 txt_loc
#ifdef synopsis
       subroutinetitle='botprescont'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef bidon

!...............................................................................
! Version date      Description des modifications
!         13/01/08: mise en service
!         18/01/08: adaptation a un cas avec maree
!         21/01/08: on utilise une echelle de temps lue dans notebook_spongelayer
!         22/02/08: securisation de l'ecart par rapport a l'ogcm
!         06/03/08: Assure que DT et DS sont de signes opposes
!         23/02/09: Modification de seuillage dans le calcul du coef de Lagrange, car
!                   lorsque 2 champs ogcm consecutifs etaient identiques ou quasi
!                   identiques on tombait sur une division par zero.
!         05/03/09: renforcement du point precedent si rap_obc<0 on sort de la
!                   routine
!         06/03/09: Annulation du point precedent mais borne max pour la fonction
!                   d'amplification de correction.
!         05-06-09: Parallelisation
! 2009.3  01-10-09: Utilisation des nouveaux facteurs d'echelle verticale
!         05-10-09: ajout d'un "ifdef parallele"
!         10-10-09: tem et sal remplacent thz et shz
! 2010.11 16-07-10  temobc & salobc renommes temobc & salobc
! 2010.12 20-09-10  Possibilite de calcul en simple precision
! 2010.22 22-04-11  - relax_ts est desormais separe de sponge_t, rendant la
!                   division de sponge_t par relax_ts inutile
!         06-05-11  Ajout d'un schema Lagrangien
!         11-07-14  rap_obc devient timeweightobc
!         10-09-16  routine deactivee
!...............................................................................

      if(nest_onoff_in.eq.1) then ! Stop debug
      write(6,*)'l''imbrication symphonie-symphonie n''est pas'
      write(6,*)'prête pour le controle de la pression de fond'
      write(6,*)'car il faut la moyenne de la sse de forcage.'
      stop ' donc dans botprescont.f'
      endif                       ! Stop debug

!-------------------------------------------------------
! Ici on calcule la moyenne de la sse, la moyenne
! de la sse de MOM, puis la difference entre les
! deux. Pourquoi? Parce qu'on ne cherche pas a corriger
! la moyenne de la sse par cette methode. L'ecart moyenne
! (CONST9) n'est donc pas pris dans les calculs qui
! suivront. Plus loin, La soustraction de CONST9
! s'explique comme cela.
! Attention que ce calcul de moyenne n'est pas forcement
! adapte aux simu avec la maree. Il faut aussi que le
! forcage soit dispo sur toute la grille, ce qui n'est
! pas toujours le cas.... Prudence donc.
!     x2=    rap_obc
!     x0=(1.-rap_obc)
      x2=timeweightobc(ssh_id) !11-07-14
      x0=1.-x2
      sum1=0.
      sum2=0.
      do j=1,jmax
      do i=1,imax
       if(  mask_t(i,j,kmax+1).eq.1) then !------>

! Poid=H*masque_parallelisation:
      x10=h_w(i,j)                                                      &
        *mask_i_w(i)*mask_j_w(j) !#MPI ne pas sommer zone fantome 05-06-09

! On controle avant tout la zone profonde (les plateaux reagissent trop
! fortement, notamment au vent). Donc la moyenne est ponderee par H
        sum1=sum1+x10                                         !05-06-09
        sum2=sum2+x10*(ssh_int_w(i,j,1)                                 &
                    -x2*sshobc_w(i,j,2)                                 &
                    -x0*sshobc_w(i,j,0))
       endif                          !------>
      enddo
      enddo

#ifdef parallele
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI 05-06-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI 05-06-09
                         mpi_sum,par%comm2d,ierr)
#else
      sum1glb=sum1
      sum2glb=sum2
#endif
      mean_ssh_loc=sum2glb/sum1glb                                   & !#MPI 05-06-09
            -sshmeatide(2) ! moins contribution de la maree            !18/01/08
                           ! a la moyenne ponderee par H

#ifdef bidouille
      stop 'L''etiquette Bidouille ne sert plus a rien!'
#endif


!..................................................................................
! NEXT STEP: A VERTICAL LAGRANGIAN APPROACH:
!..................................................................................
      if(trim(txt_loc)=='lagrangian') then !lllllllllllllllllllllllllll> !06-05-11


      const1=rho*relax_bpc*dti_fw          !

      do loop1=1,4                         ! Boundary loop
!     do loop1=5,5 ! BIDOUILLE POUR PASSER SUR TOUT LE DOMAINE
      do j=spo_j1_t(loop1),spo_j2_t(loop1) ! j loop
      do i=spo_i1_t(loop1),spo_i2_t(loop1) ! i loop

       if(mask_t(i,j,kmaxp1)==1) then !>>>>>>>>>>>>>>>>

! dz adjusted (dza) = dz+LagrangeCoef*dz*(rhp-<rhp>)
! <rhp>=sum(rhp*dz)/sum(dz)
        sum1=small1
        sum2=0.
        do k=kmin_w(i,j),kmax
         sum1=sum1+dz_t(i,j,k,after)
         sum2=sum2+dz_t(i,j,k,after)*rhp_t(i,j,k)
        enddo
        rhp_avr_loc=sum2/sum1    ! = <rhp>

!      if(par%rank==0)then
!        if(j==1) then
!         if(i==imax/2) then
!          write(*,*)'rhp_avr_loc',rhp_avr_loc
!          stop 'feu'
!         endif
!        endif
!       endif

        sum1=small1
        do k=kmin_w(i,j),kmax
         sum1=sum1+dz_t(i,j,k,after)*(rhp_t(i,j,k)-rhp_avr_loc)*rhp_t(i,j,k)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIDOUILLE:
!       i1=imax/2 ; j1=jmax/2 ; x2=0.01
!       sshobc_w(i,j,1)=    &
!            max(0.,x2*(1.-sqrt(real(i-i1)**2+real(j-j1)**2)/30.))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        lagrange_coef_loc=const1*(   ssh_int_w(i,j,1)        &
                                  -sshtotide_w(i,j,1)        & ! la maree est otee
                                     -sshobc_w(i,j,1)        &
                                       -mean_ssh_loc  )/sum1

! Note: sponge_t(i,j,1)=1 at open boundaries (ob) and =0 far from the ob
        x1=1.-sponge_t(i,j,1)*(1.-cellboxfactor1)
        x2=1.-sponge_t(i,j,1)*(1.-cellboxfactor2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIDOUILLE:
!       x1=cellboxfactor1
!       x2=cellboxfactor2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       sum0=0.
        restore_loc=0.
        do k=kmin_w(i,j),kmax

! dz of reference:
         x0=hz_w(i,j,after)*dsig_t(i,j,k)

! dznew=dz*(1.+lagrange_coef_loc*(rhp-<rhp>))
         anyv3d(i,j,k,1)=dz_t(i,j,k,after)                       &
        *( 1.+lagrange_coef_loc*(rhp_t(i,j,k)-rhp_avr_loc) )

         if( anyv3d(i,j,k,1) > x2*x0 )                           &
             restore_loc=max(restore_loc,                        &
            (anyv3d(i,j,k,1)-x2*x0)/(anyv3d(i,j,k,1)-x0) )

         if( anyv3d(i,j,k,1) < x1*x0 )                           &
             restore_loc=max(restore_loc,                        &
            (anyv3d(i,j,k,1)-x1*x0)/(anyv3d(i,j,k,1)-x0) )

!       sum0=sum0+lagrange_coef_loc*(rhp_t(i,j,k)-rhp_avr_loc)*dz_t(i,j,k,after)*rhp_t(i,j,k)

!      if(par%rank==0)then
!        if(j==1) then
!         if(i==imax/2) then
!          write(*,*)k,( 1.+lagrange_coef_loc*(rhp_t(i,j,k)-rhp_avr_loc) ),restore_loc
!         endif
!        endif
!       endif

        enddo

!      if(par%rank==0)then
!        if(j==1) then
!         if(i==imax/2) then
!          write(*,*)'sum0',sum0,rho*1.
!         endif
!        endif
!       endif

        if(restore_loc/=0.) then !------->
!       sum0=0.
!       sum1=0.
         x0=1.-restore_loc
         do k=kmin_w(i,j),kmax
             anyv3d(i,j,k,1)=                                         &
          x0*anyv3d(i,j,k,1)+restore_loc*hz_w(i,j,after)*dsig_t(i,j,k)

!         sum0=sum0+anyv3d(i,j,k,1)
!      if(par%rank==0)then
!        if(j==1) then
!         if(i==imax/2) then
!          write(*,*)'R',k,anyv3d(i,j,k,1)/hz_w(i,j,after)/dsig_t(i,j,k)
!         endif
!        endif
!       endif
         enddo
        endif                    !------->
!      if(par%rank==0)then
!        if(j==1) then
!         if(i==imax/2) then
!          write(*,*)'sum0 et hz',sum0,hz_w(i,j,after)
!         endif
!        endif
!       endif


         do k=kmin_w(i,j),kmax
          ddz_loc=anyv3d(i,j,k,1)-dz_t(i,j,k,after)
          dz_t(i,j,k,after) =dz_t(i,j,k,after) +ddz_loc
          dz_t(i,j,k,now)   =dz_t(i,j,k,now)   +ddz_loc
          dz_t(i,j,k,before)=dz_t(i,j,k,before)+ddz_loc
         enddo

       endif                          !>>>>>>>>>>>>>>>>
      enddo ! Boundary loop
      enddo ! j
      enddo ! i

      return
      endif                                !lllllllllllllllllllllllllll>

!..................................................................................
! NEXT STEP: A VERTICAL EULERIAN APPROACH:
!..................................................................................
      if(trim(txt_loc)=='eulerian') then   !eeeeeeeeeeeeeeeeeeeeeeeeeee>
!...............................................................................
! Controle de la pression de fond.
!...............................................................................

! Soit dzet ecart entre sse et sse de forçage: dzet=sse-sseF
! Imaginons une perturbation de densite proportionnelle a
! dzet via une constante alpha (que nous allons chercher
! a determiner) et une fonction arbitraire dependante de z
! uniquement: drho=alpha*f(z)*dzet
! Soit F la primitive de f
! La perturbation de la pression de fond est approximativement:
! dp= g*rho0*dzet + g*alpha*(F(0)-F(-h))*dzet
! Imaginons un instant que les perturbations de densite et de
! surface se sont ajustees l'une a l'autre de sorte que la pression
! de fond soit inchangee (dp=0) alors on peut exprimer alpha
! en fonction de dzet et F: alpha=-rho0 / (F(0)-F(-h))
! A partir de la notre raisonnement est le suivant. On constate
! un ecart entre le model_ et le forcage: dzet=sse-sseF. On modifie
! alors la densite aux conditions aux limites du model_ en comptant
! sur une reponse du model_ de type "pression de fond inchangee" pour
! provoquer un deplacement de la sse (pour rapprocher celle ci
! de la sse de forcage). En pratique on cherche a ce que sse se deplace
! de -dzet, autrement dit la valeur alpha=rho0/( F(0)-F(-h))
! est celle que nous appliquons aux frontieres.
! La correction de densite est donc:
! delta_rho=rho0/(F(0)-F(-h))*f(z)*(sse-sseF)
! En resumer, quand on constante que la sse est "trop" haute on injecte
! de l'eau plus dense pour faire baisser le niveau et inversement.

! Comment choisir f(z)? On va considerer que f(z) est la forme
! "typique" des variations de densite sur la verticale et qu'on s'inspire
! du forcage. On decide que f(z) sera donnee par la variation
! des champs de forcages (ici model_ mom) entre 2 echeances successive
! (ici 2 semaines).

! Comment reporter la correction sur les tableaux de forcages?
! Le delta_rho ci dessus est un increment de densite par rapport
! a la densite du model_, rho. La densite de forcage, aux frontieres
! ouvertes, est donc:
! rho_obc=rho+rho0/(F(0)-F(-h))*f(z)*(sse-sseF)

! Enfin une fois que la densite est corrigee il faut traduire cela
! en perturbation de temperature et de salinite. Ici on fait un
! truc sommaire qui consiste a dire que la temperature varie "5 fois
! plus" que la salinite. dS=-dT/5 le signe "moins" permettant a la
! temperature et a la salinity de "tirer dans le même sens":
! pour augmenter la densite on "augmente S" et on "diminue T"
! et inversement. Soient Tf et Sf les temperature et salinite dans le
! model_ de grande echelle (mom) et dT et dS les
! ajustements afin de retrouver rho_obc. Soient (a,b,c) les coef
! de l'equation d'etat lineaire. On a donc:
!     rho_obc=a(Tf+dT)+b(Sf+dS)+c avec dS=-dT/5
! On en tire donc dT=( rho_obc-aTf-bSf-c)/(a-b/5)
! et dans la foulee dS=-dT/5

! Comment on evite les "exces" de correction?
! Bien sûr il faut faire attention  a ne pas obtenir des
! corrections (dT,dS) irrealistes. On maintient dT dans une
! barre d'erreur qui nous est donnee par une estimation de
! la variabilite dans
! le model_ de grande echelle (mom). Cette derniere nous est
! simplement donnee par l'ecart en temperature entre 2 echeances
! du model_
! mom (2 semaines): valeur absolue de dt inferieure a
! valeur absolue de ( Tf(t+2semaines)-Tf(t) )*const5 où
! const5 est un bouton de calibration (on a pris const5=1)
! pour commencer.


! Borne max pour la fonction d'amplification de correction:
      const3=2.

! coef equation d'etat: rho=c1*T+c2*S+c7
      const1=-rho*alp_t
      const2= rho*alp_s
      const7=-const1*t0-const2*s0

! rhonouveau=rho+drho*exp(z/c3)
      const6=const1+const4*const2

! Coeficient de calibration de la barre d'erreur:
!cccccCONST5=3. ! ecart tolere (fois le delta du champs de forcage entre 2 echeances successives)

! La division par RELAX_TS permet de ramener le tableau sponge_z a 1
! (sur frontiere). La division suivante permet de donner une echelle
! de temps au transitoire. Attention on divise par HMAX car ensuite
! on multiplie par H. C'est un moyen de donner plus d'importance a la zone
! profonde:
!cccccCONST8=DTI_FW/(HMAX*RELAX_TS*86400.) ! Temps de transition: 1  jour
!     const8=relax_bpc*dti_fw/(hmax*relax_ts) ! Temps de transition: 1  jour
      const8=relax_bpc*dti_fw/hmax         !22-04-11

      do loop1=1,4                         ! boucle sur 4 frontieres
      do j=spo_j1_t(loop1),spo_j2_t(loop1) ! debut boucle j
      do i=spo_i1_t(loop1),spo_i2_t(loop1) ! debut boucle i

      if(  mask_t(i,j,kmax+1).eq.1) then !>>>>>>>>>>>>>>>>

       sum1=0.
       do k=kmin_w(i,j),kmax ! debut de boucle n°2 sur k

! Pour eviter qu'il y ait des tendances contradictoire entre le fond
! et la surface (au fond ca augmenterait quand en surface ca diminuerait)
! on impose un "sens unique":
!      X3=SIGN(UN, CONST1*(TOBC_Z(I,J,K,2)-TOBC_Z(I,J,K,0))
!    &            +CONST2*(SOBC_Z(I,J,K,2)-SOBC_Z(I,J,K,0)))
!    &        *HZ_Z(I,J,1)*DSIG_Z(I,J,K)
!      x3=hz_z(i,j,1)*dsig_c(i,j,k)                                    !06/03/08
!      x3=dz_c(i,j,k,1)                                                !01-10-09

! ANYV3D est la fonction f(z) fois dz
!      anyv3d(i,j,k,1)= x3*abs(temobc_c(i,j,k,2)-temobc_c(i,j,k,0))        !06/03/08
!    &              *(1.-exp( depth_c(i,j,k)/75.)) ! on n'aime pas trop la variabilite de surface...
!      anyv3d(i,j,k,2)=-x3*abs(salobc_c(i,j,k,2)-salobc_c(i,j,k,0))        !06/03/08
!    &              *(1.-exp( depth_c(i,j,k)/75.))
       anyv3d(i,j,k,1)= abs(temobc_t(i,j,k,2)-temobc_t(i,j,k,0))        & !06/03/08
                             *dz_t(i,j,k,1)                           & !10-10-09
                 *(1.-exp( depth_t(i,j,k)/75.)) ! on n'aime pas trop la variabilite de surface...
       anyv3d(i,j,k,2)=-abs(salobc_t(i,j,k,2)-salobc_t(i,j,k,0))        & !06/03/08
                             *dz_t(i,j,k,1)                           & !10-10-09
                 *(1.-exp( depth_t(i,j,k)/75.))

! SUM1 est l'integrale de f(z) entre le fond et la surface c.a.d. F(0)-F(-h)
       sum1=sum1+const1*anyv3d(i,j,k,1)                                 &
                +const2*anyv3d(i,j,k,2)

       enddo

! X1 c'est -dzet*rh0/(F(0)-F(-h))
! fractionnee par le rapport du pas de temps sur l'echelle de temps du transitoire
! fractionne aussi par un effet de profondeur: les zones profondes ayant
! une variabilite rapide et de forte amplitude en sse, on cherche moins a corriger
! ces zones la
      x1=(   ssh_int_w(i,j,1)                                           &
          -sshtotide_w(i,j,1) & ! la maree est otee                     !18/01/08
             -sshobc_w(i,j,1)                                           &
                             -mean_ssh_loc)*rho                               &
        /sign( max(h_w(i,j)*0.0001d0,abs(sum1)) , sum1 )                & !23/02/09
       *const8*sponge_t(i,j,1)                                          &
       *h_w(i,j)
! note: h_z*0.0001 = h_z*abs(const1)*DT avec const1=-0.2 et DT de l'ordre
! de 1E-3 representant un ecart minimum en degres.


       do k=kmin_w(i,j),kmax ! debut de boucle n°2 sur k

       x0=x1*anyv3d(i,j,k,1)            ! correction a priori
       x2=tem_t(i,j,k,1)-temobc_t(i,j,k,1)  ! Ecart a l'ogcm avant correction
!      thz_c(i,j,k,2)=thz_c(i,j,k,2)+x0
       tem_t(i,j,k,2)=tem_t(i,j,k,2)+x0/dz_t(i,j,k,2)                   & !10-10-09
       *min(max(1.-abs(x2)/(abs(3.*anyv3d(i,j,k,1))+small1) & ! fonction d'amplification
                   *sign(un,x0*x2),zero),const3)                        !06/03/09

       x0=x1*anyv3d(i,j,k,2)            ! correction a priori
       x2=sal_t(i,j,k,1)-salobc_t(i,j,k,1)  ! Ecart a l'ogcm avant correction
!      shz_c(i,j,k,2)=shz_c(i,j,k,2)+x0
       sal_t(i,j,k,2)=sal_t(i,j,k,2)+x0/dz_t(i,j,k,2)                   & !10-10-09
       *min(max(1.-abs(x2)/(abs(3.*anyv3d(i,j,k,2))+small1) & ! fonction d'amplification
                   *sign(un,x0*x2),zero),const3)                        !06/03/09

! Role de la fonction d'amplification:
! La correction est diminuee, voir annulee si on est amene a trop s'eloigner
! de l'ogcm et au contraire amplifier si on est a la fois loin de l'ogcm et que
! la correction tend a diminuer l'ecart....


       enddo                 ! fin de boucle n°2 sur k

      endif                          !>>>>>>>>>>>>>>>>
      enddo   ! fin de boucle sur i
      enddo   ! fin de boucle sur j
      enddo   ! fin de boucle sur loop1

      return
      endif                                !eeeeeeeeeeeeeeeeeeeeeeeeeee>

#endif

      stop 'botprescont: your choice is not valid'

      end subroutine botprescont
