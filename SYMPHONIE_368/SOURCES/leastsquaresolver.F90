      subroutine leastsquaresolver(ki1_sd,ki2_sd,ki3_sd)
!______________________________________________________________________
! S model
! release S.26 - last update: 13-06-14
!______________________________________________________________________

      use module_principal
      use module_systeme
      implicit none
      integer                                                         &
        ki1_sd                                                        & !04/04/03
       ,ki2_sd                                                        & !04/04/03
       ,ki3_sd                                                        & !04/04/03
       ,i_sd                                                          & !04/04/03
       ,i1_sd                                                         & !04/04/03
       ,ii_sd                                                         & !04/04/03
       ,j_sd                                                          & !04/04/03
       ,k1_sd                                                         & !04/04/03
       ,k2_sd                                                         !04/04/03
#ifdef synopsis
       subroutinetitle='leastsquaresolver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!..............................................................................
! Version date      Description des modifications
!         19/11/01: debug
!         30/01/02: amenagement pour accueuillir analyse harmonique maree pour
!                   laquelle le produit matrice transposee par matrice egalite
!                   n'est pas a faire car deja fait dans ANALYSETIDE_Z
!         08/04/02: suppression instruction inutile
!         04/04/03: les indices de boucles deviennent des variables locales
!                   avec l'extension _SD au bout de leur nom.
!         09/07/04: ajout d'un cas permettant de ne pas refaire le calcul de la
!                   transposée, ainsi qu'une memorisation des pivots qui évite
!                   la longue reorganisation de la matrice principale requise
!                   par la triangularisation.
! 2010.10 16-06-10  surdetermine est renommé leastsquaresolver
! 2010.11 03-07-10  suppression test sur dimsparse
!         18-07-10  La routine est adaptée au fait que l'archivage dans le fichier
!                   time.dat a été remplacé par un calcul en ligne dans
!                   tideanalysismatrix
! S.26    13-06-14  Introduction d'un cas test a decommenter....
!..............................................................................

!     if(nbsparse.gt.dimsparse)stop 'surdetermine.f erreur1'            !15/11/01
      if(nbinco.gt.nbincomax)then
      write(6,*)'nbinco   =',nbinco
      write(6,*)'nbincomax=',nbincomax
      stop 'donc dans surdetermine.f erreur2'                           !19/11/01
      endif



! LA METHODE EST DONNEE PAR NOUGIER
! Pour la validation, un exemple est donne par Nougier page 43
! la solution est 6/7=0.8571 et 3/7=0.4285

!------------------------------------------------------------
! Ces lignes permettent de tester le cas test NOUGIER page 43: !13-06-14
! Debut:
!------------------------------------------------------------
!     nbinco=2
!     nbsparse=6
!     call allocate_systeme(1,        & ! Flag "allocate"
!                           3,        & ! Nombre d'equations (>= nb d'inconnues)
!                           nbsparse, & ! Nombre d'elements non nuls dans la matrice creuse
!                           nbinco  )   ! Nombre d'inconnues
!     ksparse=0
!     kequation=1
!     matrix_b(kequation)=1.
!      kinconnue=1
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=2.
!      kinconnue=2
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=1.
!     kequation=2
!     matrix_b(kequation)=1.
!      kinconnue=1
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=1.
!      kinconnue=2
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=-1.
!     kequation=3
!     matrix_b(kequation)=3.
!      kinconnue=1
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=1.
!      kinconnue=2
!      ksparse=ksparse+1
!          isparse(ksparse)=kequation
!          jsparse(ksparse)=kinconnue
!      sparsevalue(ksparse)=1.
!     call leastsquaresolver(0,0,0)
!     write(6,*)'Solutions=',tab7(1:2)
!     call allocate_systeme(2,0,0,0) ! Flage "deallocate"
!------------------------------------------------------------
! Ces lignes permettent de tester le cas test NOUGIER page 43:
! Fin.
!------------------------------------------------------------


!     WRITE(6,*)'avant produit transposée'
!............................................................................
! PRODUIT DE LA MATRICE DE DEPART PAR SA TRANSPOSEE
! AFIN D'OBTENIR LA MATRICE CARREE DU SYSTEME:

! cas standart:
      if(ki1_sd.eq.0) then !***************************************>   !30/01/02

      do i_sd=1,nbinco
        tab4(i_sd)=0.
          do j_sd=1,nbinco
            tab3(i_sd,j_sd)=0.
          enddo
      enddo

      do 300 k1_sd=1,nbsparse

       tab4(jsparse(k1_sd))=tab4(jsparse(k1_sd))                        &
       +sparsevalue(k1_sd)* matrix_b(isparse(k1_sd))

        do 200 k2_sd=1,nbsparse

         if(isparse(k2_sd).eq.isparse(k1_sd))                           &
         tab3(jsparse(k1_sd),jsparse(k2_sd))=                           &
         tab3(jsparse(k1_sd),jsparse(k2_sd))                            &
         +sparsevalue(k1_sd)*sparsevalue(k2_sd)

  200   continue
  300   continue

      do j_sd=1,nbinco
      jpos1(j_sd)=j_sd
      jpos2(j_sd)=j_sd
      enddo
      do i_sd=1,nbinco
      do j_sd=1,nbinco
      tab5(i_sd,j_sd)=tab3(i_sd,j_sd)
      enddo
      tab6(i_sd)=tab4(i_sd)
      enddo


      endif             !******************************************>   !30/01/02

! cas special harmonique marees:
! 1: TAB4 déjà traité en cours de simulation par symphonie
! 2: La matrice principale est la même en tout point, ce qui signifie
!    que certains calculs, couteux, peuvent nêtre fait qu'une fois.
!    C'est le cas du calcul du produit de la matrice principale
!    par sa transposée. Ce calcul produit TAB3 qui est copié en
!    double dans TAB5. Au 2eme passage dans cette routine KI1_SD=2
!    et on recharge TAB3 avec
      if(ki1_sd.eq.1) then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>   !30/01/02

!     do i_sd=1,nbinco
!         do j_sd=1,nbinco
!           tab3(i_sd,j_sd)=0.
!         enddo
!     enddo

!     do 301 k1_sd=1,nbsparse

!TAB4(JSPARSE(K1_SD)) est calculé en amont dans la routine appelant surdetermine
!       do 201 k2_sd=1,nbsparse

!        isparse(k2_sd)=isparse(k2_sd)
!        jsparse(k2_sd)=jsparse(k2_sd)

!        if(isparse(k2_sd).eq.isparse(k1_sd))                           &
!        tab3(jsparse(k1_sd),jsparse(k2_sd))=                           &
!        tab3(jsparse(k1_sd),jsparse(k2_sd))                            &
!       +sparsevalue(k1_sd)*sparsevalue(k2_sd)

! 201   continue
! 301   continue

      do j_sd=1,nbinco
      jpos1(j_sd)=j_sd
      jpos2(j_sd)=j_sd
      enddo
      do i_sd=1,nbinco
      do j_sd=1,nbinco
      tab5(i_sd,j_sd)=tab3(i_sd,j_sd)
      enddo
      tab6(i_sd)=tab4(i_sd)
      enddo

      endif             !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>   !30/01/02

!............................................................................


      if(ki1_sd.le.1) then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>   !09/07/04

      do 20 i_sd=1,nbinco
! recherche du plus grand pivot non nul
!     JNELLA=0
      valeur1=abs(tab3(i_sd,i_sd))

       do 22 j_sd=i_sd,nbinco
       valeur2=abs(tab3(i_sd,j_sd))
       if(valeur2.ge.valeur1) then
       jpivot(i_sd)=j_sd
       valeur1=valeur2
       endif
   22 continue


      do 21 i1_sd=1,nbinco
      col1(i1_sd)=tab3(i1_sd,i_sd)
      col2(i1_sd)=tab3(i1_sd,jpivot(i_sd))
      tab3(i1_sd,i_sd)=col2(i1_sd)
      tab3(i1_sd,jpivot(i_sd))=col1(i1_sd)
  21  continue
      jpos2(i_sd)=jpos1(jpivot(i_sd))
      jpos2(jpivot(i_sd))=jpos1(i_sd)
      jpos1(i_sd)=jpos2(i_sd)
      jpos1(jpivot(i_sd))=jpos2(jpivot(i_sd))

!-------------------------------
      pivot(i_sd)=tab3(i_sd,i_sd)
      do 23 j_sd=1,nbinco
      tab3(i_sd,j_sd)=tab3(i_sd,j_sd)/pivot(i_sd)
  23  continue
      tab4(i_sd)=tab4(i_sd)/pivot(i_sd)
!-------------------------------


!-------------------------------
      do 24 ii_sd=i_sd+1,nbinco
      pilier(ii_sd,i_sd)=tab3(ii_sd,i_sd)
      tab4(ii_sd)=tab4(ii_sd)-tab4(i_sd)*pilier(ii_sd,i_sd)
      do 24 j_sd=1,nbinco
       tab3(ii_sd,j_sd)=                                                &
       tab3(ii_sd,j_sd)-tab3(i_sd,j_sd)*pilier(ii_sd,i_sd)
  24  continue
!-------------------------------


!     WRITE(6,*)'PIVOT(I_SD)',PIVOT(I_SD)
!     DO I2_SD=1,NBINCO
!     WRITE(6,100)(TAB3(I2_SD,J2),J2=1,NBINCO)
!     ENDDO
!     WRITE(6,*)

  20  continue

      else                 !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>   !09/07/04

      do 40 i_sd=1,nbinco

!-------------------------------
      tab4(i_sd)=tab4(i_sd)/pivot(i_sd)
!-------------------------------


!-------------------------------
      do 44 ii_sd=i_sd+1,nbinco
      tab4(ii_sd)=tab4(ii_sd)-tab4(i_sd)*pilier(ii_sd,i_sd)
  44  continue
!-------------------------------


  40  continue

      endif                !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>   !09/07/04

!     WRITE(6,*)'ordre des inconnues'
!     DO 26 I_SD=1,NBINCO
!     WRITE(6,*)'en ',I_SD,' inconnue ',JPOS1(I_SD)
! 26  CONTINUE

  100 format(3f8.3)
  102 format(2f8.3)
  101 format(f8.3)

!________________________________________________________!
!  resolution du systeme:
      tabx(nbinco)=tab4(nbinco)/tab3(nbinco,nbinco)
      i_sd=nbinco
!cccccWRITE(6,*)'inconnue ',JPOS1(I_SD),' =',TABX(I_SD)
      do 27 i_sd=nbinco-1,1,-1
       x1=0.
       do 28 j_sd=i_sd+1,nbinco
        x1=x1+tabx(j_sd)*tab3(i_sd,j_sd)
   28 continue
       tabx(i_sd)=tab4(i_sd)-x1
!cccccWRITE(6,*)'inconnue ',JPOS1(I_SD),' =',TABX(I_SD)
   27 continue

      do j_sd=1,nbinco
      tab7(jpos1(j_sd))=tabx(j_sd)
      enddo
!________________________________________________________!

      end subroutine leastsquaresolver
