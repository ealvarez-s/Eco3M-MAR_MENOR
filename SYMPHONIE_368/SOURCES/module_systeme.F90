      module module_systeme
      use module_principal
      implicit none
!______________________________________________________________________
! S model
! release S26 - last update: 11-11-14
!______________________________________________________________________
! modifs: 24/07/03: dimension de NAMETIDE aggrandie pour recevoir
!                   nom de fichier supplementaire (fichier de courant)
!        21/08/03: ajout d'une possibilité de selectionner les forces de
!                  marée: ondes + tout potentiel
!                          ondes + potentiel astro
!                          ondes + potentiel charge
!                          ondes + aucun potentiel
!                  c.a.d. bienvenue à TIDEFORCES
!                         bienvenue à TIDEVEL (choix courant fichier ou lin)
!        28/11/07: ajout de nouveax tableaux
! S.26   24-02-13 a cause de l'absence de implicit et de l'appel
!                 a module principal imax et jmax n'etaient pas correctement
!                 definis. allocate_systeme passe dans le module_systeme
!                 Initialiser les tableaux apres allocation
!______________________________________________________________________
!       DEBUT COMMONS - NE SURTOUT PAS TOUCHER CETTE LIGNE SVP!

!______________________________________________________________________
! LES NOMBRES EN DOUBLE PRECISION:

      double precision,dimension(:,:),allocatable ::                    &
        pilier                                                          &
       ,tab3                                                            &
       ,tab5

      double precision,dimension(:),allocatable ::                      &
        matrix_b                                                        &
       ,sparsevalue                                                     &
       ,tab4                                                            &
       ,tabx                                                            &
       ,tab6                                                            &
       ,tab7                                                            &
       ,col1                                                            &
       ,col2                                                            &
       ,pivot

      double precision ::                                               &
        valeur1                                                         &
       ,valeur2

!______________________________________________________________________
! LES NOMBRES ENTIERS:

      integer,dimension(:,:,:),allocatable ::                           &
        numsyst

      integer,dimension(:),allocatable ::                               &
        isparse                                                         &
       ,jsparse                                                         &
       ,jpos1                                                           &
       ,jpos2                                                           &
       ,jpivot                                                          &
       ,isyst                                                           &
       ,jsyst

      integer ::                                                        &
        kinconnue                                                       &
       ,kequation                                                       &
       ,ksparse

!       FIN COMMONS   - NE SURTOUT PAS TOUCHER CETTE LIGNE SVP!
!______________________________________________________________________

contains

      subroutine allocate_systeme(ichoix,dim1_,dim2_,dim3_)
      implicit none
      integer ichoix,dim1_,dim2_,dim3_
#ifdef synopsis
       subroutinetitle='allocate_systeme'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ichoix.eq.1) then !111111111111111>

      allocate(matrix_b      (dim1_))           ; matrix_b=0.
      allocate(isyst         (dim1_/2))         ; isyst=0
      allocate(jsyst         (dim1_/2))         ; jsyst=0
      allocate(numsyst       (imax+1,jmax+1,1)) ; numsyst=0
      allocate(sparsevalue   (dim2_))           ; sparsevalue=0.
      allocate(isparse       (dim2_))           ; isparse=0
      allocate(jsparse       (dim2_))           ; jsparse=0

      allocate(pilier(dim3_,dim3_)) ; pilier=0.
      allocate(tab3  (dim3_,dim3_)) ; tab3=0.
      allocate(tab5  (dim3_,dim3_)) ; tab5=0.

      allocate(tab4  (dim3_))       ; tab4=0 
      allocate(tabx  (dim3_))       ; tabx=0
      allocate(tab6  (dim3_))       ; tab6=0
      allocate(tab7  (dim3_))       ; tab7=0
      allocate(col1  (dim3_))       ; col1=0
      allocate(col2  (dim3_))       ; col2=0
      allocate(pivot (dim3_))       ; pivot=0

      allocate(jpos1 (dim3_))       ; jpos1=0
      allocate(jpos2 (dim3_))       ; jpos2=0
      allocate(jpivot(dim3_))       ; jpivot=0

      endif                !111111111111111>

      if(ichoix.eq.2) then !222222222222222>

      deallocate(matrix_b      )
      deallocate(isyst         )                                       !28/11/08
      deallocate(jsyst         )                                       !28/11/08
      deallocate(numsyst       )                                       !28/11/08
      deallocate(sparsevalue   )
      deallocate(isparse       )
      deallocate(jsparse       )

      deallocate(pilier)
      deallocate(tab3  )
      deallocate(tab5  )

      deallocate(tab4  )
      deallocate(tabx  )
      deallocate(tab6  )
      deallocate(tab7  )
      deallocate(col1  )
      deallocate(col2  )
      deallocate(pivot )

      deallocate(jpos1 )
      deallocate(jpos2 )
      deallocate(jpivot)

      endif                !222222222222222>

      end subroutine allocate_systeme

      end module module_systeme
