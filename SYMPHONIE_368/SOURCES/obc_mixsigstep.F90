      subroutine obc_mixsigstep(ichoix)
!______________________________________________________________________
! S model
! release S.26 - last update: 31-07-14
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='obc_mixsigstep'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         09/03/09: mis en service
!         20/04/09: Pas possible de coller une C.L. en i=2 etc dans le cas de la
!                   parallelisation. On opte pour une simple C.L. en i=1,...
! 2009.3  05-10-09  ajout d'un "ifdef parallele"
! 2010.10 17-06-10  C.L. de type z1 sur kmin
! 2010.18 23-02-11  ajout C.L. de type z2 sur kmin_w
! 2010.25 28-03-12  echange compatible avec pgf Stelios
! S.26    31-07-14  nouveaux echanges
!...............................................................................

!----------------------------------------------------------------------!10/05/04
! TRAITEMENT SPECIAL DES FRONTIERES POUR QUE LES CONDITIONS DE GRADIENT
! NUL NE SOIENT PAS GENEES PAR UNE MARCHE D'ESCALIER:
! DEBUT:
!----------------------------------------------------------------------!10/05/04
! Même kmin (le + grand) à la frontière
! On n'oublie de modifier la bathy en consequence sinon le dz du bas
! risque d'être énorme et creer des pb pour les obc

       do i=1,imax

        j0=0
        j1=1
        j2=2
        if(kmin_w(i,j1).lt.kmin_w(i,j2)) then !>>>>>>>
           k=kmin_w(i,j2)
             kmin_w(i,j1)=k
                h_w(i,j1)=- depth_w(i,j1,k)
        endif                                 !>>>>>>>
        kmin_w(i,j0)=kmin_w(i,j1)             !17-06-10

        j0=jmax+1
        j1=jmax
        j2=jmax-1
        if(kmin_w(i,j1).lt.kmin_w(i,j2)) then !>>>>>>>
           k=kmin_w(i,j2)
             kmin_w(i,j1)=k
                h_w(i,j1)=- depth_w(i,j1,k)
        endif                                 !>>>>>>>
        kmin_w(i,j0)=kmin_w(i,j1)             !17-06-10

       enddo

       do j=1,jmax

        i0=0
        i1=1
        i2=2
        if(kmin_w(i1,j).lt.kmin_w(i2,j)) then !>>>>>>>
           k=kmin_w(i2,j)
             kmin_w(i1,j)=k
                h_w(i1,j)=- depth_w(i1,j,k)
        endif                                 !>>>>>>>
        kmin_w(i0,j)=kmin_w(i1,j)             !17-06-10

        i0=imax+1
        i1=imax
        i2=imax-1
        if(kmin_w(i1,j).lt.kmin_w(i2,j)) then !>>>>>>>
           k=kmin_w(i2,j)
             kmin_w(i1,j)=k
                h_w(i1,j)=- depth_w(i1,j,k)
        endif                                 !>>>>>>>
        kmin_w(i0,j)=kmin_w(i1,j)             !17-06-10

       enddo

! Coins:
       kmin_w(0     ,0     )=kmin_w(1   ,1   )               !23-02-11
       kmin_w(imax+1,jmax+1)=kmin_w(imax,jmax)               !23-02-11
       kmin_w(imax+1,0     )=kmin_w(imax,1   )               !23-02-11
       kmin_w(0     ,jmax+1)=kmin_w(1   ,jmax)               !23-02-11

! C.L. type z2
      do i=0,imax+1                          !23-02-11
       kmin_w(i,-1    )=kmin_w(i,0     )
       kmin_w(i,jmax+2)=kmin_w(i,jmax+1)
      enddo
      do j=-1,jmax+2
       kmin_w(-1    ,j)=kmin_w(0     ,j)
       kmin_w(imax+2,j)=kmin_w(imax+1,j)
      enddo

!----------------------------------------------------------------------!10/05/04
! TRAITEMENT SPECIAL DES FRONTIERES POUR QUE LES CONDITIONS DE GRADIENT
! NUL NE SOIENT PAS GENEES PAR UNE MARCHE D'ESCALIER:
! FIN.
!----------------------------------------------------------------------!10/05/04

#ifdef parallele
      call get_type_echange('z0','h_w_z0'   ,h_w   ,lbound(h_w)   ,ubound(h_w)   ,i1)
      call get_type_echange('zc','kmin_w_zc',kmin_w,lbound(kmin_w),ubound(kmin_w),i2)
      do loop1=1,subcycle_exchange
        call echange_voisin(h_w   ,i1,mpi_neighbor_list(loop1)) !31-07-14
        call echange_voisin(kmin_w,i2,mpi_neighbor_list(loop1))
      enddo
      call loc_wait()
#endif

      end subroutine obc_mixsigstep
