      subroutine obc_bio_upd
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
      use module_parallele !#mpi
      implicit none
      real ox_A0,ox_A1,ox_A2,ox_A3,ox_A4,ox_A5,ox_B0,ox_B1,ox_B2,  &
           ox_B3,ox_A,ox_B,ox_C,ox_D,ox_C0,Ts,A
      integer kobc 


! modifs: 22/05/06: correction bug sur K1 (KR au lieu de KR)
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.20 18-04-11  Calculs sur la base d'un temps en secondes


!----------------------------------------------------------------------
!     BUT:
!     ----
!     Ce programme lit les concentrations des nutriments dans un fichier
!     et met a jour à chaque pas de temps
!
!
!     OUTPUT:
!     ------
!     BIO_t
!
!
!
!     MODIFS:
!     -------
!
!-----------------------------------------------------------------------
! mise en service: 15/04/14
!-----------------------------------------------------------------------



! LE CAS PARTICULIER DE L'ETAT INITIAL
! 1. Definir les parametres d'avancee dans le temps (date repere & echantillonage)
! Debut:
      if(iteration3d==kount0) then !.................!       !18-04-11

        do kobc=1,4 ! debut boucle sur les rivieres
         i1=obc_bio_date(1,kobc)
         i2=obc_bio_date(2,kobc)
         i3=obc_bio_date(3,kobc)
         i4=obc_bio_date(4,kobc)
         i5=obc_bio_date(5,kobc)
         i6=obc_bio_date(6,kobc)
         write(6,*)'appel à date_to_kount dans obc_bio_upd'
!        call date_to_kount
         call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!        bioriv_dt(kr,2)=xdtk_out                     ! date repere NC=1 en kount
!        bioriv_dt(kr,1)=bioriv_info(kr)*3600./dti_fw ! duree, en kount, entre 2 infos
         obc_bio_dt(kobc,2)=elapsedtime_out              !18-04-11
         obc_bio_dt(kobc,1)=obc_bio_info(kobc)*3600.        !18-04-11

        enddo          ! fin boucle sur les nutriments

      endif                    !.................!
! Fin.


      do 1000 kobc=1,4 ! boucle sur les nutriments


      ksecu=0 ! A priori on n'ouvre pas le fichier de "donnee riviere"

! Est on de part et d'autre d'une nouvelle echeance?
      x2=(elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)     !18-04-11
      x1=(elapsedtime_bef-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
! Si la reponse est oui alors il faut lire le fichier:
      if(j2-j1.eq.1)         ksecu=1  ! indique que nous lirons le fichier
      if(iteration3d==kount0)ksecu=1  ! L'initialisation impose de lire le fichier !18-04-11

!      if(i==11.and.j==32) then
!      if(kobc==1) then
!      if(mask_t(i,j,kmax)==1) then
!      print*,'test',elapsedtime_now,elapsedtime_bef,x2,x1,j2,j1
!      print*,'test',j2-j1,ksecu,obc_bio_data(kmax,kobc,0),obc_bio_data(kmax,kobc,2)
!      endif
!      endif 
!      endif
      

      if(ksecu.eq.1) then !*************************************>

      if(iteration3d==kount0)then !---------------->        !18-04-11
          k4=0 ! à l'initialisation on lit 2 champs
      else                        !---------------->                       !27/04/04
! sinon l'ancienne echeance est sauvée avant de lire la suivante:
          k4=2
          do k=1,kmax
           obc_bio_data(k,kobc,0)=obc_bio_data(k,kobc,2)
          enddo
      endif                       !---------------->


      print*,'obc_bio_file(kobc)',obc_bio_file(kobc),kobc

! Ouverture du fichier:
      open(unit=3,file=obc_bio_file(kobc)                                &
                 ,access='direct'                                      &
                 ,recl=4*kmax                                          &
                 ,form='unformatted')


! lecture du fichier:
      do k3=k4,2,2
         x2=(elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)      !18-04-11
         nc=int(1.+x2)+k3/2
         read(3,rec=nc)(obc_bio_data(k,kobc,k3),k=1,kmax)
      enddo

      close(3)
      endif               !*************************************>

      if(kobc==1) vb=initrate
      if(kobc==2) vb=iphosphate
      if(kobc==3) vb=isilice
      if(kobc==4) vb=ioxygen
! Interpolation temporelle entre 2 echeances du fichier:
      rap=      (elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)  & !18-04-11
          -int( (elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1) )


      do j=1,jmax
      do i=1,imax

      if(lon_t(i,j)<-5.8*deg2rad) then


      do k=1,kmax 
      if(mask_t(i,j,k)==1) then
       bio_t(i,j,k,vb)=(1.-rap)*obc_bio_data(k,kobc,0)                  &
                      +    rap *obc_bio_data(k,kobc,2)
      endif
      enddo


      if(i==11.and.j==32) then
      if(kobc==1) then
      if(mask_t(i,j,kmax)==1) then
      print*,elapsedtime_now,obc_bio_dt(kobc,2),obc_bio_dt(kobc,1)  
      print*,'obc_bio_upd',i,j,bio_t(i,j,kmax,vb),kobc
      print*,rap,obc_bio_data(kmax,kobc,0),obc_bio_data(kmax,kobc,2)
!     stop
      endif
      endif
      endif

      endif ! lon

      enddo !i
      enddo !j 

 1000 continue ! fin de boucle sur les nutriments

      return
      end
