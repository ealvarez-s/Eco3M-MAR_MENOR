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
      real depth_l(100),obc_bio_data_l(100) !,obc_bio_data_li(100)
      integer kobc,K_pr 


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

        do kobc=1,6 ! debut boucle sur les rivieres
         i1=obc_bio_date(1,kobc)
         i2=obc_bio_date(2,kobc)
         i3=obc_bio_date(3,kobc)
         i4=obc_bio_date(4,kobc)
         i5=obc_bio_date(5,kobc)
         i6=obc_bio_date(6,kobc)
!         write(6,*)'appel à date_to_kount dans obc_bio_upd'
!        call date_to_kount
         call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
!        bioriv_dt(kr,2)=xdtk_out                     ! date repere NC=1 en kount
!        bioriv_dt(kr,1)=bioriv_info(kr)*3600./dti_fw ! duree, en kount, entre 2 infos
         obc_bio_dt(kobc,2)=elapsedtime_out              !18-04-11
         obc_bio_dt(kobc,1)=obc_bio_info(kobc)*3600.        !18-04-11

        enddo          ! fin boucle sur les nutriments

       if(flag_nemoffline.ne.1) then
       open(unit=3,file=obc_bio_filedepth)
       read(3,*)kobcmax
       read(3,*)(depth_l(kobc),kobc=1,kobcmax) 
       close(3)
       do k=1,kobcmax
       depth_biobc(k)=-depth_l(kobcmax-k+1)
!       print*,'depth simed',k,depth_biobc(k)
       enddo
       endif

      endif                    !.................!
! Fin.


!     do 1000 kobc=1,4 ! boucle sur les nutriments
      do 1000 kobc=1,6 ! boucle sur les nutriments + dic et alkalinity

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

!      ig=i+par%timax(1)      ! ds le repere global
!      jg=j+par%tjmax(1)

!      if(ig.eq.61.and.jg.eq.561) then
!      if(par%rank==1) then
!      if(kobc==1) then
!      print*,'testw1',elapsedtime_now,elapsedtime_bef,x2,x1,j2,j1
!      endif 
!      endif
      

      if(ksecu.eq.1) then !*************************************>

      if(iteration3d==kount0)then !---------------->        !18-04-11
          k4=0 ! à l'initialisation on lit 2 champs
      else                        !---------------->                       !27/04/04
! sinon l'ancienne echeance est sauvée avant de lire la suivante:
          k4=2
          if(flag_nemoffline==1) then
            do k=1,kmax
              obc_bio_data(k,kobc,0)=obc_bio_data(k,kobc,2)
            enddo       

          else
            do k=1,kobcmax
             obc_bio_data_li(k,kobc,0)=obc_bio_data_li(k,kobc,2)
            enddo

!      if(par%rank==1) then
!      if(kobc==1) then
!      print*,'testw2',obc_bio_data_li(1,kobc,0),obc_bio_data_li(1,kobc,2),k
!      endif
!      endif

           endif

      endif                       !---------------->

!      if(par%rank==1) then
!      if(kobc==1) then
!      print*,'testw2',obc_bio_data(kmax,kobc,0),obc_bio_data(kmax,kobc,2)
!      endif
!      endif


!      print*,'obc_bio_file(kobc)',obc_bio_file(kobc),kobc

! Ouverture du fichier:
      if(flag_nemoffline==1) then 
      open(unit=3,file=obc_bio_file(kobc)                                &
                 ,access='direct'                                      &
                 ,recl=4*kmax                                          &
                 ,form='unformatted')

      else
!       print*,'passe ici read'
      open(unit=3,file=obc_bio_file(kobc)                             &
                 ,access='direct'                                      &
                 ,recl=4*kobcmax                                          &
                 ,form='unformatted')
      endif

! lecture du fichier:
      do k3=k4,2,2
         x2=(elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)      !18-04-11
         nc=int(1.+x2)+k3/2
         if(flag_nemoffline==1) read(3,rec=nc)(obc_bio_data(k,kobc,k3),k=1,kmax)
         if(flag_nemoffline.ne.1) then 
!        print*,'lire kobcmax= ',kobcmax,nc
         read(3,rec=nc)(obc_bio_data_l(k),k=1,kobcmax)
         do k=1,kobcmax 
         obc_bio_data_li(k,kobc,k3)=obc_bio_data_l(kobcmax-k+1)
!         print*,'obc_bio_data',obc_bio_data_li(k,kobc,k3),kobc
         enddo
         endif
      enddo

      close(3)
      endif               !*************************************>

      if(kobc==1) vb=initrate
      if(kobc==2) vb=iphosphate
      if(kobc==3) vb=isilice
      if(kobc==4) vb=ioxygen
      if(kobc==5) vb=idic
      if(kobc==6) vb=ialkalinity
! Interpolation temporelle entre 2 echeances du fichier:
      rap=      (elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1)  & !18-04-11
          -int( (elapsedtime_now-obc_bio_dt(kobc,2))/obc_bio_dt(kobc,1) )


      do j=1,jmax
      do i=1,imax

      if(mask_t(i,j,kmax)==1) then
!     if(lon_t(i,j)<-5.77*deg2rad) then
      if(lon_t(i,j)<-5.52*deg2rad) then !2019/3/7
!      if(lon_t(i,j)<-6.87*deg2rad) then !2021/12/22

      if(flag_nemoffline==1) then
      do k=1,kmax 
       bio_t(i,j,k,vb)=(1.-rap)*obc_bio_data(k,kobc,0)                  &
                      +    rap *obc_bio_data(k,kobc,2)

       if(vb==iAlkalinity.or.vb==idic.or.vb==ioxygen) bio_t(i,j,k,vb)=bio_t(i,j,k,vb) & ! in micromol/kg
                                                      *(rhp_t(i,j,k)+rho) /1000.
       if(vb==iAlkalinity) bio_t(i,j,k,vb)=bio_t(i,j,k,vb)  & ! total alkalinity
                                  -bio_t(i,j,k,iAmmonium)    &
                                  +bio_t(i,j,k,iPhosphate)   &
                                  +bio_t(i,j,k,iNitrate)

      enddo

!      if(i+par%timax(1)==69.and.j+par%tjmax(1)==565) then
!      if(kobc==1) then
!!      if(mask_t(i,j,kmax)==1) then
!      print*,'interp',bio_t(i,j,kmax,vb),obc_bio_data(kmax,kobc,0),obc_bio_data(kmax,kobc,2),rap
!!     endif
!      endif
!      endif

      else

      
     
       do k=1,kmax 
             if(-depth_t(i,j,k)<=depth_biobc(1))then
               bio_t(i,j,k,vb) =  (1.-rap)*obc_bio_data_li(1,kobc,0) &
                                 +    rap *obc_bio_data_li(1,kobc,2)
!               print*,'cas 1',-depth_t(i,j,k),k,depth_biobc(1)
!               print*,'cas 1',obc_bio_data_li(1,kobc,0),obc_bio_data_li(1,kobc,2)



             elseif(-depth_t(i,j,k)>=depth_biobc(kobcmax))then
               bio_t(i,j,k,vb) = (1.-rap)*obc_bio_data_li(kobcmax,kobc,0) &
                                +    rap *obc_bio_data_li(kobcmax,kobc,2)



!               print*,'cas 2',-depth_t(i,j,k),k,depth_biobc(kobcmax)
!               print*,'cas 2',obc_bio_data_li(kobcmax,kobc,0),obc_bio_data_li(kobcmax,kobc,2)

             else

              K_pr=1
              do while(-depth_t(i,j,k)>depth_biobc(K_pr))
               K_pr=K_pr+1
              enddo

              bio_t(i,j,k,vb)=(                              &
              ( depth_biobc(K_pr  )+depth_t(i,j,k))*((1.-rap)*obc_bio_data_li(K_pr-1,kobc,0) &
                                                    +    rap *obc_bio_data_li(K_pr-1,kobc,2)) &
             +(-depth_biobc(K_pr-1)-depth_t(i,j,k))*((1.-rap)*obc_bio_data_li(K_pr,kobc,0) &
                                                    +    rap*obc_bio_data_li(K_pr,kobc,2)) &
             )/(depth_biobc(K_pr  )-depth_biobc(K_pr-1))


!              print*,'cas 3',-depth_t(i,j,k),depth_biobc(K_pr-1),depth_biobc(K_pr  )
!              print*,'cas 3',obc_bio_data_li(K_pr-1,kobc,0),obc_bio_data_li(K_pr,kobc,0)
!              print*,'cas 3',obc_bio_data_li(K_pr-1,kobc,2),obc_bio_data_li(K_pr,kobc,2)
!              print*,'cas 3',rap
!              print*,'cas 3',bio_t(i,j,k,vb)

            endif !prof3D
            
            if(vb==iAlkalinity.or.vb==idic.or.vb==ioxygen) bio_t(i,j,k,vb)=bio_t(i,j,k,vb) & ! in micromol/kg
                                                           *(rhp_t(i,j,k)+rho) /1000.
            if(vb==iAlkalinity) bio_t(i,j,k,vb)=bio_t(i,j,k,vb)  & ! total alkalinity
                                  -bio_t(i,j,k,iAmmonium)    &
                                  +bio_t(i,j,k,iPhosphate)   &
                                  +bio_t(i,j,k,iNitrate)



!      if(i+par%timax(1)==69.and.j+par%tjmax(1)==565) then
!      if(kobc==1) then
!!      if(mask_t(i,j,kmax)==1) then
!      print*,'interp',bio_t(i,j,k,vb),obc_bio_data_li(k,kobc,0),obc_bio_data_li(k,kobc,2),rap,k
!!     endif
!      endif
!      endif

       enddo 

!      if(i+par%timax(1)==69.and.j+par%tjmax(1)==565) then
!      if(kobc==1) then
!!      if(mask_t(i,j,kmax)==1) then
!      print*,'interp',bio_t(i,j,kmax,vb),obc_bio_data_li(1,kobc,0),obc_bio_data_li(1,kobc,2),rap
!!     endif
!      endif
!      endif


      endif 
      endif ! mask

!     if(i==11.and.j==32) then
!      if(kobc==1) then
!!      if(mask_t(i,j,kmax)==1) then
!      print*,elapsedtime_now,obc_bio_dt(kobc,2),obc_bio_dt(kobc,1)  
!      print*,'obc_bio_upd',i,j,bio_t(i,j,kmax,vb),kobc
!      print*,rap,obc_bio_data(kmax,kobc,0),obc_bio_data(kmax,kobc,2)
!!     stop
!!      endif
!      endif
!      endif

      endif ! lon

      enddo !i
      enddo !j 

 1000 continue ! fin de boucle sur les nutriments


      call mpi_barrier(par%comm2d,k_out)      ! synchro processes 16-09-09

! BUFFER ZONE MER MARMARA
      call scalars_spongelayer_marmara_bio

      return
      end
!---------------------------------------------------------------------------------------------------------

      subroutine scalars_spongelayer_marmara_bio
      use module_principal ; use module_parallele
      use ModuleComBlock
      use UserDeclaration
      use ModuleDeclaration
      implicit none
      real*4 :: c1_,c2_,c3_
      integer :: istr_,jstr_,iend_,jend_

      if(iglb/=1120.or.jglb/=865) &
      stop 'scalars_spongelayer_marmara_bio: iglb/=1120.or.jglb/=865'

      istr_=max(1048-par%timax(1),1)  !avant 2020-09-11: 1056
      iend_=min(1102-par%timax(1),imax) !avant 2020-09-11: 1099
      jstr_=max( 246-par%tjmax(1),1)  !avant 2020-09-11: 245
      jend_=min( 306-par%tjmax(1),jmax) !avant 2020-09-11: 299

      if(elapsedtime_now<864000.)then  ! temps relax court les 10 premiers jours
      c1_= 1.-dti_lp/(86400.*1) !avant 2020-09-11: 3
      else
      c1_= 1.-dti_lp/(86400.*1.)     ! !avant 2020-09-11: ensuite 20 jours 
      endif
      c2_=(1.-c1_)*    timeweightobc(trc_id)
      c3_=(1.-c1_)*(1.-timeweightobc(trc_id))

!     ksecu=0
      do k=1,kmax
      do j=jstr_,jend_
      do i=istr_,iend_

       if (depth_t(i,j,k)>-15.)then
       bio_t(i,j,k,initrate)=bio_t(i,j,k,initrate)*c1_           &
                  +0.24*(1.-c1_)
!      print*,'passe scalar'
       else
       bio_t(i,j,k,initrate)=bio_t(i,j,k,initrate)*c1_           &
                  +1.03*(1.-c1_)           
       endif

       if (depth_t(i,j,k)>-15.)then
       bio_t(i,j,k,iphosphate)=bio_t(i,j,k,iphosphate)*c1_           &
                  +0.06*(1.-c1_)
       else
       bio_t(i,j,k,iphosphate)=bio_t(i,j,k,iphosphate)*c1_           &
                  +0.05*(1.-c1_)
       endif

!     ksecu=1
      enddo
      enddo
      enddo
!     if(ksecu==1)write(10+par%rank,*)'Je sponge marmara'

      end subroutine scalars_spongelayer_marmara_bio

!---------------------------------------------------------------------------------------------------------
