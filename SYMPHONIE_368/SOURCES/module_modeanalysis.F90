      module module_modeanalysis
!______________________________________________________________________
! S model
! release S26.1 - last update: 14-04-16
!______________________________________________________________________
! S.26  19-04-13  Accorder un temps de spin-up avant de commencer l'analyse
!       14-01-15  - adaptation aux nouveaux tableau de tridiagonal solver...
!                 - archiver quelques profils dans des fichiers classes 
!       10-02-16  plus d'infos dans les fichiers comme les vitesses de phases
!                 A la frequence de M2 
!       12-04-16  timewindow_3dwaves fenetre temporelle d'analyse
!       14-04-16  ajout possibilite d'inclure la pression 3D dans l'analyse 
!                 harmonique
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none
      real cwave_guess
      integer :: countmodemax=3          &
                ,flag_3dwaves_modes      &
                ,flag_3dwaves_harmonics  &
                ,flag_3dwaves_pharmonics &
                ,analysiscount_3dwaves=0
      double precision :: timewindow_3dwaves=2.*86400. ! start analyse xxx seconds before simulation end 

contains

      subroutine modeanalysis_roadmap
      implicit none
#ifdef synopsis
       subroutinetitle='modeanalysis_roadmap'
       subroutinedescription= &
       'Sturm-Liouville 3D Waves Modes Analysis Driver'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! En attendant la possible creation d'un pilotage depuis un notebook, on definit
! les parametres d'entree en dur dans la presente routine (voir lignes en suivant)

! Si flag=1 calculer les modes propres baroclines
      flag_3dwaves_modes=0
! Si flag=1 calculer l'amplitude et la phase des modes pour les periodes de la maree
      flag_3dwaves_harmonics=0  ! Pour les courants
      flag_3dwaves_pharmonics=0 ! Pour la pression (de densitE) !14-04-16

! Ne pas faire l'analyse 3D si kmax=1: !13-04-16
      if(kmax==1) then  !pmx>
           flag_3dwaves_modes=0 
           flag_3dwaves_harmonics=0 
           flag_3dwaves_pharmonics=0
      endif             !pmx>

      if(flag_3dwaves_modes==1) then  !--------->

       allocate(c_wave_mode(imax,jmax     ,countmodemax)) ; c_wave_mode=0
       allocate(ucoefmode_t(imax,jmax     ,countmodemax)) ; ucoefmode_t=0
       allocate(vcoefmode_t(imax,jmax     ,countmodemax)) ; vcoefmode_t=0
       allocate( uv_wmode_t(imax,jmax,kmax,countmodemax)) ; uv_wmode_t=0

       if(flag_3dwaves_pharmonics==1) then
         allocate(pcoefmode_t(imax,jmax,countmodemax)) ; pcoefmode_t=0. !14-04-16
       endif

       call modeanalysis_verticalmodes

      endif                           !--------->

      if(flag_3dwaves_harmonics==1) then  !////////>

        allocate(  analysis_u3d_t(imax,jmax,countmodemax,kmaxtide*2+1));  analysis_u3d_t=0.
        allocate(  analysis_v3d_t(imax,jmax,countmodemax,kmaxtide*2+1));  analysis_v3d_t=0.
        allocate(   u3dmode_cos_t(imax,jmax,countmodemax,kmaxtide+1));   u3dmode_cos_t=0.
        allocate(   u3dmode_sin_t(imax,jmax,countmodemax,kmaxtide  ));   u3dmode_sin_t=0.
        allocate(   v3dmode_cos_t(imax,jmax,countmodemax,kmaxtide+1));   v3dmode_cos_t=0.
        allocate(   v3dmode_sin_t(imax,jmax,countmodemax,kmaxtide  ));   v3dmode_sin_t=0.
        allocate(analysis3dmatrix(kmaxtide*2+1,kmaxtide*2+1        ));analysis3dmatrix=0.

      endif                               !////////>

      if(flag_3dwaves_pharmonics==1) then !////////>
        allocate(  analysis_p3d_t(imax,jmax,countmodemax,kmaxtide*2+1));  analysis_p3d_t=0.
        allocate(   p3dmode_cos_t(imax,jmax,countmodemax,kmaxtide+1));   p3dmode_cos_t=0.
        allocate(   p3dmode_sin_t(imax,jmax,countmodemax,kmaxtide  ));   p3dmode_sin_t=0.
      endif                               !////////>

      end subroutine modeanalysis_roadmap

!....................................................................

      subroutine modeanalysis_verticalmodes
      implicit none
      integer :: i1_,j1_,countmode_,flag_screenoutput_
#ifdef synopsis
       subroutinetitle='modeanalysis_verticalmodes'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Equation de depart:

! Methode 1: 
! Equation 6.11.18 de Gill:
! d/dz(dX/dz)-grav/(rho*(c**2))*drhp/dz*X=0
! que l'on discretise pour un obtenir un systeme d'equation lineaire
! dont on cherche les valeurs de c qui annule le determinant du systeme
! Les C.L. au fond et à la surface sont X=0

! Methode 2:
! Equation 6.11.17 de Gill 1982:
! (1./rhp)*d/dz[rhp*(dX/dz)]-grav/(rhp*(c**2))*drhp/dz*X=0
! c.a.d.    d/dz[rhp*(dX/dz)]-grav/(c**2)*drhp/dz*X=0 

! Rappel de l'alternance des variables le long de l'axe vertical (pointant vers le haut):

! -------------   z_w(k+1)
!  - - - - - -    rhp_t(k) , z_t(k)
! -------------   z_w(k  )
!  - - - - - -    rhp_t(k-1) , z_t(k-1)
! -------------   z_w(k-1)
       


      i1_=1  ; j1_=1

      do 100 j1=1,jmax ! j1_,j1_
      do 100 i1=1,imax ! i1_,i1_
!     do 100 j1=j1_,j1_
!     do 100 i1=i1_,i1_

! Afficher des resultats si flag_screenoutput_=1
      flag_screenoutput_=0
      if(par%rank==0.and.mod(i1,40)==0         &
                    .and.mod(j1,40)==0         &
                    .and.mask_t(i1,j1,kmax)==1)flag_screenoutput_=1

      sum2=0.
      countmode_=0
      cwave_guess=4. ! reset 3m/s

  202 continue

!     if(par%rank==0)write(6,*)cwave_guess

!     const2=-grav/(rho*cwave_guess**2) ! Methode 1
      const2=-grav/(    cwave_guess**2) ! Methode 2

      do 200 k=2,kmax

!     x1=1./( depth_w(i1,j1,k  )- depth_w(i1,j1,k-1))  ! Methode 1
!     x3=1./( depth_w(i1,j1,k+1)- depth_w(i1,j1,k  ))  ! Methode 1
      x1=(rho+rhp_t(i1,j1,k-1))/( depth_w(i1,j1,k  )- depth_w(i1,j1,k-1))  ! Methode 2
      x3=(rho+rhp_t(i1,j1,k  ))/( depth_w(i1,j1,k+1)- depth_w(i1,j1,k  ))  ! Methode 2
      x2=1./( depth_t(i1,j1,k  )- depth_t(i1,j1,k-1))

! La fonction "min" blinde le cas mélangé:
!     x0=const2*min((rhp_t(i1,j1,k)-   rhp_t(i1,j1,k-1)),-small1)      &
!               /( depth_t(i1,j1,k)- depth_t(i1,j1,k-1))
! Attention il est important que rhp_t soit bien la densite potentielle (et non
! pas la densite vraie)
      x0=const2*min(    (rhp_t(i1,j1,k)-   rhp_t(i1,j1,k-1))          &
                    /( depth_t(i1,j1,k)- depth_t(i1,j1,k-1)) ,-small1)

      tridia_in(i1,j1,k,1)=x1*x2*  mask_t(i1,j1,k-1)
      tridia_in(i1,j1,k,3)=x3*x2*  mask_t(i1,j1,k-1)
      tridia_in(i1,j1,k,2)=x0-tridia_in(i1,j1,k,1)-tridia_in(i1,j1,k,3)

      x4=sqrt((tridia_in(i1,j1,k,1)**2                                 &
              +tridia_in(i1,j1,k,2)**2                                 &
              +tridia_in(i1,j1,k,3)**2)/3.)

      tridia_in(i1,j1,k,1)=tridia_in(i1,j1,k,1)/x4
      tridia_in(i1,j1,k,2)=tridia_in(i1,j1,k,2)/x4
      tridia_in(i1,j1,k,3)=tridia_in(i1,j1,k,3)/x4

       tridia_in(i1,j1,k,4)=0./x4

  200 continue

       tridia_in(i1,j1,kmax+1,1)=0.
       tridia_in(i1,j1,kmax+1,2)=1.
       tridia_in(i1,j1,kmax+1,3)=0.
       tridia_in(i1,j1,kmax+1,4)=0.

       tridia_in(i1,j1,kmin_w(i1,j1),1)=0.
       tridia_in(i1,j1,kmin_w(i1,j1),2)=1.
       tridia_in(i1,j1,kmin_w(i1,j1),3)=0.
       tridia_in(i1,j1,kmin_w(i1,j1),4)=0.

      call tridiagonalsolver(1,i1,i1,j1,j1,kmax+1)                     !30/03/06

!     sum1=anyv1d(kmax+1,1) ! anyv1d(1) = 1/tridia_in(0)
      sum1=tridia_in(i1,j1,kmax+1,0)        !14-01-15
      do 201 k=kmax,kmin_w(i1,j1),-1
!     sum1=sum1*anyv1d(k,1)
      sum1=sum1*tridia_in(i1,j1,k,0)
  201 enddo

        if(sum1*sum2<0.)then !----------------------------->
          countmode_=countmode_+1
          c_wave_mode(i1,j1,countmode_)=cwave_guess

          if(flag_screenoutput_==1) then !>>>
           if(countmode_==1)write(6,*)'indices=',i1,j1
           write(6,*)'celerite du mode ',countmode_ &
                    ,' =',c_wave_mode(i1,j1,countmode_)
          endif                          !>>>

          anyv3d(i1,j1,kmin_w(i1,j1)  ,countmode_)=0.
          anyv3d(i1,j1,kmin_w(i1,j1)+1,countmode_)=0.1

           do k=kmin_w(i1,j1)+1,kmax

            anyv3d(i1,j1,k+1,countmode_)=                         &
        (tridia_in(i1,j1,k,4)                                     &
        -tridia_in(i1,j1,k,2)*anyv3d(i1,j1,k  ,countmode_)        &
        -tridia_in(i1,j1,k,1)*anyv3d(i1,j1,k-1,countmode_))       &
        /tridia_in(i1,j1,k,3)

          enddo

! Sous le fond w est zéro:
          do k=1,kmin_w(i1,j1)-1
           anyv3d(i1,j1,k,countmode_)=0.
          enddo
! Ici on chipote: ces quelques lignes servent à s'assurer que w(kmax+1)=0 exactement.
          do k=kmin_w(i1,j1),kmax+1
!         rap=real(k-kmin_w(i1,j1))/real(kmax+1-kmin_w(i1,j1))
          rap=(depth_w(i1,j1,k     )-depth_w(i1,j1,kmin_w(i1,j1))) &
             /(depth_w(i1,j1,kmax+1)-depth_w(i1,j1,kmin_w(i1,j1)))   !14-01-15
          anyv3d(i1,j1,k,countmode_)=anyv3d(i1,j1,k,countmode_)-rap*anyv3d(i1,j1,kmax+1,countmode_)
          enddo

        endif               !-------------------------------->

      cwave_guess=cwave_guess-0.01
      sum2=sum1
      if(cwave_guess>0.01  &   ! En dessous de 0.01m/s on arrête
         .and.countmode_<countmodemax)goto 202

! Si couche mélangée l'identification des modes est inaboutie (et on impose des valeurs bidon):
      if(countmode_<countmodemax) then !fffffffff>
          c_wave_mode(i1,j1,:)=0.
          do countmode_=1,countmodemax
          do k=1,kmax+1
           anyv3d(i1,j1,k,countmode_)=sin(pi*real(k-1)/real(kmax)*real(countmode_))
          enddo
          enddo
      endif                             !ffffffff>

  100 continue

      do j1=1,jmax
      do i1=1,imax

      flag_screenoutput_=0
      if(i1==imax/2.and.j1==jmax/2.and.mask_t(i1,j1,kmax)==1)flag_screenoutput_=1

      if(flag_screenoutput_==1) then !>>>>>
       k2=index(directory,'/',back=.true.)
       write(texte60,'(a,a,i0)')trim(directory(1:k2)) &
                               ,'3DMODES/phasespeed_modes_dom',par%rank
       open(unit=3,file=texte60)
        write(3,'(a,i0,a,i0,1x,i0)')'Phase speed in m*s-1 of the ',countmodemax &
       ,' first baroclinic modes in the nonrotating frame at i,j=',i1+par%timax(1),j1+par%tjmax(1)
         write(3,'(6(e13.6,1x))') c_wave_mode(i1,j1,1:countmodemax)

        write(3,*)'DT & LF DTmax at this point',dti_fw,0.5*min(dx_t(i1,j1),dy_t(i1,j1))/max(c_wave_mode(i1,j1,1),small1)

       write(3,'(a)')'----------------------------------------'
       write(3,'(a,1x,f7.3)')'Coriolis period in hours:',2.*pi/coriolis_t(i1,j1)/3600. !10-02-16
       write(3,'(a)')        'M2 period in hours      :  12.42'
       write(3,'(a,1x,e13.6)')'1/sqrt(1-(f/M2freq)**2) factor='   &
       ,1./sqrt(1.-(12.42*coriolis_t(i1,j1)*3600./2./pi)**2)
       write(3,'(a,a)')'Corresponding phase speed in the rotating' &
                      ,' frame at M2freq:'
       write(3,'(6(e13.6,1x))') c_wave_mode(i1,j1,1:countmodemax) &
                 /sqrt(1.-(12.42*coriolis_t(i1,j1)*3600./2./pi)**2)


       close(3)

       write(texte60,'(a,a,i0)')trim(directory(1:k2)) &
                               ,'3DMODES/3d_w_tmp_modes_dom',par%rank
       open(unit=3,file=texte60)
       write(3,'(a,i0,a,i0,1x,i0)')'depth and ',countmodemax &
      ,' first w modes before normalisation i,j=',i1+par%timax(1),i1+par%tjmax(1)
        do k=1,kmax+1
         write(3,'(6(e13.6,1x))')depth_w(i1,j1,k),anyv3d(i1,j1,k,1:countmodemax)
        enddo
       close(3)
      endif                          !>>>>>

! ETAPE 1: premiere estimation des modes:
      do loop1=1,countmodemax
       do k=1,kmax
         uv_wmode_t(i1,j1,k,loop1)=                               &
            (anyv3d(i1,j1,k+1,loop1)-anyv3d(i1,j1,k,loop1))       &
         /( depth_w(i1,j1,k+1)    - depth_w(i1,j1,k))             &
          *  mask_t(i1,j1,k)
       enddo
      enddo

      if(flag_screenoutput_==1) then !>>>>>
       k2=index(directory,'/',back=.true.)
       write(texte60,'(a,a,i0)')trim(directory(1:k2)) &
                               ,'3DMODES/3d_uv_tmp_modes_dom',par%rank
       open(unit=3,file=texte60)
       write(3,'(a,i0,a,i0,1x,i0)')'depth and ',countmodemax &
        ,' first horizontal velocity modes before normalisation i,j=' &
        ,i1+par%timax(1),j1+par%tjmax(1)
        do k=1,kmax
        write(3,'(6(e13.6,1x))')depth_t(i1,j1,k),uv_wmode_t(i1,j1,k,1:countmodemax)
        enddo
       close(3)
      endif                          !>>>>>

! ETAPE 2: assurer l'orthogonalité des modes

! assurer l'orthogonalité des modes baroclines entre eux
      do k1=2,countmodemax
      do k2=1,k1-1
!     if(par%rank==0)write(6,*)'k1 & k2=',k1,k2

      sum1=0.
      sum2=0.
      do k=kmin_w(i1,j1),kmax
       sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k2)     &
                                *uv_wmode_t(i1,j1,k,k1)
       sum2=sum2+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k2)**2
      enddo
      sum2=max(sum2,small1) !15-01-15

!     if(par%rank==0)write(6,*)sum1,sum2
      do k=kmin_w(i1,j1),kmax
       uv_wmode_t(i1,j1,k,k1)=uv_wmode_t(i1,j1,k,k1)         &
                             -uv_wmode_t(i1,j1,k,k2)*sum1/sum2
      enddo

      enddo
      enddo


! assurer l'orthogonalité des modes baroclines avec le mode barotrope
      do k1=1,countmodemax
       sum1=0.
       sum2=0.
       do k=kmin_w(i1,j1),kmax
        sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k1)
        sum2=sum2+dz_t(i1,j1,k,1)
       enddo
       do k=kmin_w(i1,j1),kmax
        uv_wmode_t(i1,j1,k,k1)=uv_wmode_t(i1,j1,k,k1)-sum1/sum2
       enddo
      enddo

! ETAPE 3: normalisation des modes
      do loop1=1,countmodemax
      sum1=0.
      do k=kmin_w(i1,j1),kmax
       sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,loop1)**2
      enddo
      x1=sqrt(1./max(sum1,small2))
      sum1=0.
      do k=kmin_w(i1,j1),kmax
      uv_wmode_t(i1,j1,k,loop1)=x1*uv_wmode_t(i1,j1,k,loop1)
          sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,loop1)**2
      enddo
!     if(par%rank==0)write(6,*)'verif norme',loop1,sum1
      enddo

      if(flag_screenoutput_==1) then !>>>>>
       k2=index(directory,'/',back=.true.)
       write(texte60,'(a,a,i0)')trim(directory(1:k2)) &
                               ,'3DMODES/3d_uv_modes_dom',par%rank
       open(unit=3,file=texte60)
       write(3,'(a,i0,a,i0,1x,i0)')'depth and ',countmodemax &
        ,' first horizontal velocity modes after post-treatment i,j=' &
        ,i1+par%timax(1),j1+par%tjmax(1)
        do k=1,kmax
         write(3,'(6(e13.6,1x))')depth_t(i1,j1,k),uv_wmode_t(i1,j1,k,1:countmodemax)
        enddo
       close(3)
      endif                          !>>>>>
!     if(i1==i1_.and.j1==j1_) then !>>>>>
!      do k=1,kmax
!      write(100+par%rank,'(6(e13.6,1x))') depth_t(i1,j1,k)          &
!        ,uv_wmode_t(i1,j1,k,1)                                     &
!        ,uv_wmode_t(i1,j1,k,2)                                     &
!        ,uv_wmode_t(i1,j1,k,3)                                     &
!        ,uv_wmode_t(i1,j1,k,4)                                     &
!        ,uv_wmode_t(i1,j1,k,5)
!      enddo
!     endif                     !>>>>>

! Verif baroclines entre eux:
      do k1=2,countmodemax
      do k2=1,k1-1
      sum1=0.
      sum2=0.
      do k=1,kmax
      sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k2)    &
                               *uv_wmode_t(i1,j1,k,k1)
      sum2=sum2+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k2)**2
      enddo
!     if(i1==imax/2.and.j1==jmax/2) &
!     write(10+par%rank,*)'k1 k2 sum1 sum2 ',k1,k2,sum1,sum2
      enddo
      enddo
! verif baroclines avec barotrope
      do k1=1,countmodemax
      sum1=0.
      do k=1,kmax
      sum1=sum1+dz_t(i1,j1,k,1)*uv_wmode_t(i1,j1,k,k1)
      enddo
!     if(par%rank==0)write(6,*)'k1 sum1 ',k1,sum1
      enddo

      enddo
      enddo ! fin de boucles sur i et j

!     call graph_out
!     stop 'coco dingo'

      end subroutine modeanalysis_verticalmodes

!......................................................................

      subroutine modeanalysis_modeprojection
      implicit none
      integer loopm_
#ifdef synopsis
       subroutinetitle='modeanalysis_modeprojection'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_3dwaves_modes==0)return

      do loopm_=1,countmodemax

       do j=1,jmax
       do i=1,imax
        ucoefmode_t(i,j,loopm_)=0.
        vcoefmode_t(i,j,loopm_)=0.
        do k=1,kmax
         ucoefmode_t(i,j,loopm_)=ucoefmode_t(i,j,loopm_)+dz_t(i,j,k,1)*0.5*(vel_u(i,j,k,1)+vel_u(i+1,j,k,1))*uv_wmode_t(i,j,k,loopm_)
         vcoefmode_t(i,j,loopm_)=vcoefmode_t(i,j,loopm_)+dz_t(i,j,k,1)*0.5*(vel_v(i,j,k,1)+vel_v(i,j+1,k,1))*uv_wmode_t(i,j,k,loopm_)
        enddo
       enddo
       enddo

      enddo


      end subroutine modeanalysis_modeprojection

!......................................................................

      subroutine modeanalysis_harmonics(case_)
      use module_systeme
      implicit none
      integer loop,loopm_,case_
#ifdef synopsis
       subroutinetitle='modeanalysis_harmonics'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      if(kmaxtide==0)return                                        !28/07/04
!     if(elapsedtime_now<10.*86400.)return ! temps de spin-up avant de commencer l'analyse !19-04-13
      if(elapsedtime_now<elapsedtime_end-timewindow_3dwaves)return !12-04-16

      if(case_==1) then !11111111111111111111111111>

      analysiscount_3dwaves=analysiscount_3dwaves+1

! Projeter le courant sur les modes baroclines verticaux:
      call modeanalysis_modeprojection

! Analyse en harmoniques du signal
! pour eviter un stokage volumineux et gagner de la cpu le produit de la
! transposée de la matrice principale par la matrice egalite est réalisé
! progressivement au cours de la simulation.

      kinconnue=0
      do 426 ktide=1,kmaxtide+1

       time2=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))    &  !15-04-16
             +v0tide(ktide)+utide(ktide,1)

      if(ktide==kmaxtide+1) then !ooo>
       k2=1 ! background frequence nulle
      else                       !ooo>
       k2=2 ! cas ondes de marees
      endif                      !ooo>

      do 426 k1=1,k2 ! k2=2 sauf pour background oU k2=1
       kinconnue=kinconnue+1

       x1=(     mod(k1,2) *cos(time2)                    &
           +(1.-mod(k1,2))*sin(time2) )*ftide(ktide,1)

       do 425 loopm_=1,countmodemax
       do 425 j=1,jmax
       do 425 i=1,imax

! comme la routine est appelee avant external_mode on peut considerer
! que velbar(1) et ssh_int(2) sont la solution au t "now" d'où
! le calcul de time2 basé sur kount (et non pas kount+1 ou autre...)
        analysis_u3d_t(i,j,loopm_,kinconnue)=     &
        analysis_u3d_t(i,j,loopm_,kinconnue)      &
       +x1*ucoefmode_t(i,j,loopm_)

        analysis_v3d_t(i,j,loopm_,kinconnue)=     &
        analysis_v3d_t(i,j,loopm_,kinconnue)      &
       +x1*vcoefmode_t(i,j,loopm_)

  425  continue

! Eventuellement analyser egalement la pression 3D
      if(flag_3dwaves_pharmonics==1) then !pmxpmx>
       do loopm_=1,countmodemax ; do j=1,jmax ; do i=1,imax
        analysis_p3d_t(i,j,loopm_,kinconnue)=     &
        analysis_p3d_t(i,j,loopm_,kinconnue)      &
       +x1*pcoefmode_t(i,j,loopm_)
       enddo ; enddo ; enddo
      endif                               !pmxpmx>

  426 continue

! Calcul progressif du produit de la matrice principale par sa transposee: !18-07-10
      do kinconnue=1,2*kmaxtide+1
      do kequation=1,2*kmaxtide+1

      k1=int(real(kequation+1)/2.)
      k2=int(real(kinconnue+1)/2.)

      time1=frqtide(k1)*(elapsedtime_now-ti0tide(k1))            &  !16-04-11
            +v0tide(k1)+utide(k1,1)

      time2=frqtide(k2)*(elapsedtime_now-ti0tide(k2))            &  !16-04-11
            +v0tide(k2)+utide(k2,1)

      analysis3dmatrix(kequation,kinconnue)=                         &
      analysis3dmatrix(kequation,kinconnue)+                         &
       (mod(kequation,2)*cos(time1)+(1.-mod(kequation,2))*sin(time1))*ftide(k1,1) &
      *(mod(kinconnue,2)*cos(time2)+(1.-mod(kinconnue,2))*sin(time2))*ftide(k2,1)


      enddo
      enddo

      endif             !11111111111111111111111111>


      if(case_==2) then !22222222222222222222222222>

! Final step at the end of the simulation:
      if(par%rank==0)write(6,*) &
       '3dwaves analysis based on ',analysiscount_3dwaves &
      ,' analysed fields'

      call allocate_systeme(1,1,1,2*kmaxtide+1)            !21/01/08

      do loopm_=1,countmodemax

      do kinconnue=1,2*kmaxtide+1
      do kequation=1,2*kmaxtide+1
       tab3(kequation,kinconnue)=analysis3dmatrix(kequation,kinconnue)
      enddo
      enddo

      ksecu=0                                                          !09/07/04
      nbinco=2*kmaxtide+1

      if(par%rank==0)write(6,'(a,i0)') &
      'Analyse harmonique pour u3d mode ',loopm_
      do j3=1,jmax
      do i3=1,imax
      if(  mask_t(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtide
            do k1=1,2
               kinconnue=kinconnue+1
               tab4(kinconnue)=analysis_u3d_t(i3,j3,loopm_,kinconnue)
            enddo
            enddo
            ktide=kmaxtide+1 ! background freq nulle
            kinconnue=kinconnue+1
            tab4(kinconnue)=analysis_u3d_t(i3,j3,loopm_,kinconnue)
!.....................................................................
! Au premier coup (KSECU=1) on calcule le produit par la transposée
! Au deuxieme coup (KSECU.GE.2) on relie le resultat:
                           ksecu=ksecu+1                               !09/07/04
                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04
!.....................................................................
                           do ktide=1,kmaxtide
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            u3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)
                            u3dmode_sin_t(i3,j3,loopm_,ktide)=tab7(k2)
                           enddo
                           ktide=kmaxtide+1 ! background freq nulle
                           k1=(ktide-1)*2+1
                           u3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)

      endif                            !££££££££££££££££££££££££>
      enddo ! fin de boucles
      enddo ! sur i3 et j3


      if(par%rank==0)write(6,'(a,i0)') &
      'Analyse harmonique pour v3d mode ',loopm_
      do j3=1,jmax
      do i3=1,imax
      if(  mask_t(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtide
            do k1=1,2
               kinconnue=kinconnue+1
               tab4(kinconnue)=analysis_v3d_t(i3,j3,loopm_,kinconnue)
            enddo
            enddo
            ktide=kmaxtide+1 ! background freq nulle
            kinconnue=kinconnue+1
            tab4(kinconnue)=analysis_v3d_t(i3,j3,loopm_,kinconnue)
!.....................................................................
! Au premier coup (KSECU=1) on calcule le produit par la transposée
! Au deuxieme coup (KSECU.GE.2) on relie le resultat:
                           ksecu=ksecu+1                               !09/07/04
                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04
!.....................................................................
                           do ktide=1,kmaxtide
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            v3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)
                            v3dmode_sin_t(i3,j3,loopm_,ktide)=tab7(k2)
                           enddo
                           ktide=kmaxtide+1 ! background freq nulle
                           k1=(ktide-1)*2+1
                           v3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)

      endif                            !££££££££££££££££££££££££>
      enddo ! fin de boucles
      enddo ! sur i3 et j3


      if(flag_3dwaves_pharmonics==1) then !pmxpmx> !14-04-16

      if(par%rank==0)write(6,'(a,i0)') &
      'Analyse harmonique pour p3d mode ',loopm_
      do j3=1,jmax
      do i3=1,imax
      if(  mask_t(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtide
            do k1=1,2
               kinconnue=kinconnue+1
               tab4(kinconnue)=analysis_p3d_t(i3,j3,loopm_,kinconnue)
            enddo
            enddo
            ktide=kmaxtide+1 ! background freq nulle
            kinconnue=kinconnue+1
            tab4(kinconnue)=analysis_p3d_t(i3,j3,loopm_,kinconnue)
!.....................................................................
! Au premier coup (KSECU=1) on calcule le produit par la transposée
! Au deuxieme coup (KSECU.GE.2) on relie le resultat:
                           ksecu=ksecu+1                               !09/07/04
                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04
!.....................................................................
                           do ktide=1,kmaxtide
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            p3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)
                            p3dmode_sin_t(i3,j3,loopm_,ktide)=tab7(k2)
                           enddo
                           ktide=kmaxtide+1 ! background freq nulle
                           k1=(ktide-1)*2+1
                           p3dmode_cos_t(i3,j3,loopm_,ktide)=tab7(k1)

      endif                            !££££££££££££££££££££££££>
      enddo ! fin de boucles
      enddo ! sur i3 et j3

      endif                               !pmxpmx>

      enddo ! fin de boucle sur loopm_
! LIBERER LA MEMOIRE DU SOLVER:
      call allocate_systeme(2,0,0,0) ! desallouer

      endif             !22222222222222222222222222>

      end subroutine modeanalysis_harmonics

!......................................................................

      subroutine modeanalysis_pmodeprojection !14-04-16
      implicit none
      integer loopm_
#ifdef synopsis
       subroutinetitle='modeanalysis_presmodeprojection'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(rhrefmode_t)) then !pmxpmx>
! on analyse une anomalie de pression par rapport A une pression de reference
! prise ici A la pression au depart du cycle d'analyse:
               allocate(rhrefmode_t(imax,jmax,kmax))
               do k=1,kmax ; do j=1,jmax ; do i=1,imax
                       rhrefmode_t(i,j,k)=rhp_t(i,j,k)
               enddo       ; enddo       ; enddo
      endif                               !pmxpmx>

      if(elapsedtime_now<elapsedtime_end-timewindow_3dwaves)return

      id_prs=0 ! attention anyv3d(:,:,:,1 a 4 et de 6 a 7) est pris par les coef de l'EOS
! verifier la disponibilite de anyv3d(:,:,:,id_prs)
      if(anyv3d(-1,-1,0,id_prs)==-9999.) &
      stop 'Err anyv3d id_prs not available'

! Calculer la pression hydro A partir de l'anomalie de densitE:
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,1)=0.
      enddo ; enddo
      do k=kmax,1,-1 ; do j=1,jmax ; do i=1,imax
!$ Hydrostatic pressure:
       x2=grav*(rhp_t(i,j,k)-rhrefmode_t(i,j,k))*dsig_t(i,j,k)*h_w(i,j) ! dsig*h pour oter la composante 2D
       xy_t(i,j,1)=xy_t(i,j,1)+x2               
       anyv3d(i,j,k,id_prs)=xy_t(i,j,1)-0.5*x2 
      enddo ; enddo ; enddo

! Cette routine s'appelle depuis la routine de gradient
      do loopm_=1,countmodemax
       do j=1,jmax
       do i=1,imax

        pcoefmode_t(i,j,loopm_)=0.
        do k=1,kmax
         pcoefmode_t(i,j,loopm_)=pcoefmode_t(i,j,loopm_)  &
               +dz_t(i,j,k,1)*anyv3d(i,j,k,id_prs)        &
         *uv_wmode_t(i,j,k,loopm_)
        enddo

       enddo
       enddo
      enddo

! Marquer "disponibles" les tableaux generiques 
      anyv3d(-1,-1,0,id_prs)=0.

      end subroutine modeanalysis_pmodeprojection

!......................................................................

      end module module_modeanalysis
