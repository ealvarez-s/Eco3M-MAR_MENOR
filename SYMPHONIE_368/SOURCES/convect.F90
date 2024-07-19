      subroutine convect
!______________________________________________________________________
!
! S model
! release 2010.11  - last update: 23-07-10
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________

      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='convect'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! 02/06/06: mise en service par claude
! 16/06/06: correction deroulement des boucles
!...............................................................................
! Version date      Description des modifs
!         02/06/06: mise en service par claude
!         16/06/06: correction deroulement des boucles
! 2009.3  30-09-09: utilisation des nouveaux facteurs d'echelle verticale
! 2010.11 23-07-10  routine density renommee equation_of_state
!...............................................................................

      const1=-rho*alp_t
      const2= rho*alp_s

! En entrée RHP_Z contient un profil de densite fixé à priori.
! En sortie RHP_Z contient le même profil mais ajusté convectivement.

      call equation_of_state('potential density',after)                !23-07-10

! traitement de l'instabilite convective
      do j=1,jmax                                                      !16/06/06
      do i=1,imax

      if(  mask_t(i,j,kmax+1).eq.1) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



!.............................................>
 50   continue
      k1=0
      do k=kmax,kmin_w(i,j)+1,-1

       k5=k
       do k2=k-1,kmin_w(i,j),-1
        if(rhp_t(i,j,k2).gt.rhp_t(i,j,k2+1))then !--------->
         k5=k2+1
         goto 60
        endif                                          !--------->
        if(k2.eq.kmin_w(i,j))k5=kmin_w(i,j)   ! cas instable jusqu'au fond
       enddo

   60 if(k.gt.k5.and.                                                   &
         abs(rhp_t(i,j,k)-rhp_t(i,j,k5)).gt.1.e-5)then !*******>

       x3=0.
       x2=0.
       x1=0.
       x0=0.
       do k2=k,k5,-1
        x0=x0+dz_t(i,j,k2,1)                                          !30-09-09
        x1=x1+dz_t(i,j,k2,1)*rhp_t(i,j,k2)
        x2=x2+dz_t(i,j,k2,1)*tem_t(i,j,k2,2)                           !10-10-09
        x3=x3+dz_t(i,j,k2,1)*sal_t(i,j,k2,2)                           !10-10-09
       enddo
       do k2=k,k5,-1
        rhp_t(i,j,k2)  =x1/x0
        tem_t(i,j,k2,2)=x2/x0                                          !10-10-09
        sal_t(i,j,k2,2)=x3/x0
       enddo
       k1=1

      endif                                                  !*******>


      enddo

      if(k1.eq.1) goto 50    ! Test : il restait de l'instabilite à la precedente iteration


!.............................................>


      endif                          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      enddo
      enddo


! fin du traitement de l'instabilite convective
      end subroutine convect
