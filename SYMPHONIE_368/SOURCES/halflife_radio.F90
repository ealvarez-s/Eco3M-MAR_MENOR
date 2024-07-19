      subroutine halflife_radio
!______________________________________________________________________
! S model
! release S.26 - last update: 01-12-16
!______________________________________________________________________

      use module_principal
      implicit none
!..............................................................................
! Version date      Description des modifications
! S.26    01-12-16  Creation de la subroutine halflife_radio
!                   Ajout d'une composante "radioactif" dans les traceurs via: 
!                   radionucleide(vb), alpha_radio(vb) et radio_coef(vb)
!..............................................................................
#ifdef synopsis
       subroutinetitle='halflife_radio'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do vb=1,vbmax 

! Pour les radionuleide s
        if (radionucleide(vb)==1) then
           radio_coef(vb)= 1.-alpha_radio(vb)*dti_fw
           bio_t(:,:,:,vb)=bio_t(:,:,:,vb)*radio_coef(vb) ! coefificient de radioactivité
        endif
      enddo

      end subroutine halflife_radio
