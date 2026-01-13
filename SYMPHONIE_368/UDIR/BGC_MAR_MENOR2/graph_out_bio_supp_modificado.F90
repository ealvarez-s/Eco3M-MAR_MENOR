        subroutine graph_out_bio_supp
!______________________________________________________________________
! S model
! release S26 - last update: 02-04-15
!______________________________________________________________________
!......................................................................
! Version date      Description des modifications:
! S.26    02-11-14  mise en service
!         08-11-14  filval deplace
!         20-12-14  masquage avec masque de surface
!         02-04-15  modif sur nfmpi a la demande de malek
!......................................................................
      use module_principal ; use module_parallele 
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_bio_supp'
       subroutinedescription='Produces a netcdf output file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! si variable tracer non selectionnee dans notebook_rivers sortir:
       if(grh_out_var(54)/=1)return

! returns texte250 the bio netcdf file:
       call graph_out_bio_file_name_supp 

! writes bio_t in the bionetcdf file:
       call graph_out_bio_write_var_supp 

      end  subroutine graph_out_bio_supp

!.......................................................................

      subroutine graph_out_bio_file_name_supp
      use module_principal ; use module_parallele 
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_bio_file_name_supp'
       subroutinedescription='returns texte250 the bio netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!--------------------------------------------
! DEFINIR LE NOM DU FICHIER:
!.....calcul de la date                                            !27/01/03
      call elapsedtime2date(elapsedtime_now,i5,i6,i7,i3,i2,i1)     !01-04-13
      write(texte80(1),                            &
       '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')  &
       i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

!.....écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04


! Nom des fichiers de variables
      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //'.ECO3M-S.LA.tlse_supp.nc'

!.....puis afficher la date à l'écran:                                     !27/01/03
      if(par%rank==0) then !00000>
       write(6,*)'--------------------------------'
       write(6,*)'Biogeochemical netcdf file date='
       write(6,'(a33)')texte80(1)(1:33)
       write(6,'(a,a)')'File name: ',trim(texte250)
      endif                !00000>

      end subroutine graph_out_bio_file_name_supp

!.......................................................................

      subroutine graph_out_bio_write_var_supp
      use module_principal
      use module_parallele !#MPI
      use module_modeanalysis
      use module_my_outputs
      use ModuleDeclaration 
      use pnetcdf
      implicit none
      integer ichoix,looproc_,looprocmax_,loopm_,loop1_,loop2_
      character*6 type_
      integer(kind=MPI_OFFSET_KIND) len_
      integer indice_bio(37),kk,np1,np2      !bio
      real  &
       orga_picop_400(imax,jmax), &
       orga_nanop_400(imax,jmax), &
       orga_microp_400(imax,jmax), &
       orga_nanoz_400(imax,jmax), &
       orga_microz_400(imax,jmax), &
       orga_mesoz_400(imax,jmax)
      real  &
       orga_phy_sw(imax,jmax), &
       orga_phy_iw(imax,jmax), &
       orga_phy_dw(imax,jmax), &
       orga_hetero_sw(imax,jmax), &
       orga_hetero_iw(imax,jmax), &
       orga_hetero_dw(imax,jmax), &
       orga_poc_sw(imax,jmax), &
       orga_poc_iw(imax,jmax), &
       orga_poc_dw(imax,jmax), &
       orga_doc_sw(imax,jmax), &
       orga_doc_iw(imax,jmax), &
       orga_doc_dw(imax,jmax)

  
!      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='graph_out_bio_write_var_supp'
       subroutinedescription='writes bio_t in the bionetcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(bio_t))return

      filval=-9999.

      do loop_netcdf=0,1

      count_netcdfvar=0

      if(loop_netcdf==0) then  !§§§§§§§>
       status=nfmpi_create(par%comm2d,texte250 ,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      else                     !§§§§§§§>
       status=nfmpi_open(par%comm2d,texte250 ,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      endif                    !§§§§§§§>
      if(status/=0)stop ' stop module_offline erreur 1'

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      call kount_to_date(0)                          !16-11-09  time origin corresponds to kount=0

      write(texte80(2),'(a14,i4,a1,i2,4(a1,i2))')                    & !units
      'seconds since ',i5,'-',i6,'-',i7,' ',i3,':',i2,':',i1
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0'
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0'
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'
      texte80(8)=texte80(2)(14:33)                                      ! time_origin : kount=0
      texte80(9)='gregorian'                                            ! calendar
      texte80(3:4)='time'                                               ! long_name
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'
      call netcdf_main('_t')

! write variables dimensions:
      call graph_out_trueaxis_nfmpi

      texte80(10)='none'

!.....
      goto 200
! write bio_t
      do vb=1,vbmax !VBVBVBVBVBVB>

!        data indice_bio /iZooNanoC,iZooMicroC,iZooMesoC           &
!                     ,iSyneC,iNanoC,iDiaC                         &
!                      ,iSyneN,iNanoN,iDiaN                         &
!                      ,iSyneP,iNanoP,iDiaP,iDiaSi                  &
!                      ,iSyneChl,iNanoChl,iDiaChl                  &
!                      ,iBactC                                       &
!                      ,iSMOPC                                       &
!                      ,iSMOPN,iSMOPP,iSMOPSI,ISMOPCHL                &
!                      ,iLMOPC                                        &
!                       ,iLMOPN,iLMOPP,iLMOPSI                &
!                      ,iMODC                                        &
!                      ,iMODN,iMODP                           &
!                      ,iNitrate,iAmmonium,iPhosphate,iSilice       &
!                      ,iOxygen/   ! MODIF LE 26 03 2013
!        data indice_bio /iZooNanoC,iZooMicroC,iZooMesoC           &
!                     ,iSyneC,iNanoC,iDiaC/
!         data indice_bio /iZooNanoC,iZooMicroC,iZooMesoC           &
!                      ,iSyneC,iNanoC,iDiaC                         &
!                       ,iBactC, &
!                          iSMOPC,iLMOPC, &
!                           iMODC,iMODN,iMODP, &
!                          iNitrate,iAmmonium,iPhosphate, &
!                          iSilice,iOxygen/           
!                     ,iSyneC,iNanoC,iDiaC/
!        data indice_bio /iNitrate,iAmmonium,iPhosphate, &
!                          iSilice,idic,ialkalinity/



       if (vbmax.gt.0) then !<><><><><>->

!     do vb=1,vbmax !- vb loop ->
!      do kk=1,23 !- vb loop ->
!      do kk=1,36 !- vb loop ->
!      do kk=1,15 !- vb loop ->
!     do kk=1,34 !- vb loop ->
!      do kk=30,30 !- vb loop ->

!      goto 200

!     do kk=1,17 !- vb loop ->
!      do kk=1,6 !- vb loop ->

!          vb=indice_bio(kk)

             if(vb.eq.1)texte80(1)='zoonanoc'
             if(vb.eq.2)texte80(1)='zoomicroc'
             if(vb.eq.3)texte80(1)='zoomesoc'
             if(vb.eq.4)texte80(1)='synec'
             if(vb.eq.5)texte80(1)='synen'
             if(vb.eq.6)texte80(1)='synep'
             if(vb.eq.7)texte80(1)='synechl'
             if(vb.eq.8)texte80(1)='nanoc'
             if(vb.eq.9)texte80(1)='nanon'
             if(vb.eq.10)texte80(1)='nanop'
             if(vb.eq.11)texte80(1)='nanochl'
             if(vb.eq.12)texte80(1)='diac'
             if(vb.eq.13)texte80(1)='dian'
             if(vb.eq.14)texte80(1)='diap'
             if(vb.eq.15)texte80(1)='diachl'
             if(vb.eq.16)texte80(1)='diasi'
             if(vb.eq.17)texte80(1)='bactc'
             if(vb.eq.18)texte80(1)='smopc'
             if(vb.eq.19)texte80(1)='smopn'
             if(vb.eq.20)texte80(1)='smopp'
             if(vb.eq.21)texte80(1)='smopchl'
             if(vb.eq.22)texte80(1)='smopsi'
             if(vb.eq.23)texte80(1)='lmopc'
             if(vb.eq.24)texte80(1)='lmopn'
             if(vb.eq.25)texte80(1)='lmopp'
             if(vb.eq.26)texte80(1)='lmopsi'
             if(vb.eq.27)texte80(1)='modc'
             if(vb.eq.28)texte80(1)='modn'
             if(vb.eq.29)texte80(1)='modp'
             if(vb.eq.30)texte80(1)='nitrate'
             if(vb.eq.31)texte80(1)='ammonium'
             if(vb.eq.32)texte80(1)='phosphate'
             if(vb.eq.33)texte80(1)='silice'
             if(vb.eq.34)texte80(1)='oxygen'
             if(vb.eq.35)texte80(1)='odu'
             if(vb.eq.36)texte80(1)='dic'
             if(vb.eq.37)texte80(1)='talk'

      filval=-9999. !08-11-14

      anyvar2d=-9999 ! pour supprimer les valeurs aberrantes ncview
      anyvar3d=-9999 !

      if(loop_netcdf==1) then !=======>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax

        anyvar3d(i,j,k)=bio_t(i,j,k,vb)*mask_t(i,j,k)          & !20-12-14
                                    +(1-mask_t(i,j,k))*filval
!         anyvar3d(i,j,k)=bio_t(i,j,k,vb)*mask_t(i,j,kmax)          & !20-12-14
!                                     +(1-mask_t(i,j,kmax))*filval
                if(vb==ioxygen) then 
                        anyvar3d(i,j,k)=bio_t(i,j,max0(k,kmin_w(i,j)),vb) &
                        /(rhp_t(i,j,max0(k,kmin_w(i,j)))+rho)*1000. !pour avoir des umol/kg     
                endif
 
          if(vb==idic) &
          anyvar3d(i,j,k)=bio_t(i,j,k,vb)/(rhp_t(i,j,k)+rho)*1000.
          if(vb==ialkalinity) &
           anyvar3d(i,j,k)=( bio_t(i,j,k,vb)         &
                            +bio_t(i,j,k,iAmmonium)  &
                            -bio_t(i,j,k,iNitrate)   &
                            -bio_t(i,j,k,iPhosphate)) &
                            /(rhp_t(i,j,k)+rho)*1000.


        enddo ; enddo ; enddo


      endif                   !=======>
!      write(texte80(1),'(a,i0)')'bio',vb
      texte80(2)='none'         ! variable ; units
      write(texte80(3),'(a,i0)')'biogeochemical_variable_number_',vb
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')

      endif !<><><><><>->

      enddo         !VBVBVBVBVBVB>

!      endif
 
  200 continue

!! EA: added pH and PAR to supp output
! PH---------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=0,jmax+1 !30-07-14
         do i=0,imax+1
!          if (mask_t(i,j,kmax  ).ne.0) then
!          anyvar3d(i,j,k)=spH(i,j,max0(k,kmin_w(i,j))) 
          if (mask_t(i,j,k  ).ne.0) then
          anyvar3d(i,j,k)=sPH(i,j,k)
          else
           anyvar3d(i,j,k)=-9999.
          endif
         enddo
         enddo
         enddo
      endif                  !=======>
      texte80(1)='pH'     ; texte80(2)='unit'    ! variable;units
      write(texte80(3),'(a)')'pH'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')

!--------------------------------------------------------

! PAR---------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=0,jmax+1
         do i=0,imax+1
          if (mask_t(i,j,k  ).ne.0) then
          anyvar3d(i,j,k)=sPAR(i,j,k)
          else
           anyvar3d(i,j,k)=-9999.
          endif
         enddo
         enddo
         enddo
      endif                  !=======>
      texte80(1)='PAR'     ; texte80(2)='w/m2'      ! variable;units
      write(texte80(3),'(a)')'PAR'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')

!--------------------------------------------------------

! sal---------------------------------------------------
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=0,jmax+1 !30-07-14
!         do i=0,imax+1
!!          if (mask_t(i,j,kmax  ).ne.0) then
!!          anyvar3d(i,j,k)=spH(i,j,max0(k,kmin_w(i,j))) 
!          if (mask_t(i,j,k  ).ne.0) then
!          anyvar3d(i,j,k)=sal_t(i,j,k,1)
!          else
!           anyvar3d(i,j,k)=-9999.
!          endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======>
!      texte80(1)='sal'     ; texte80(2)='unit'                       !variable;units
!      call netcdf_main('_t')

!--------------------------------------------------------

! sal---------------------------------------------------
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=0,jmax+1 !30-07-14
!         do i=0,imax+1
!!          if (mask_t(i,j,kmax  ).ne.0) then
!!          anyvar3d(i,j,k)=spH(i,j,max0(k,kmin_w(i,j))) 
!          if (mask_t(i,j,k  ).ne.0) then
!          anyvar3d(i,j,k)=tem_t(i,j,k,1)
!          else
!           anyvar3d(i,j,k)=-9999.
!          endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======>
!      texte80(1)='tem'     ; texte80(2)='unit'                       !variable;units
!      call netcdf_main('_t')
!
!--------------------------------------------------------

!! 18/03/2017 Test change chltot Alex
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax  !

!          anyvar3d(i,j,k)= (bio_t(i,j,k,idiachl) &
!                         +bio_t(i,j,k,isynechl) &
!                         +bio_t(i,j,k,inanochl))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='chl_tot_prueba'     ; texte80(2)='mg/m3'
!      write(texte80(3),'(a)')'chlo tot'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!---------------------------------------------------------
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_2d>0) then
!          anyvar3d(i,j,k)= nitrif3d(i,j,k)/tps_ppb_2d
!         endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='nitrif3d'     ; texte80(2)='mmolN/m3/day'
!      write(texte80(3),'(a)')'nitrif3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!------------------------------------------------------------
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_2d>0) then
!          anyvar3d(i,j,k)= uptnit3d(i,j,k)/tps_ppb_2d
!         endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!           
!      texte80(1)='uptnit3d'     ; texte80(2)='mmolN/m3/day'
!      write(texte80(3),'(a)')'uptnit3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!---------------------------------------------------------------      
! Alex ajout variables 2D d'apres graph_out_bon.F90 dans PELAGIC/ 
!      if(loop_netcdf==1) then !=======>
!       k=kmax
!           do j=-1,jmax+2
!           do i=-11,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!         do j=1,jmax !30-07-14
!         do i=1,imax  
!
!          anyvar2d(i,j)= (bio_t(i,j,k,idiachl) &
!                         +bio_t(i,j,k,isynechl) &
!                         +bio_t(i,j,k,inanochl))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='chl_tot'     ; texte80(2)='mg/m3'
!      write(texte80(3),'(a)')'chlo tot'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!----------------------------------------------------------------------


!! EA: 3d fluxes in supp output
!------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= ppb3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='gpp3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'gpp3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!-------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= resp3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
         enddo
      endif                  !=======> 
!
      texte80(1)='resp3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'resp3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!
!----------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= nitrif3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='nitrif3d'     ; texte80(2)='mgN/m3/day'
      write(texte80(3),'(a)')'nitrif3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= uptnit3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
         enddo
      endif                  !=======> 
           
      texte80(1)='uptnit3d'     ; texte80(2)='mgN/m3/day'
      write(texte80(3),'(a)')'uptnit3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!--------------------------------------------------------------------------      
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= ppbp3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='gpppico3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'gpppico3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= ppbn3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 
!
      texte80(1)='gppnano3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'gppnano3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!----------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= ppbm3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 
!
      texte80(1)='gppmicro3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'gppmicro3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!-----------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= resppp3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 
!
      texte80(1)='respPpico3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respPpico3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!--------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= resppn3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 
!
      texte80(1)='respPnano3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respPnano3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!---------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= resppm3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 
!
      texte80(1)='respPmicro3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respPmicro3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!-------------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= uptnitMol3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='uptnitMol3d'     ; texte80(2)='mmolN/m3/day'
      write(texte80(3),'(a)')'uptnitMol3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!----------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= nitrifMol3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='nitrifMol3d'     ; texte80(2)='mmolN/m3/day'
      write(texte80(3),'(a)')'nitrifMol3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!-----------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= respP3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='respP3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respP3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!--------------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= respZ3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='respZ3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respZ3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!---------------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax !30-07-14
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= respB3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='respB3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'respB3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!----------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax !30-07-14
         do i=1,imax
         if(tps_ppb_3d>0) then
          anyvar3d(i,j,k)= exuc3d(i,j,k)/tps_ppb_3d
         endif
         enddo
         enddo
        enddo
      endif                  !=======> 

      texte80(1)='exuc3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
      write(texte80(3),'(a)')'exuc3d'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------      

!*****************************************************************
! production primaire
! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
!! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=ppb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
            enddo
        endif
        texte80(1)='ppb' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'gross primary production'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'

        call netcdf_main('_t')
!        print*,'ppb tracee'
!------------------------------------------------------------------------
!! Total N, P and Si standing stocks
      if(loop_netcdf==1) then !=====
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=TOT_COL_N_t(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
            enddo
        endif
        texte80(1)='tot_N' ; texte80(2)='mmolN/m2'
      write(texte80(3),'(a)')'total N'
      texte80(4)=texte80(3)
     texte80(5)='TYX' ; texte80(7)='real'

        call netcdf_main('_t')
!--------------------------------------------------------------------------
      if(loop_netcdf==1) then !=====
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=TOT_COL_P_t(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
            enddo
        endif
        texte80(1)='tot_P' ; texte80(2)='mmolP/m2'
      write(texte80(3),'(a)')'total P'
      texte80(4)=texte80(3)
     texte80(5)='TYX' ; texte80(7)='real'

        call netcdf_main('_t')
!--------------------------------------------------------------------------
      if(loop_netcdf==1) then !=====
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=TOT_COL_Si_t(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
            enddo
        endif
        texte80(1)='tot_Si' ; texte80(2)='mmolSi/m2'
      write(texte80(3),'(a)')'total Si'
      texte80(4)=texte80(3)
     texte80(5)='TYX' ; texte80(7)='real'

        call netcdf_main('_t')
!--------------------------------------------------------------------------
!! Euphotic layer depth
      if(loop_netcdf==1) then !=====
! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=EuphoticLayerDepth_t(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
            enddo
        endif
        texte80(1)='z_eu' ; texte80(2)='m'
      write(texte80(3),'(a)')'depth of the euphotic layer'
      texte80(4)=texte80(3)
     texte80(5)='TYX' ; texte80(7)='real'

        call netcdf_main('_t')

!!*****************************************************************
! production primaire nette
! mise ?zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
               anyvar2d(i,j)=netppb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='netppb' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'net primary production'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
        call netcdf_main('_t')
!!        print*,'ppb tracee'
!!*****************************************************************
! production nouvelle
! mise  zero pour masquer couronne peripherique qd MBIO1 NE 1 etc...

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=npb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='npb' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'production npb'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
      call netcdf_main('_t')

!*****************************************************************
! production regeneree

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=rpb2d(i,j,1)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif

        texte80(1)='rpb' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'production rpb'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
         
        call netcdf_main('_t')
!!        print*,'rpb tracee'
!!!***************************************************************
!! nitrification
        
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=nitrif2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
          texte80(1)='nitrif' ; texte80(2)='mgN/m2/d'
      write(texte80(3),'(a)')'nitrification'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
          
        call netcdf_main('_t')
!!        print*,'nitrif tracee'
!!!*****************************************************************
!! respiration total

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=resp2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
          texte80(1)='resp' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'total respiration'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
          
        call netcdf_main('_t')
!       print*,'resp tracee'
!!!*****************************************************************
! respiration zoo

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=rzoo2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
          texte80(1)='respZoo' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'zoo respiration'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
          
        call netcdf_main('_t')
!!       print*,'resp tracee'
!!!*****************************************************************
! respiration bacteria

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=rbac2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
          texte80(1)='respBact' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'bacteria respiration'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
          
        call netcdf_main('_t')
!!       print*,'resp tracee'
!-------------------------------------------------------------------------
! respiration phyto

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=rphyto2d(i,j)/tps_ppb_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='respPhyto' ; texte80(2)='mgC/m2/d'
      write(texte80(3),'(a)')'phyto respiration'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
         
        call netcdf_main('_t')
!!        print*,'rpb tracee'        
!!!*****************************************************************
!!! CDepo  

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax

              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=CDepo(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='CDepo' ; texte80(2)='mmol/m2/d'
        call netcdf_main('_t')
!        print*,'nitrif tracee'
!!!*****************************************************************
!!! NDepo  

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax

              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=NDepo(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='NDepo' ; texte80(2)='mmol/m2/d'
        call netcdf_main('_t')
!        print*,'nitrif tracee'
!!!*****************************************************************
!! PDepo  

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=PDepo(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='PDepo' ; texte80(2)='mmol/m2/d'
        call netcdf_main('_t')
!!        print*,'nitrif tracee'
!!!*****************************************************************
!! SiDepo  

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=SiDepo(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
          endif
        texte80(1)='SiDepo' ; texte80(2)='mmol/m2/d'
        call netcdf_main('_t')
!!        print*,'nitrif tracee'
!!!!*****************************************************************
!
!!!!*****************************************************************
!!! Benthic flux NO3
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=NO3efflux2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='NO3efflux2d' ; texte80(2)='mmolN/m2/d'
        call netcdf_main('_t')
!!        print*,'export_mopc tracee'
!!!!*****************************************************************
!!! Benthic flux NH4
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=NH4efflux2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='NH4efflux2d' ; texte80(2)='mmolN/m2/d'
        call netcdf_main('_t')
!!        print*,'export_mopc tracee'
!!!*****************************************************************
!!! Benthic flux P
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=Pefflux2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='Pefflux2d' ; texte80(2)='mmolP/m2/d'
        call netcdf_main('_t')
!!!        print*,'export_mopc tracee'
!!!!*****************************************************************
!!! Benthic flux Si
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax          
              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=Siefflux2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='Siefflux2d' ; texte80(2)='mmolSi/m2/d'
        call netcdf_main('_t')
!!!        print*,'export_mopc tracee'
!!!!*****************************************************************
!! Benthic flux DIC        
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax

              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=DICefflux2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='DICefflux2d' ; texte80(2)='mmolC/m2/d'
        call netcdf_main('_t')
!!!*****************************************************************
!! Benthic flux O2        
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax

              if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
                anyvar2d(i,j)=O2influx2d(i,j)/tps_benth_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
         endif
        texte80(1)='02influx2d' ; texte80(2)='mmol/m2/d'
        call netcdf_main('_t')
!!!*****************************************************************
!
!
!!!!*****************************************************************
!!------------------------------------------------------------------
!! Atmospheric deposition
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
!               anyvar2d(i,j)=dd_w(i,j,1)
!               anyvar2d(i,j)=fluxbio_w(i,j,iphosphate,2)
                anyvar2d(i,j)=AtmDepPho(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='AtmDepPhoDust' ; texte80(2)='mmol/m2/s'
        call netcdf_main('_t')
!
!!------------------------------------------------------------------
!! Atmospheric deposition
!      if(loop_netcdf==1) then !=======>
!           do j=1,jmax
!           do i=1,imax
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!! reduction des sorties a la zone de calcul
!           do j=1,jmax
!           do i=1,imax
!
!!             do j=nbio1,nbio2
!!             do i=mbio1,mbio2
!              if (mask_t(i,j,kmax+1)==1) then
!!               anyvar2d(i,j)=dd_w(i,j,1)
!!                 anyvar2d(i,j)=fluxbio_w(i,j,imodp,2)
!               anyvar2d(i,j)=AtmDepDOP(i,j)
!              else
!                anyvar2d(i,j)=-9999.
!              endif
!             enddo
!             enddo
!        endif
!        texte80(1)='AtmDepDOP' ; texte80(2)='mmol/m2/s'
!        call netcdf_main('_t')
!
!!------------------------------------------------------------------
!! Atmospheric deposition
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax

!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
!               anyvar2d(i,j)=dd_w(i,j,1)
!               anyvar2d(i,j)=fluxbio_w(i,j,initrate,2)
               anyvar2d(i,j)=AtmDepNit(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='AtmDepNit' ; texte80(2)='mmol/m2/s'
        call netcdf_main('_t')
!
!!------------------------------------------------------------------
!! Atmospheric deposition
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo

! reduction des sorties a la zone de calcul
           do j=1,jmax
           do i=1,imax

!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
!               anyvar2d(i,j)=dd_w(i,j,1)
!               anyvar2d(i,j)=fluxbio_w(i,j,iammonium,2)
                anyvar2d(i,j)=AtmDepAmmo(i,j)
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='AtmDepAmmo' ; texte80(2)='mmol/m2/s'
        call netcdf_main('_t')
!
!!------------------------------------------------------------------
!! Atmospheric deposition
!      if(loop_netcdf==1) then !=======>
!           do j=1,jmax
!           do i=1,imax
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!! reduction des sorties a la zone de calcul
!           do j=1,jmax
!           do i=1,imax
!
!!             do j=nbio1,nbio2
!!             do i=mbio1,mbio2
!              if (mask_t(i,j,kmax+1)==1) then
!!               anyvar2d(i,j)=dd_w(i,j,1)
!!                anyvar2d(i,j)=fluxbio_w(i,j,imodn,2)
!               anyvar2d(i,j)=AtmDepDON(i,j)
!              else
!                anyvar2d(i,j)=-9999.
!              endif
!             enddo
!             enddo
!        endif
!        texte80(1)='AtmDepDON' ; texte80(2)='mmol/m2/s'
!        call netcdf_main('_t')

!------------------------------------------------------------------

!! Atmospheric deposition
!      if(loop_netcdf==1) then !=======>
!           do j=1,jmax
!           do i=1,imax
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!! reduction des sorties a la zone de calcul
!           do j=1,jmax
!           do i=1,imax
!
!!             do j=nbio1,nbio2
!!             do i=mbio1,mbio2
!              if (mask_t(i,j,kmax+1)==1) then
!!               anyvar2d(i,j)=dd_w(i,j,1)
!!                anyvar2d(i,j)=fluxbio_w(i,j,imodc,2)
!               anyvar2d(i,j)=AtmDepDOC(i,j)
!              else
!                anyvar2d(i,j)=-9999.
!              endif
!             enddo
!             enddo
!        endif
!        texte80(1)='AtmDepDOC' ; texte80(2)='mmol/m2/s'
!        call netcdf_main('_t')

!!!*****************************************************************
!!!!*****************************************************************
! O2_Flux
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1.and.tps_o2f_2d>0) then
                anyvar2d(i,j)=O2flux2d(i,j)/tps_o2f_2d
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='O2Flux' ; texte80(2)='mmolO2/m2/s'
        call netcdf_main('_t')
!!!!*****************************************************************
! O2_Flux
      call equation_of_state_potzref_jmfwg('grp',1)

      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                anyvar2d(i,j)=Oxysat(i,j)/(rhp_t(i,j,kmax)+rho)*1000.
              else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='Oxygen sat' ; texte80(2)='µmol kg-1'
        call netcdf_main('_t')

!!!!*****************************************************************
! pCO2 water
      if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                 anyvar2d(i,j)=pCO2W(i,j)
                else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='pCO2w' ; texte80(2)='µatm'
        call netcdf_main('_t')

!!!!*****************************************************************
! pCO2 air
     if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1) then
                 anyvar2d(i,j)=atm_bio(1)
!                 anyvar2d(i,j)=rho
                else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='pCO2a' ; texte80(2)='µatm'
        call netcdf_main('_t')

!*****************************************************************
! CO2 air-sea flux
     if(loop_netcdf==1) then !=======>
           do j=1,jmax
           do i=1,imax
            anyvar2d(i,j)=-9999.
           enddo
           enddo
           do j=1,jmax
           do i=1,imax
!             do j=nbio1,nbio2
!             do i=mbio1,mbio2
              if (mask_t(i,j,kmax+1)==1.and.tps_asf_2d>0) then
                 anyvar2d(i,j)=asf2d(i,j)/tps_asf_2d
                else
                anyvar2d(i,j)=-9999.
              endif
             enddo
             enddo
        endif
        texte80(1)='CO2_Flux' ; texte80(2)='mmolC/m2/s'
        call netcdf_main('_t')

!------------------------------------------------------------------
!!*****************************************************************
! phychl surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! instantannée
!          anyvar2d(i,j)= (bio_t(i,j,k,idiachl) &
!                         +bio_t(i,j,k,isynechl) &
!                         +bio_t(i,j,k,inanochl)   )*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= chlsurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='CHL'     ; texte80(2)='mgChl/m3'
      write(texte80(3),'(a)')'CHL'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')
!
!!*****************************************************************
! phychl surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= chlpsurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='CHLpico'     ; texte80(2)='mgChl/m3'
      write(texte80(3),'(a)')'CHLpico'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')
!
!!*****************************************************************
! phychl surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= chlnsurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='CHLnano'     ; texte80(2)='mgChl/m3'
      write(texte80(3),'(a)')'CHLnano'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')
!
!!*****************************************************************
! phychl surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= chlmsurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='CHLmicro'     ; texte80(2)='mgChl/m3'
      write(texte80(3),'(a)')'CHLmicro'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')
!
!!*****************************************************************
!! nitrate surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= nitsurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='surfnit'     ; texte80(2)='mmol/m3'
      write(texte80(3),'(a)')'surfnit'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')
!
!!*****************************************************************
!! phosphate surface
      if(loop_netcdf==1) then !=======>

           do j=-1,jmax+2
           do i=-1,imax+2
            anyvar2d(i,j)=-9999.
           enddo
           enddo
     k= kmax
         do j=1,jmax !30-07-14
         do i=1,imax
! moyennée sur la période de sorties
         if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
          anyvar2d(i,j)= phosurf(i,j) / tps_ppb_2d
         endif
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='surfpho'     ; texte80(2)='mmol/m3'
      write(texte80(3),'(a)')'surfpho'
      texte80(4)=texte80(3)
      texte80(5)='TYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')

!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!     k= kmax
!
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar2d(i,j)= bio_t(i,j,k,initrate) *mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='NO3'     ; texte80(2)='mmolN/m3'
!      write(texte80(3),'(a)')'NO3'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!!*****************************************************************
!! phosphate surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!     k= kmax
!
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar2d(i,j)= bio_t(i,j,k,iphosphate) *mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='PO4'     ; texte80(2)='mmolP/m3'
!      write(texte80(3),'(a)')'PO4'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!
!!*****************************************************************
!! oxygen surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!     k= kmax
!
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar2d(i,j)= bio_t(i,j,k,ioxygen) *mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='O2'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'O2'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!!*****************************************************************
! oxygen surface
!      call equation_of_state_potzref_jmfwg('grp',1)
!
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!
!    k= kmax
!
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar2d(i,j)= (bio_t(i,j,k,ioxygen) *mask_t(i,j,k) &
!                        /(rhp_t(i,j,k)+rho)*1000.)   &  !pour avoir des umol/kg
!                                    +(1-mask_t(i,j,k))*filval
!
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='O2kg'     ; texte80(2)='µmol/kg'
!      write(texte80(3),'(a)')'O2kg'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!*****************************************************************
! PP net surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!     k= kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!!        if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
!!          anyvar2d(i,j)= SurfNetPP2D(i,j) * mask_t(i,j,k) &
!!                         +(1-mask_t(i,j,k))*filval
!
!!         if (mask_t(i,j,kmax+1)==1) then
!#ifdef bilanbio
!          anyvar2d(i,j)= zone1_mask(i,j) 
!#endif
!!         endif
!         enddo
!         enddo
!      endif                  !=======> 
!!     texte80(1)='PP'     ; texte80(2)='mgC/m3'
!      texte80(1)='zone1_mask_bilan'     ; texte80(2)='-'
!!     write(texte80(3),'(a)')'PP'
!      write(texte80(3),'(a)')'zone1_mask_bilan'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!      call netcdf_main('_t')

! PP net surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!     k= kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!        if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
!          anyvar2d(i,j)= SurfNetPP2D(i,j) * mask_t(i,j,k) &
!                         +(1-mask_t(i,j,k))*filval
!
!         if (mask_t(i,j,kmax+1)==1) then
!#ifdef bilanbio
!          anyvar2d(i,j)= zone2_mask(i,j)
!#endif

!         endif

!         enddo
!         enddo
!      endif                  !=======> 

!!     texte80(1)='PP'     ; texte80(2)='mgC/m3'
!      texte80(1)='zone2_mask_bilan'     ; texte80(2)='-'
!!     write(texte80(3),'(a)')'PP'
!      write(texte80(3),'(a)')'zone2_mask_bilan'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!      call netcdf_main('_t')

!! PP net surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!     k= kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax

!        if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
!          anyvar2d(i,j)= SurfNetPP2D(i,j) * mask_t(i,j,k) &
!                         +(1-mask_t(i,j,k))*filval

!         if (mask_t(i,j,kmax+1)==1) then
!#ifdef bilanbio
!          anyvar2d(i,j)= zone3_mask(i,j)
!#endif

!         endif

!         enddo
!         enddo
!      endif                  !=======> 

!!     texte80(1)='PP'     ; texte80(2)='mgC/m3'
!      texte80(1)='zone3_mask_bilan'     ; texte80(2)='-'
!!     write(texte80(3),'(a)')'PP'
!      write(texte80(3),'(a)')'zone3_mask_bilan'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!      call netcdf_main('_t')
! PP net surface
!      if(loop_netcdf==1) then !=======>
!
!           do j=-1,jmax+2
!           do i=-1,imax+2
!            anyvar2d(i,j)=-9999.
!           enddo
!           enddo
!     k= kmax
!
!         do j=1,jmax !30-07-14
!         do i=1,imax
!        if (mask_t(i,j,kmax+1)==1.and.tps_ppb_2d>0) then
!          anyvar2d(i,j)= SurfNetPP2D(i,j) * mask_t(i,j,k) &
!                         +(1-mask_t(i,j,k))*filval
!
!!         if (mask_t(i,j,kmax+1)==1) then
!#ifdef bilanbio
!          anyvar2d(i,j)= zone4_mask(i,j)
!#endif
!
!         endif

!         enddo
!         enddo
!      endif                  !=======> 

!!     texte80(1)='PP'     ; texte80(2)='mgC/m3'
!      texte80(1)='zone4_mask_bilan'     ; texte80(2)='-'
!!     write(texte80(3),'(a)')'PP'
!      write(texte80(3),'(a)')'zone4_mask_bilan'
!      texte80(4)=texte80(3)
!      texte80(5)='TYX' ; texte80(7)='real'
!      call netcdf_main('_t')

!!!!*****************************************************************
!!*****************************************************************

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

!     status=nfmpi_put_att_text(ncid,nf_global,'Production',27        &
      status=nfmpi_put_att_text(ncid,nf_global,'Production',int(30,8) & !02-04-15
      ,'ECO3M-S LEGOS-SIROCCO Toulouse')
!     if(status/=0)stop ' Err1 graph_out_bio nfmpi_put_att_text'

! Definition of variables: done.
      status=nfmpi_enddef(ncid)

      endif                  !>>>>>>>>>>>>>>>>>>>


      status=nfmpi_close(ncid)

      enddo  ! fin de boucle sur loop_netcdf

!******************************************************************************
!! impression var 2d fin
! remise a zero
!        print*,'remise a zero des var 2D'
      tps_ppb_2d=0
      tps_benth_2d=0
      tps_strada_2d=0
      tps_asf_2d=0
      tps_o2f_2d=0
        do j=nbio1,nbio2 ! debut boucle sur j
        do i=mbio1,mbio2 ! debut boucle sur i
        do k=1,4
          ppb2d(i,j,k)=0.
          npb2d(i,j,k)=0.
          rpb2d(i,j,k)=0.
          netppb2d(i,j,k)=0.
!         chl2d(i,j,k)=0. 
        enddo
        do k=1,22
          exp2d(i,j,k)=0.
        enddo
          nitrif2d(i,j)=0.
          resp2d(i,j)=0.
          rzoo2d(i,j)=0.
          rbac2d(i,j)=0.
!          exuctot2d(i,j)=0.
          rpb2d(i,j,5)=0.
          CDepo(i,j)=0.
          NDepo(i,j)=0.
          PDepo(i,j)=0.
          SiDepo(i,j)=0.
          NO3efflux2d(i,j)=0.
          NH4efflux2d(i,j)=0.
          Pefflux2d(i,j)=0.
          Siefflux2d(i,j)=0.
          DICefflux2d(i,j)=0.
          O2influx2d(i,j)=0.
          graz2d(i,j)=0.
          ExcbactNH4NTOPLAYER2D(i,j)=0.
          ExcZooAmmoNTOPLAYER2D(i,j)=0.
          ExcZooPO4PTOPLAYER2D(i,j)=0.
          ExcbactPO4PTOPLAYER2D(i,j)=0.
          ExuSiTOPLAYER2D(i,j)=0.
          UptNitTOPLAYER2D(i,j)=0.
          UptAmmoTOPLAYER2D(i,j)=0.
          UptPTOPLAYER2D(i,j)=0.
          UptSiTOPLAYER2D(i,j)=0.
          NitrifTOPLAYER2D(i,j)=0.
          RemSMOPSiTOPLAYER2D(i,j)=0.
          RemLMOPSiTOPLAYER2D(i,j)=0.
          ExcbactNH4NINT2D(i,j)=0.
          ExcZooAmmoNINT2D(i,j)=0.
          ExcZooPO4PINT2D(i,j)=0.
          ExcbactPO4PINT2D(i,j)=0.
          ExuSiINT2D(i,j)=0.
          UptNitINT2D(i,j)=0.   
          UptAmmoINT2D(i,j)=0.
          UptPINT2D(i,j)=0.
          UptSiINT2D(i,j)=0.
          NitrifINT2D(i,j)=0.   
          RemSMOPSiINT2D(i,j)=0.
          RemLMOPSiINT2D(i,j)=0.
          ExcbactNH4NDEEP2D(i,j)=0.
          ExcZooAmmoNDEEP2D(i,j)=0.
          ExcZooPO4PDEEP2D(i,j)=0.
          ExcbactPO4PDEEP2D(i,j)=0.
          ExuSiDEEP2D(i,j)=0.
          UptNitDEEP2D(i,j)=0.
          UptAmmoDEEP2D(i,j)=0.
          UptPDEEP2D(i,j)=0.
          UptSiDEEP2D(i,j)=0. 
          NitrifDEEP2D(i,j)=0. 
          RemSMOPSiDEEP2D(i,j)=0.
          RemLMOPSiDEEP2D(i,j)=0.
          UptBactPTOPLAYER2D(i,j)=0.
          UptBactNTOPLAYER2D(i,j)=0.
          UptBactPINT2D(i,j)=0.
          UptBactNINT2D(i,j)=0.
          UptBactPDEEP2D(i,j)=0.
          UptBactNDEEP2D(i,j)=0.
          O2flux2d(i,j)=0.
          asf2d(i,j)=0
          chlsurf(i,j)=0.
          dzsw(i,j)=0.
          dziw(i,j)=0.
          chlsw(i,j)=0.
          gppsw(i,j)=0.
          gppiw(i,j)=0.
          crpsw(i,j)=0.
          crpiw(i,j)=0.
          crzsw(i,j)=0.
          crziw(i,j)=0.
          crbsw(i,j)=0.
          crbiw(i,j)=0.
          gbsw(i,j)=0.
          gbiw(i,j)=0.
          nitrifsw(i,j)=0.
          nitrifiw(i,j)=0.
          graztotsw(i,j)=0.
          graztotsw(i,j)=0.
          grazherbisw(i,j)=0.
          grazherbisw(i,j)=0.
          nitsurf(i,j)=0.
          phosurf(i,j)=0.
          nitsw(i,j)=0.
          phosw(i,j)=0.
          gpppsw(i,j)=0.
          gppnsw(i,j)=0.
          gppmsw(i,j)=0.
          crppsw(i,j)=0.
          crpnsw(i,j)=0.
          crpmsw(i,j)=0.
          chlpsurf(i,j)=0.
          chlnsurf(i,j)=0.
          chlmsurf(i,j)=0.
          chlpsw(i,j)=0.
          chlnsw(i,j)=0.
          chlmsw(i,j)=0.
          uptnitsw(i,j)=0.
          uptammosw(i,j)=0.
          uptpsw(i,j)=0.
          uptpbsw(i,j)=0.
          TOT_COL_N_t(i,j)=0.
          TOT_COL_P_t(i,j)=0.
          TOT_COL_Si_t(i,j)=0.
          EuphoticLayerDepth_t(i,j)=0.
!          SurfNetPP2D(i,j)=0.
        do k=1,kmax
         nitrif3d(i,j,k)=0.
         uptnit3d(i,j,k)=0.
         sPH(i,j,k)=0.
         sPAR(i,j,k)=0.
         nitrifMol3d(i,j,k)=0.
         uptnitMol3d(i,j,k)=0.
         ppb3d(i,j,k)=0.
         resp3d(i,j,k)=0.
         ppbp3d(i,j,k)=0.
         ppbn3d(i,j,k)=0.
         ppbm3d(i,j,k)=0.
         respPp3d(i,j,k)=0.
         respPn3d(i,j,k)=0.
         respPm3d(i,j,k)=0.
         respP3d(i,j,k)=0.
         respZ3d(i,j,k)=0.
         respB3d(i,j,k)=0.
         exuc3d(i,j,k)=0.
        enddo
        enddo    ! fin de boucle i1
        enddo    ! fin de boucle j1
! fin remise a zero
!******************************************************************************

      end  subroutine graph_out_bio_write_var_supp
