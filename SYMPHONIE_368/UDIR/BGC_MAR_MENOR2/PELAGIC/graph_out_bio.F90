        subroutine graph_out_bio
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
       subroutinetitle='graph_out_bio'
       subroutinedescription='Produces a netcdf output file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! si variable tracer non selectionnee dans notebook_rivers sortir:
       if(grh_out_var(54)/=1)return

! returns texte250 the bio netcdf file:
       call graph_out_bio_file_name 

! writes bio_t in the bionetcdf file:
       call graph_out_bio_write_var 

      end  subroutine graph_out_bio

!.......................................................................

      subroutine graph_out_bio_file_name
      use module_principal ; use module_parallele 
      implicit none
#ifdef synopsis
       subroutinetitle='graph_out_bio_file_name'
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

!.....�crire ann�e ann�e mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....�crire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de s�paration, "_"
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04


! Nom des fichiers de variables
      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //'.ECO3M-S.LA.tlse.nc'

!.....puis afficher la date � l'�cran:                                     !27/01/03
      if(par%rank==0) then !00000>
       write(6,*)'--------------------------------'
       write(6,*)'Biogeochemical netcdf file date='
       write(6,'(a33)')texte80(1)(1:33)
       write(6,'(a,a)')'File name: ',trim(texte250)
      endif                !00000>

      end subroutine graph_out_bio_file_name

!.......................................................................

      subroutine graph_out_bio_write_var
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
      integer indice_bio(38),kk      !bio
!      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='graph_out_bio_write_var'
       subroutinedescription='writes bio_t in the bionetcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(bio_t))return

      filval=-9999.

      do loop_netcdf=0,1

      count_netcdfvar=0

      if(loop_netcdf==0) then  !�������>
       status=nfmpi_create(par%comm2d,texte250 ,nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      else                     !�������>
       status=nfmpi_open(par%comm2d,texte250 ,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
      endif                    !�������>
      if(status/=0) then
      print*,'texte250 1',status,nf_clobber,MPI_INFO_NULL,ncid
      print*,'texte250 2',status,texte250,nf_write, MPI_INFO_NULL,ncid
      print*,'Error at rank ', nfmpi_strerror(status)
      stop ' stop module_offline erreur 1'
      endif  

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
! write bio_t
!      do vb=1,vbmax !VBVBVBVBVBVB>

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
!        data indice_bio /iSyneChl,iNanoChl,iDiaChl &
!                       ,iSyneC,iNanoC,iDiaC                         &
!                      ,iSyneN,iNanoN,iDiaN                         &
!                      ,iSyneP,iNanoP,iDiaP,iDiaSi                  &
!                       ,ibactC,iMODC,iMODN,iMODP                                        &
!                      ,iNitrate,iPhosphate,iOxygen &
!                      ,idic,iAlkalinity /  
!        data indice_bio /iZooNanoC,iZooMicroC,iZooMesoC           &
!                        ,iSyneChl,iNanoChl,iDiaChl &
!                        ,iSyneC,iNanoC,iDiaC                         &
!                      ,iOxygen / !&
        data indice_bio /iOxygen /

!                      ,idic,iAlkalinity /
!        data indice_bio /ibactC,iMODC                                        &
!                     ,iOxygen /
!      goto 200

!      call equation_of_state_potzref_jmfwg('grp',1)
       call equation_of_state('potential density',1)
!       do k=1,kmax ; do j=1,jmax ; do i=1,imax
!       rhp_t(i,j,k)=mask_t(i,j,k)*(rhp_t(i,j,k)+rho-1000.)
!       enddo ; enddo ; enddo


       if (vbmax.gt.0) then !<><><><><>->

     do vb=1,vbmax !- vb loop ->
!      do kk=1,23 !- vb loop ->
!      do kk=1,36 !- vb loop ->
!      do kk=1,15 !- vb loop ->
!     do kk=1,34 !- vb loop ->
!      do kk=30,30 !- vb loop ->
!      do kk=1,7 !- vb loop ->
!      do kk=1,1 !- vb loop ->

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
                        anyvar3d(i,j,k)=bio_t(i,j,k,vb) &
                        /(rhp_t(i,j,k)+rho)*1000.*mask_t(i,j,k)          & !20-12-14
                                    +(1-mask_t(i,j,k))*filval !pour avoir des umol/kg     
                endif
          if(vb==idic) &
          anyvar3d(i,j,k)=bio_t(i,j,k,vb)/(rhp_t(i,j,k)+rho)*1000.*mask_t(i,j,k) & !20-12-14
                                    +(1-mask_t(i,j,k))*filval
          if(vb==ialkalinity) &
           anyvar3d(i,j,k)=( bio_t(i,j,k,vb)         &
                            +bio_t(i,j,k,iAmmonium)  &
                            -bio_t(i,j,k,iNitrate)   &
                            -bio_t(i,j,k,iPhosphate)) &
                            /(rhp_t(i,j,k)+rho)*1000.*mask_t(i,j,k) & !20-12-14
                                    +(1-mask_t(i,j,k))*filval 

        enddo ; enddo ; enddo


      endif                   !=======>
!      write(texte80(1),'(a,i0)')'bio',vb
      texte80(2)='none'         ! variable ; units
      write(texte80(3),'(a,i0)')'biogeochemical_variable_number_',vb
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      call netcdf_main('_t')


      enddo         !VBVBVBVBVBVB>


      endif

!  200 continue

! 18/03/2017 Test change chltot Alex 
 
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax  
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,idiachl) &
!                         +bio_t(i,j,k,isynechl) &
!                         +bio_t(i,j,k,inanochl))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='chl_tot'     ; texte80(2)='mg/m3'
!      write(texte80(3),'(a)')'chlo tot'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')

! PH---------------------------------------------------
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=0,jmax+1 !30-07-14
!         do i=0,imax+1
!!          if (mask_t(i,j,kmax  ).ne.0) then
!!          anyvar3d(i,j,k)=spH(i,j,max0(k,kmin_w(i,j))) 
!          if (mask_t(i,j,k  ).ne.0) then
!          anyvar3d(i,j,k)=spH(i,j,k)
!          else
!           anyvar3d(i,j,k)=-9999.
!          endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======>
!      texte80(1)='pH'     ; texte80(2)='unit'                       ! variable ;units
!      call netcdf_main('_t')

!--------------------------------------------------------

!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,idiac) &
!                         +bio_t(i,j,k,isynec) &
!                         +bio_t(i,j,k,inanoc))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='phytoC'     ; texte80(2)='mg/m3'
!      write(texte80(3),'(a)')'phytoC'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
 !        do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,izoonanoc) &
!                         +bio_t(i,j,k,izoomicroc) &
!                         +bio_t(i,j,k,izoomesoc))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='zooC'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'zooC'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,ilmopc) &
!                         +bio_t(i,j,k,ismopc))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='POC'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'POC'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')


!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,initrate))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='nitrate'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'nitrate'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,iphosphate))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='phosphate'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'phosphate'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,isilice))*mask_t(i,j,k) &
!                                    +(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='silicate'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'silicate'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!
!          anyvar3d(i,j,k)= (bio_t(i,j,k,ioxygen))/(rhp_t(i,j,max0(k,kmin_w(i,j)))+rho)*1000. &
!                           *mask_t(i,j,k)+(1-mask_t(i,j,k))*filval
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='oxygen'     ; texte80(2)='mmol/m3'
!      write(texte80(3),'(a)')'oxygen'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')

!     call equation_of_state_potzref_jmfwg('grp',1)
       call equation_of_state('potential density',1)

      if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax !30-07-14
         do i=1,imax

          anyvar3d(i,j,k)=(rhp_t(i,j,k)+rho) &
                           *mask_t(i,j,k)+(1-mask_t(i,j,k))*filval
         enddo
         enddo
         enddo
      endif                  !=======> 

      texte80(1)='density'     ; texte80(2)='kg/m3'
      write(texte80(3),'(a)')'density'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
! variable ; units
      call netcdf_main('_t')

!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= ppb3d(i,j,k)/tps_ppb_3d
!!         endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='gpp3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'gpp3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= resp3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 
!
!      texte80(1)='resp3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'resp3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')
!

!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= nitrif3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!         enddo
!      endif                  !=======> 

!      texte80(1)='nitrif3d'     ; texte80(2)='mmolN/m3/day'
!      write(texte80(3),'(a)')'nitrif3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!!! variable ; units
!      call netcdf_main('_t')

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



           
!      texte80(1)='uptnit3d'     ; texte80(2)='mmolN/m3/day'
!      write(texte80(3),'(a)')'uptnit3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!! variable ; units
!      call netcdf_main('_t')

!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= ppbp3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='gpppico3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'gpppico3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= ppbn3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='gppnano3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'gppnano3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= ppbm3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='gppmicro3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'gppmicro3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= resppp3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respPpico3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respPpico3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= resppn3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respPnano3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respPnano3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= respPm3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respPmicro3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respPmicro3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptnitp3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptnitp3d'     ; texte80(2)='mmolN/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptnitp3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptnitn3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptnitn3d'     ; texte80(2)='mmolN/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptnitn3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptnitm3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptnitm3d'     ; texte80(2)='mmolN/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptnitm3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptammop3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptammop3d'     ; texte80(2)='mmolN/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptammop3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptammon3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptammon3d'     ; texte80(2)='mmolN/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptammon3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= uptammom3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='uptammom3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'uptammom3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= respZn3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respZn3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respZn3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= respZmi3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respZmi3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respZmi3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= respZme3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respZme3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respZme3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= grazn3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='grazn3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'grazn3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= grazmi3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='grazmi3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'grazmi3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= grazme3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='grazme3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'grazme3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= messyfeed3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='messyfeed3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'messyfeed3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= exuc3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='exuc3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'exuc3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= gbac3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='gbac3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'gbac3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= (resp3d(i,j,k) - &
!                            (respPp3d(i,j,k)+respPn3d(i,j,k)+respPm3d(i,j,k) &
!                            +respZn3d(i,j,k)+respZmi3d(i,j,k)+respZme3d(i,j,k))) /tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='respB3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'respB3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= mortbact3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='mortbact3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'mortbact3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')
!
!      if(loop_netcdf==1) then !=======>
!         do k=1,kmax
!         do j=1,jmax !30-07-14
!         do i=1,imax
!         if(tps_ppb_3d>0) then
!          anyvar3d(i,j,k)= remmopc3d(i,j,k)/tps_ppb_3d
!         endif
!         enddo
!         enddo
!        enddo
!      endif                  !=======> 
!
!      texte80(1)='remmopc3d'     ; texte80(2)='mgC/m3/day' !'mmolC/m3/day'
!      write(texte80(3),'(a)')'remmopc3d'
!      texte80(4)=texte80(3)
!      texte80(5)='TZYX' ; texte80(7)='real'
!      call netcdf_main('_t')


! Alex ajout variables 2D d'apres graph_out_bon.F90 dans PELAGIC/  

!*****************************************************************
!        call netcdf_main('_t')


!!!!*****************************************************************
!!*****************************************************************

      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

!     status=nfmpi_put_att_text(ncid,nf_global,'Production',27        &
      status=nfmpi_put_att_text(ncid,nf_global,'Production',int(27,8) & !02-04-15
      ,'ECO3M-S LA-SIROCCO Toulouse')
!     if(status/=0)stop ' Err1 graph_out_bio nfmpi_put_att_text'

! Definition of variables: done.
      status=nfmpi_enddef(ncid)

      endif                  !>>>>>>>>>>>>>>>>>>>


      status=nfmpi_close(ncid)

      enddo  ! fin de boucle sur loop_netcdf

!******************************************************************************
!******************************************************************************
!! impression var 2d fin
! remise a zero
!        print*,'remise a zero des var 2D'

      tps_ppb_3d=0
        do j=nbio1,nbio2 ! debut boucle sur j
        do i=mbio1,mbio2 ! debut boucle sur i
        do k=1,kmax
          ppb3d(i,j,k)=0.
          resp3d(i,j,k)=0.
          nitrif3d(i,j,k)=0.
      uptnit3d(i,j,k)=0.
      uptnitp3d(i,j,k)=0.
      uptnitn3d(i,j,k)=0.
      uptnitm3d(i,j,k)=0.
      uptammo3d(i,j,k)=0.
      uptammop3d(i,j,k)=0.
      uptammon3d(i,j,k)=0.
      uptammom3d(i,j,k)=0.
      ppbp3d(i,j,k)=0.
      ppbn3d(i,j,k)=0.
      ppbm3d(i,j,k)=0.
      respPp3d(i,j,k)=0.
      respPn3d(i,j,k)=0.
      respPm3d(i,j,k)=0.
      respZn3d(i,j,k)=0.
      respZmi3d(i,j,k)=0.
      respZme3d(i,j,k)=0.
      grazn3d(i,j,k)=0.
      grazmi3d(i,j,k)=0.
      grazme3d(i,j,k)=0.
      messyfeed3d(i,j,k)=0.
      exuc3d(i,j,k)=0.
      gbac3d(i,j,k)=0.
      mortbact3d(i,j,k)=0.
      remmopc3d(i,j,k)=0.
        enddo
        enddo
        enddo

      end  subroutine graph_out_bio_write_var
