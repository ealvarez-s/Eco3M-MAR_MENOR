










      subroutine initial_bio
!______________________________________________________________________
! S model
! release S.26 - last update: 12-02-17
!______________________________________________________________________
      use module_principal
      use module_parallele !#mpi
      use module_biobc
      use module_dust 
      use module_depnit
      use module_depammo
      use ModuleComBlock
      use UserDeclaration
      use ModuleDeclaration
      implicit none
      integer kobc


!______________________________________________________________________
! Version date      Description des modifications
!         29/08/01: si strada n'est pas utilisé on ne lit pas le
!                   notebook_tracer. Ce qui évite de "planter" le modèle
!                   si l'utilisateur n'a pas pris le soin de bien
!                   construire son notebook_tracer
!         05/09/01: modification de l'ordre de lecture des profondeurs
!                   min et max de l'immersion initial du traceur (l'odre
!                   croissant est plus logique)
!         25/09/01: amenagements sur CALL DATE_TO_KOUNT
!         14/12/01: Bienvenue à ITIMEBIO
!         17/12/01: amenagements pour cas forward pour economie de RAM
!         22/11/02: bug dans asselin_upd quand ITIMEBIO=1. LE CAS EST STOPPE
!         24/12/02: mise en place filiere restart pour la biologie
!         26/12/02: remise en service du schema forward
!         27/01/03: debuggeur taille memoire
!         11/04/06: modif sur dimensions de RIVER_BIO
!         18/04/06: fonctions compatibles avec double precision
!         10/05/06: fonctions compatibles avec double precision
!         21/06/06: option model_ bio 1D
! 2009.3  01-10-09: Utilisation des nouveaux facteurs d'echelle verticale
! 2010.8  03-05-10  terminologie
!         09-05-10  seul le proc 0 ecrit interrupteur
! 2010.10 23-06-10  hot_restart renommé dyn_restart et bio_restart
! S.26    28-03-13  downscaling bio
!         09-12-14  initpelagic 'coquille vide'
!         12-02-17  wseb dimension verticale
!______________________________________________________________________



!*******************************************************************************
! LECTURE du notebook_tracer
! DEBUT:
!*******************************************************************************

      itimebio=dim_timebio ! 0=leapfrog 1=forward                      !26/12/02
! section debug debut:                                                 !26/12/02
      if(dim_timebio.ne.0.and.dim_timebio.ne.1) then !%%%%%%>
      write(6,*)'erreur dans parameter,'
      write(6,*)'dim_timebio doit etre egal a 0 ou a 1'
      stop ' dans initial_bio.f!'
      endif                                          !%%%%%%>
      if(itimebio.ne.dim_timebio) then !§§§§§§§§§§§§§§§§§§§§>
      write(6,*)'erreur dans parameter,'
      write(6,*)'mettre dim_timebio égal à ',itimebio,' puis relancer.'
      stop ' dans initial_bio.f!'
      endif                            !§§§§§§§§§§§§§§§§§§§§>
! section debug fin.                                                   !26/12/02



!cccccIF(ITIMEBIO.EQ.1)stop 'OPTION NON VALIDE. BUG DANS ASSELIN_UPD'   !26/12/02
! Explication du bug: si itimebio=1, couple_modes est appelé ds strada !22/11/02
! et vlx_x(i,j,k,1) est modifié. Cette modification interfere à tort   !22/11/02
! dans les calculs de asselin_upd.F. A corriger à l'avenir si on veut  !22/11/02
! remettre en route le schema forward.                                 !22/11/02


      open(unit=3,file=nomfichier(12)) ! lecture notebook_bio

!     print*,'notebook_bio:',nomfichier(12)

      read(3,'(a80)')texte90
         if(texte90(80:80).ne.'1')then
         write(6,*)'erreur de lecture du notebook_tracer section 1'
         stop 'dans initial_bio.f'
         endif
      read(3,*)imodelbio

! si strada n'est pas utilisé on ne lit pas le reste du notebook
      if(imodelbio.eq.0) then
          close(3)
          return                                         !29/08/01
      endif

      read(3,*)vbmax
      vbmax=vbmax*imodelbio
      call biology_allocate
      biobc_type(:)=1

      read(3,*)dti_bio
!      dti_bio=dti_bio*3600.

      read(3,'(a80)')texte90
         if(texte90(80:80).ne.'2')then
         write(6,*)'erreur de lecture du notebook_tracer section 2'
         stop 'dans initial_bio.f'
         endif
      read(3,*) ! ligne pour rien
      k0=-10
      do 44 vb=1,vbmax
      read(3,*)wsed(1,vb),texte250    !28-03-13!12-02-17
      biovarname(vb)=trim(texte250)
! securité au cas où les étourdis se trompent de signe:
      wsed(1,vb)=-abs(wsed(1,vb)) ! on est donc bien sur que wsed < 0
   44 continue

      read(3,'(a80)')texte90
         if(texte90(80:80).ne.'3')then
         write(6,*)'erreur de lecture du notebook_tracer section 3'
         write(6,*)'car le nbre de vitesses de chute ne correspond'
         write(6,*)'probablement pas au nbre de variable. verifier!'
         stop 'dans initial_bio.f'
         endif
      read(3,*) ! ligne pour rien
      read(3,*) ! ligne pour rien

      do 60 kr=1,nriver
      read(3,*) ! ligne pour rien
!      print*,'river kr:',kr
      do    vb=1,vbmax
!      print*,'vb:',vb
      read(3,*)river_bio(vb,kr,1)                                      !11/04/06
               river_bio(vb,kr,2)=river_bio(vb,kr,1)                   !11/04/06
               river_bio(vb,kr,0)=river_bio(vb,kr,1)                   !11/04/06
      enddo

      read(3,*)ild_to_sd(kr)
      read(3,*)bioriv_info(kr)
      read(3,'(A)')bioriv_file(kr)
      read(3,*)bioriv_date(1,kr)  &
              ,bioriv_date(2,kr)  &
              ,bioriv_date(3,kr)  &
              ,bioriv_date(4,kr)  &
              ,bioriv_date(5,kr)  &
              ,bioriv_date(6,kr)

!      print*,'bioriv kr',kr

   60 continue


      read(3,'(a80)')texte90
         if(texte90(80:80).ne.'4')then
         write(6,*)'erreur de lecture du notebook_tracer section 4'
         write(6,*)'car le nbre de fleuve fois le nbre '
         write(6,*)'de concentration ne correspond pas au nbre de'
         write(6,*)'fleuves et au nbre de variables. verifier!'
         stop 'dans initial_bio.f'
         endif
      read(3,*) ! ligne pour rien
      read(3,*)ksomax
      call biology_allocate_socard

      read(3,*) ! ligne pour rien

      do 88 kso=1,ksomax

      read(3,*) ! ligne pour rien
      read(3,*)i1,x1,x2
      if(i1.eq.1) then  ! @@@@@@@@@@>
        latit1=x1*pi/180.
        longi1=x2*pi/180.
!       call latlon_to_ij(1)
        call latlon_to_ij('glb')
        socard(1,kso)=min(un*(imax+1),max(zero,deci))
        socard(2,kso)=min(un*(jmax+1),max(zero,decj))
      else              ! @@@@@@@@@@>
        socard(1,kso)=x1
        socard(2,kso)=x2
      endif             ! @@@@@@@@@@>

      if(  mask_t(nint(socard(1,kso)),nint(socard(2,kso)),kmax+1).eq.0)then
      write(6,*)'attention source n°',kso,' en terre'
      write(6,*)'ses indices i j sont'                                  &
       ,nint(socard(1,kso)),nint(socard(2,kso))
      stop 'donc dans initial_bio.f'
      endif

      read(3,*)i1,x1,x2
      if(i1.eq.1) then  ! ---------->
        i2=nint(socard(1,kso))
        j2=nint(socard(2,kso))
        k1= 9999
        k2=-9999
        do 104 k=1,kmax
         if( depth_t(i2,j2,k).ge.x1)k1=min0(k1,k)                      !05/09/01
         if( depth_t(i2,j2,k).le.x2)k2=max0(k2,k)                      !05/09/01
  104   continue
        k1=min0(k1,k2)
        socard(3,kso)=min0(kmax,max0(1,k1))
        socard(4,kso)=min0(kmax,max0(1,k2))
      else              ! ---------->
        socard(3,kso)=x1
        socard(4,kso)=x2
      endif             ! ---------->

      read(3,*)socard(5,kso),socard(6,kso)

      read(3,*)socard(7,kso)

      read(3,*)i1,i2,i3,i4,i5,i6
      call datetokount
      x3=xdtk_out                                                      !25/09/01
      socard(8,kso)=x3

      read(3,*)i1,i2,i3,i4,i5,i6
      call datetokount
      x3=xdtk_out                                                      !25/09/01
      socard(9,kso)=x3

! DESCRIPTION DE SOCARD
! 1er argument:
!              SOCARD(1,KSO): indice i de la source
!              SOCARD(2,KSO): indice j de la source
!              SOCARD(3,KSO): indice kmin de la source
!              SOCARD(4,KSO): indice kmax de la source
!              SOCARD(5,KSO): numero classe min de la variable
!              SOCARD(6,KSO): numero classe max de la variable
!              SOCARD(7,KSO): concentration à la source
!              SOCARD(8,KSO): iteration d'apparition de la source
!              SOCARD(9,KSO): iteration d'extinction de la source
! 2eme argument: numéro de la source
   88 continue

      read(3,'(a80)')texte90
         if(texte90(80:80).ne.'5')then
         write(6,*)'erreur de lecture du notebook_bio section 5'
         write(6,*)'car le nbre d''infos sur les sources ne'
         write(6,*)'correspond pas au nbre de sources. verifier!'
         stop 'dans initial_bio.f'
         endif

! Definir les conditions aux limites
! reset:
      do vb=1,vbmax
       biobc_type(vb)=0
      enddo
! lecture
      read(3,*) ! une ligne pour rien
      k0=0
      do vb=1,vbmax
       read(3,*)biobc_type(vb)
       if(biobc_type(vb)==7) k0=1
      enddo

!     if(biobc_type(initrate)==7) then
      if(k0==1) then

      do kobc=1,6
      read(3,*)obc_bio_info(kobc)
      read(3,'(A)')obc_bio_file(kobc)
      read(3,*)obc_bio_date(1,kobc)  &
              ,obc_bio_date(2,kobc)  &
              ,obc_bio_date(3,kobc)  &
              ,obc_bio_date(4,kobc)  &
              ,obc_bio_date(5,kobc)  &
              ,obc_bio_date(6,kobc)
      enddo
!     print*,'passe 1'
      read(3,'(A)')obc_bio_filedepth
!      print*,'obc_bio_filedepth',obc_bio_filedepth
      endif
!      print*,'flag_nemoffline',flag_nemoffline
      if(flag_nemoffline==1) then
      read(3,'(a)')filename_runoff
      endif

      read(3,*)
! ajout pour DIC debut
      read(3,*)bioatm_info
      read(3,'(A)')bioatm_file
      read(3,*)bioatm_date(1)  &
              ,bioatm_date(2)  &
              ,bioatm_date(3)  &
              ,bioatm_date(4)  &
              ,bioatm_date(5)  &
              ,bioatm_date(6)
! ajout pour DIC fin
      read(3,*)idust
!      print*,'idust:',idust
      read(3,*)ichoicedust
!      print*,'ichoicedust:',ichoicedust
      if(ichoicedust.le.2) then
      ndust=3 ! aladin2015 (P. Nabat) ou regcm2015 (M. Mallet) (stage aude)
      elseif(ichoicedust.eq.3) then
      ndust=2 ! aladin2018
      endif
      ndepnit=3
      ndepammo=3
!      print*,'ndust',ndust
      do k=1,ndust
      read(3,'(a)')dust_file(k)
!      print*,'dust_file(k)',dust_file(k),k
!      read(3,*)
      enddo
      read(3,*)dustpfr
      read(3,*)psolublefr
      read(3,*)flag_dust_average
      read(3,*)idepnit
      do k=1,ndepnit
      read(3,'(A)')depnit_file(k)
      enddo
      read(3,*)flag_depnit_average
      read(3,*)idepammo
      if(par%rank==1) print*,'idepammo=',idepammo
      do k=1,ndepammo
      read(3,'(A)')depammo_file(k)
      if(par%rank==1) print*,'depammo_file(k)',k,depammo_file(k)  
      enddo
      read(3,*)flag_depammo_average
!      print*,'flag_depammo_average',flag_depammo_average

      read(3,'(A80)')texte90
      if(texte90(80:80).NE.'6')THEN
         write(6,*)'texte90:',texte90
         write(6,*)'Erreur de lecture du notebook_bio section 6'
         write(6,*)'car le nbre d''infos sur les conditions limites'
         write(6,*)'ne correspond pas au nbre de variables. Verifier!'
         stop 'dans initial_bio.F'
      endif

! Par defaut model_ 3D:
      i1dbio=0
      iptbio=1
      jptbio=1
      read(3,*,end=239)i1dbio,iptbio,jptbio                            !21/06/06

      read(3,'(A80)')texte90
      if(texte90(80:80).NE.'7')THEN
         WRITE(6,*)'Erreur de lecture du notebook_bio section 7'
         stop 'dans initial_bio.F'
         write(6,*)'Texte90',texte90
      endif

! Par defaut
      mbio1=1
      mbio2=imax
      nbio1=1
      nbio2=jmax

      read(3,*)gridbio
      if(gridbio.eq.1) THEN
      read(3,*)mbio1_glob,mbio2_glob,nbio1_glob,nbio2_glob
      if(par%nbdom>1)then
      if(mbio1_glob/=1.and.mbio2_glob/=iglb.and.nbio1_glob/=1.and.nbio2_glob/=jglb) then
      write(*,*)'init_bio arret car job parallele et pour le moment '
      write(*,*)'mbio1, mbio2, nbio1, nbio2 dans ce cas doivent etre égaux à 1,imax,1,jmax'
      stop 'init_bio parallele1'
      endif
      endif
      mbio2=min(mbio2,imax)
      nbio2=min(nbio2,jmax)
      if(par%rank==0)then
       write(6,*)'Attention !!!!!!!!!!  Grille du modele bio :'
       write(6,*)'mbio1,mbio2,nbio1,nbio2',mbio1,mbio2,nbio1,nbio2
      endif
      endif

  239 continue
      close(3)

      if(imodelbio.eq.1) then !*************************************>


      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=4,file=trim(tmpdirname)//'messages',position='append')
      write(4,*)'-----------------------------------------------------'
      write(4,*)'subroutine initial_bio:'
      write(4,*)
      write(4,*)'model_ traceurs strada (bio, sed..) activé'
      write(4,*)'nbre de variables: ',vbmax
      do kso=1,ksomax
      write(4,*)'source test dispersion n°',kso
      write(4,*)'i & j                 : ',socard(1,kso),socard(2,kso)
      write(4,*)'kmin et max           : ',socard(3,kso),socard(4,kso)
      write(4,*)'variables min et max  : ',socard(5,kso),socard(6,kso)
      write(4,*)'concentration         : ',socard(7,kso)
      write(4,*)'iterations min et max : ',socard(8,kso),socard(9,kso)
      enddo
      do vb=1,vbmax
         WRITE(4,*)'OBC Variable ',vb,' est de type: ',biobc_type(vb)
         ENDDO

      close(4)
      endif                !#mpi-->>-->                       !09-05-10
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10

      endif                   !*************************************>
!*******************************************************************************
! LECTURE du notebook_bio
! FIN.
!*******************************************************************************
!*******************************************************************************
! INITIALISATION DU MODELE PELAGIQUE ET DES VARIABLES BIO_T
! DEBUT
!*******************************************************************************

      call biobc_driver !19-02-24 ! moved up to test for Mar Menor
      call InitPelagic
      !call biobc_driver !28-03-13 
      if(idust==1)       call dust_driver(1)  
      if(idepnit==1)   call depnit_driver(1)
      if(idepammo==1) call depammo_driver(1)

!*******************************************************************************
! INITIALISATION DU MODELE PELAGIQUE
! FIN
!*******************************************************************************


!*******************************************************************************
! INITIALISATION DU MODELE BENTHIQUE
! DEBUT
!*******************************************************************************

!#ifdef 1
      CALL InitBenthic
!      print*,'Ibenth',IBenthic     
!#endif

!*******************************************************************************
! INITIALISATION DU MODELE BENTHIQUE
! FIN
!*******************************************************************************


!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! DEBUT:
!*******************************************************************************

      if(initial.eq.0) then
      do vb=1,vbmax
      do i=1,imax
      do j=1,jmax
        fluxbio_w(i,j,vb,1)=zero
        fluxbio_w(i,j,vb,2)=zero
      enddo
      enddo
      enddo

!     do 173 vb=1,vbmax
!     do 173  k=1,kmax
!     do 173  i=1,imax
!     do 173  j=1,jmax
!     bio_t(i,j,k,vb)=1.
! 173 continue
      endif

! Cas d'un départ d'une simulation précédente:

      if(initial.eq.1) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      write(6,*)
      write(6,*)'....................................................'
      write(6,*)'dans initial_bio.f: êtes vous sur que le fichier'
      write(6,*)'restart du model_ bio (chanelbio_in) existe?'
      write(6,*)'sinon commenter l`appel a hot_restart(11)'

!     call hot_restart(11)                                             !24/12/02
      call bio_restart('r') !23-06-10

      write(6,*)'la lecture du fichier restart bio s`est bien passée,'
      write(6,*)'on continue...'
      write(6,*)'....................................................'

      endif

      if(initial.eq.2) then
      call bio_restart('r') !23-06-10


            print*,'dans initial_bio apres bio_restart'

!       call InitNutDoxy

!       IBenthic=2



       

! Mise a zero des diagnostiques relus dans le bio_restart
! Specifique au cas "initial=2"
       do vb=1,vbmax

        EXPORT_200(vb)=0.
        SUM_EXPORT_200(vb)=0.
        EXPORT_200(vbmax+vb)    =0.
        SUM_EXPORT_200(vbmax+vb)=0.
        EXPORT_200(2*vbmax+vb)    =0.
        SUM_EXPORT_200(2*vbmax+vb)=0.
        EXPORT_plat(vb)     =0.
        SUM_EXPORT_plat(vb) =0.
        EXPORTGDL_200(vb)=0.
        SUM_EXPORTGDL_200(vb)=0.
        EXPORTGDL_200(vbmax+vb)    =0.
        SUM_EXPORTGDL_200(vbmax+vb)=0.
        EXPORTGDL_200(2*vbmax+vb)    =0.
        SUM_EXPORTGDL_200(2*vbmax+vb)=0.
        EXPORTGDL_plat(vb)     =0.
        SUM_EXPORTGDL_plat(vb) =0.
        EXPORTDYF_200(vb)=0.
        SUM_EXPORTDYF_200(vb)=0.
        EXPORTDYF_200(vbmax+vb)    =0.
        SUM_EXPORTDYF_200(vbmax+vb)=0.
        EXPORTDYF_200(2*vbmax+vb)    =0.
        SUM_EXPORTDYF_200(2*vbmax+vb)=0.
        EXPORT_200_1000m(vb)=0.
        SUM_EXPORT_200_1000m(vb)=0.

! Apport rivieres
        BIL_FLEUVES(vb)=0.
          do kr=1,nriver
           BIL_FLEUVE(kr,vb)=0.
          enddo

! Export frontieres
        SUM_OB_SOUTH(vb)=0.
        SUM_OB_WEST(vb)=0.
        SUM_OB_EAST(vb)=0.
        FLUX_OB_SOUTH(vb)=0.
        FLUX_OB_WEST(vb)=0.
        FLUX_OB_EAST(vb)=0.

          do k=1,kmax
           SUM_OB_SOUTH_VERT(k,vb)=0.
          SUM_OB_WEST_VERT(k,vb)=0.
           FLUX_OB_SOUTH_VERT(k,vb)=0.
           FLUX_OB_WEST_VERT(k,vb)=0.
!           do J=nbio1+1,nbio2
!            SUM_OB_WEST_PT(J-nbio1+1,k,vb)=0.
!            FLUX_OB_WEST_PT(J-nbio1+1,k,vb)=0.
!           enddo
!           do I=mbio1+1,mbio2
!            SUM_OB_SOUTH_PT(I-mbio1+1,k,vb)=0.
!            FLUX_OB_SOUTH_PT(I-mbio1+1,k,vb)=0.
!           enddo
          enddo

      do k=1,kmax-1
      do j=1,jmax
      do i=1,imax
      if(bio_t(i,j,k,vb).le.1.d-10) bio_t(i,j,k,vb)=1.d-10
      enddo
      enddo
      enddo

      enddo !vb


      TPS_PPB=0
      TPS_STRADA=0
      do k=1,4
      SUMPPBTOTAL(k)   = 0.
      PPBTOTAL(k)      = 0.
      SUMNPBTOTAL(k)   = 0.
      NPBTOTAL(k)      = 0.
      SUMRPBTOTAL(k)   = 0.
      RPBTOTAL(k)      = 0.
      SUMEXPTOTAL(k)   = 0.
      SUMEXPTOTAL(k+4) = 0.
      EXPTOTAL(k)      = 0.
      EXPTOTAL(k+4)    = 0.
      SUMPPBLIG(k)     = 0.
      SUMPPBDYF(k)     = 0.
      SUMPPBGDL(k)     = 0.
      PPBLIG(k)        = 0.
      PPBDYF(k)        = 0.
      PPBGDL(k)        = 0.
      SUMNETPPBTOTAL(k)= 0.
      NETPPBTOTAL(k)   = 0.
      SUMRESPPHYTOTOTAL(k)= 0.
      RESPPHYTOTOTAL(k)= 0.
      SUMNETPPBLIG(k)  = 0.
      SUMNETPPBDYF(k)  = 0.
      SUMNETPPBGDL(k)  = 0.
      NETPPBLIG(k)     = 0.
      NETPPBDYF(k)     = 0.
      NETPPBGDL(k)     = 0.
      SUMPPBTOTAL_plat(k) = 0.
      PPBTOTAL_plat(k) = 0.
      enddo

      DO K=1,3
      SUMPPBGDL_platcl(K)   = 0.
      PPBGDL_platcl(K)      = 0.
      ENDDO

      SUMPPBGDL_plat   = 0.
      SUMRPBGDL_plat   = 0.
      SUMNPBGDL_plat   = 0.
      SUMNETPPBGDL_plat= 0.
      PPBGDL_plat      = 0.
      RPBGDL_plat      = 0.
      NPBGDL_plat      = 0.
      NETPPBGDL_plat   = 0.
      SUMGBACGDL_plat  = 0.
      SUMUPTBACTDONGDL_plat    = 0.
      SUMEXCHETERONGDL_plat    = 0.
      SUMREMPOCGDL_plat = 0.
      SUMREMPONGDL_plat = 0.
      SUMLOSTPOCGDL_plat = 0.
      SUMLOSTPONGDL_plat = 0.
      SUMRESPPHYTOGDL_plat = 0.
      SUMRESPTOTGDL_plat = 0.
      SUMRESPZOOGDL_plat = 0.
      SUMRPBTOTAL(5)   = 0.
      RPBTOTAL(5)      = 0.
      SUMNITRifTOTAL   = 0.
      NITRifTOTAL      = 0.
      SUMRESPTOTTOTAL  = 0.
      RESPTOTTOTAL     = 0.
      SUMEXUCTOTTOTAL  = 0.
      EXUCTOTTOTAL     = 0.
      SUMGBACTOTAL     = 0.
      SUMGBACDYF       = 0.
      SUMGBACGDL       = 0.
      SUMGRAZCTOTAL    = 0.
      GRAZCTOTAL       = 0.

! Exportation flux
      SUMEXPDYF   = 0.
      SUMEXPGDL   = 0.

      SUM_EXPORTC_BOT = 0.
      SUM_EXPORTN_BOT = 0.
      SUM_EXPORTP_BOT = 0.
      SUM_EXPORTSI_BOT = 0.
      EXPORTC_BOT = 0.
      EXPORTN_BOT = 0.
      EXPORTP_BOT = 0.
      EXPORTSI_BOT = 0.
      SUM_SEDC    = 0.
      SUM_SEDP    = 0.
      SUM_SEDN    = 0.
      SUM_SEDSI   = 0.
          SEDC    = 0.
          SEDP    = 0.
          SEDN    = 0.
          SEDSI   = 0.

! mise a zero des diagnostics 2D

      TPS_PPB_2D=0
      TPS_STRADA_2D=0
      do j1=1,jmax ! debut boucle sur j
      do i1=1,imax ! debut boucle sur i
      do k=1,4
        PPB2D(i1,j1,k)   =0.
        NPB2D(i1,j1,k)   =0.
        RPB2D(i1,j1,k)   =0.
        CHL2D(i1,j1,k)   =0.
      enddo
        do k=1,27 ! 
        EXP2D(i1,j1,k)   =0.
        enddo
!        NITRif2D(i1,j1)  =0.
        RESP2D(i1,j1)    =0.
        EXUCTOT2D(i1,j1) =0.
        RPB2D(i1,j1,5)   =0.
        NO3EFFLUX2D(i1,j1) =0. ! Alex 
        NH4EFFLUX2D(i1,j1) =0.
        PEFFLUX2D(i1,j1) =0. 
        SIEFFLUX2D(i1,j1)=0. 
        CDEPO(i1,j1)     =0. 
        NDEPO(i1,j1)     =0. 
        PDEPO(i1,j1)     =0. 
        SIDEPO(i1,j1)    =0. 
        ExcbactNH4NTOPLAYER2D(i1,j1)=0.
        ExcZooAmmoNTOPLAYER2D(i1,j1)=0.
        ExcZooPO4PTOPLAYER2D(i1,j1)=0.
        ExcbactPO4PTOPLAYER2D(i1,j1)=0.
        ExuSiTOPLAYER2D(i1,j1)=0.
        UptNitTOPLAYER2D(i1,j1)=0.
        UptAmmoTOPLAYER2D(i1,j1)=0.
        UptPTOPLAYER2D(i1,j1)=0.
        UptSiTOPLAYER2D(i1,j1)=0.
        NitrifTOPLAYER2D(i1,j1)=0.
        RemSMOPSiTOPLAYER2D(i1,j1)=0.
        RemLMOPSiTOPLAYER2D(i1,j1)=0.
        ExcbactNH4NINT2D(i1,j1)=0.
        ExcZooAmmoNINT2D(i1,j1)=0.
        ExcZooPO4PINT2D(i1,j1)=0.
        ExcbactPO4PINT2D(i1,j1)=0.
        ExuSiINT2D(i1,j1)=0.
        UptNitINT2D(i1,j1)=0.
        UptAmmoINT2D(i1,j1)=0.
        UptPINT2D(i1,j1)=0.
        UptSiINT2D(i1,j1)=0.
        NitrifINT2D(i1,j1)=0.
        RemSMOPSiINT2D(i1,j1)=0.
        RemLMOPSiINT2D(i1,j1)=0.
        ExcbactNH4NDEEP2D(i1,j1)=0.
        ExcZooAmmoNDEEP2D(i1,j1)=0.
        ExcZooPO4PDEEP2D(i1,j1)=0.
        ExcbactPO4PDEEP2D(i1,j1)=0.
        ExuSiDEEP2D(i1,j1)=0.
        UptNitDEEP2D(i1,j1)=0.
        UptAmmoDEEP2D(i1,j1)=0.
        UptPDEEP2D(i1,j1)=0.
        UptSiDEEP2D(i1,j1)=0.
        NitrifDEEP2D(i1,j1)=0.
        RemSMOPSiDEEP2D(i1,j1)=0.
        RemLMOPSiDEEP2D(i1,j1)=0.
 
!claude verifier qu'il faut bien cidessous sum_ devant export caroline n'a pas
        SUM_EXPORTC_BOT_2D(i1,j1)    =0.
        SUM_EXPORTN_BOT_2D(i1,j1)    =0.
        SUM_EXPORTP_BOT_2D(i1,j1)    =0.
        SUM_EXPORTSI_BOT_2D(i1,j1)    =0.
      enddo    ! fin de boucle i1
      enddo    ! fin de boucle j1

      endif !if initial.eq.2

!      call InitDOxy



!*******************************************************************************
! INITIALISATION DES TABLEAUX:
! FIN.
!*******************************************************************************

      return
      end
