!----------------------------------------------------------------------
! release S26  - last update: 14-06-18
!---------------------------------------------------------------------
! version date        Description des modifications:
!         2009-04-10  Bienvenue a GTMECO,GTNECO
!         2009-09-02  Bienvenue a PAR_GATHERALL_2D_Z
! 2009.2  2009-09-03  debug
!         07-09-09    re-activer interface par_gatherall_2d_zr8
!         28-09-09    re-activer interface par_gatherall_2d_zr4
! 2009.3  05-10-09:   ajout d'un "ifdef parallele"
! 2009.3  21-10-09:   debug, par_gatherall si nombre de domaine==1
!         08-11-09:   mpi_real remplace mpi_double_precision
!         16-11-09:   debug, par_gatherall inversion orde des dimensions de globdim
!         16-11-09:   Bienvenue a initialise_parallel_manu
! 2010.8  27-04-10:   Premiere optimisation: Passage aux communications non bloquantes
!         05-05-10:   debug dans dans les cummunication non bloquante echange nord-sud
! 2010.9  28-05-10:   modif echange x2 et y2
! 2010.11 04-07-10:   Bienvenue aux echanges non-bloquant 
! 2010.13 15-10-10    mises à jours et domaines periodiques
!         29-10-10    unique -> single
! 2010.22 27-11-11    Cure d'amaigrissement : Plus d'echange Bloquant
! ??????? 03-03-12    Bienvenue aux echanges directionnels -> echange_voisin,get_type_echange
!         04-11-12    Nouveaux echanges v1,v2, u1,u2
!         29-05-13    v2, u2 dans les echanges directionnels
!         10-10-13    Bienvenue a set_position_echange, nouvelle definition des bornes des echanges
! S.26    14-02-14    mpi_neighbor_list passe dans module_principal a cause de faire_hotrestart
!         01-05-14    version qui permettra le mpi subcycling et l'eliminationdes procs
!                     inutiles
!         22-01-15    ajout get_type_echange2d_r4
!         06-02-15    Bug fix Cyril
!         07-02-15    Debug echange R (boucles internes 2:imax & 2:jmax)
!         12-04-15    ajout echange z3
!         16-02-16    norme gfortran
!         06-04-17    ajout echange 4d r4
!         19-06-17    ajout get_type_echange2d_i4_itps
!         09-09-17    mask devient integer kind=1
!         08-05-18    ajout par_gatherall_2d_i1
!         14-06-18    codage dom_c pour grand nombre de domaines
!---------------------------------------------------------------------

#ifdef parallele
      module module_parallele
      implicit none
      include  'mpif.h'
      integer :: mpi_comm_symp
      integer, parameter :: ndims=2
      integer, parameter :: ouest=1,est=2,nord=3,sud=4,haut=5,bas=6
      integer, parameter :: sudouest=7,sudest=8,nordouest=9,nordest=10,ouest2=11,est2=12 
      integer, parameter :: ouestest=1,nordsud=2
      integer, parameter :: nbvoisin=10
!      integer,dimension(10),parameter :: tagSend = (/ 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000 /)
!      integer,dimension(10),parameter :: tagRecv = (/ 2200,2100,2400,2300,2600,2500,3000,2900,2800,2700 /)
       integer,dimension(12),parameter :: tagSend = (/ 210,220,230,240,250,260,270,280,290,300, 310, 320 /)
       integer,dimension(12),parameter :: tagRecv = (/ 220,210,240,230,260,250,300,290,280,270, 320, 310 /)
!      integer,dimension(10),parameter :: tagSend = (/ 21,22,23,24,25,26,27,28,29,30 /)
!      integer,dimension(10),parameter :: tagRecv = (/ 22,21,24,23,26,25,30,29,28,27 /)
!---------------------------------------------------------------------
!     Informations paralleles contenues dans une structure
!--------------------------------------------------------------------
!
      type infopar
         integer                                    ::  rank            !RANG de PROCESSUS
         integer                                    ::  subcycle=1      !Nombre de cycles dans un pas de temps principal
         integer                                    ::  nbdom           !NOMBRE DE DOMAINES
         integer                                    ::  dimax,djmax     !NOMBRE DE DOMAINES
         integer                                    ::  iglb,jglb,imax,jmax,kmaxp1       !TAILLE GLOBALE
         integer,dimension(2)                       ::  coords          !COORDONNEES DANS LA GRILLE DE PROCESSUS
         integer,dimension(2)                       ::  timax           !OFFSET DEBUT,FIN
         integer,dimension(2)                       ::  tjmax           !DANS LA NUMEROTATION GLOBALE
         integer,dimension(:,:),allocatable         ::  gtimax          !INDICES DEBUT,FIN
         integer,dimension(:,:),allocatable         ::  gtjmax          !DANS LA NUMEROTATION GLOBALE
         integer,dimension(12)                      ::  tvoisin         !LE NUMERO DES VOISINS DANS L'ORDRE(O,E,S,N)
         integer,dimension(:,:),allocatable         ::  gtvoisin        !LE NUMERO DES VOISINS DANS L'ORDRE(O,E,S,N)
         integer                      ::  comm2d          !COMMUNICATEUR GLOBAL
         integer,dimension(2,4)       ::  posechangex     !LES INDICES DES POSITIONs DES ECHANGES
         integer,dimension(2,4)       ::  posechangey     !
         integer,dimension(2,4)       ::  posechangep     !
         integer,dimension(2,4)       ::  posechanger     !
         integer,dimension(2,4)       ::  posechanger1    !  !2014-04-18
         integer,dimension(2,4)       ::  posechangez0    !
         integer,dimension(2,4)       ::  posechangez1    !
         integer,dimension(2,4)       ::  posechangez2    !
         integer,dimension(2,4)       ::  posechangez3    !
         integer,dimension(2,4)       ::  posechangex2    !
         integer,dimension(2,4)       ::  posechangey2    !
         integer,dimension(2,4)       ::  posechangev1    !
         integer,dimension(2,4)       ::  posechangev2    !
         integer,dimension(2,4)       ::  posechangeu1    !
         integer,dimension(2,4)       ::  posechangeu2    !
      end type infopar
      integer                         :: ierr
      character(len=6)                :: dom_c            !CHAINE DE CARACTER REPRODUISANT LES COORDONNEES
      type (infopar)                  :: par
      integer preq(1000), pst(1000)
      
!_________________________________________________________________

!---------------------------------------------------------------------
!     Gestion des echanges 
!--------------------------------------------------------------------
     type subarray
       integer :: nb
       integer,dimension(:),pointer :: imin,imax,jmin,jmax
     end	type
     type borneechange
       character(len=10) :: name
       character(len=10) :: type
 !      integer :: datatype
       type(subarray),dimension(12) :: Send
       type(subarray),dimension(12) :: Recv
     end type

     type all_echange
          type (borneechange) :: X,X2,Y,Y2,Z0,Z1,Z2,Z3,Z0Z1,Z1Z2,Z0Z1Z2,R,R1 !2014-04-18
     end type
    
     type (all_echange) :: all_borne_echange

     type st_echange
          integer :: nb =0
          integer,dimension(1000) :: req 
          integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
          integer,dimension(1000) :: sender
          integer,dimension(1000) :: recv
          integer,dimension(1000) :: tagsend
          integer,dimension(1000) :: tagrecv
       end type st_echange
       type (st_echange) :: l_echange ! liste des echanges

      ! Gestion des echanges, un type par variables échangées
      type sub_array_echange
         integer :: ndim
         integer :: imin
         integer :: imax
         integer :: jmin
         integer :: jmax
         integer :: kmin
         integer :: kmax
         integer :: typedata
         character(len=2) :: typechange
         integer,dimension(12) :: subtype_s, subtype_r
      end type sub_array_echange
      type(sub_array_echange),dimension(1000) :: l_subtype_echange
      integer,dimension(1000) :: indice_echange
      character(len=30),dimension(1000) :: name_echange
      integer :: nb_subtype_echange=0

      integer,dimension(:),allocatable :: oldtype, lblock
      integer,dimension(:),allocatable :: rcounts,displs

      integer, parameter   :: ki4=selected_int_kind(4)
      integer, parameter   :: ki8=selected_int_kind(8)
      integer, parameter   :: kr4=selected_real_kind(6,37)
      integer, parameter   :: kr8=selected_real_kind(15,307)

      logical :: debug=.FALSE.
      character(len=20) :: txtdebug,typinit

      character(len=100) :: filedbg
      character(len=3) :: ckount

      integer,dimension(1000) :: partabreq
      integer :: parnbreq=0,parnbappel=1
      integer,dimension(MPI_STATUS_SIZE,1000) :: partstatus
      integer :: partag1,partag2

      !temps
      double precision :: tpscom=0,tpscomtot=0
      !  Subcycling
!     double precision,dimension(:),allocatable :: glob_dte_lp, glob_dte_lp_tmp
!     ! liste des voisins au temps t
!     integer,dimension(10) :: mpi_neighbor_list ! Passe dans module_principal a cause de faire_hotrestart !14-02-14
!     integer ::  nbvoi=8

!*********************************************************************
!                       INTERFACES
!*********************************************************************

! NON BLOQUANT
!................................................................................
      interface echange
         module procedure echange4dr8_nonblock_itps                     !04-07-10
         module procedure echange4dr8_nonblock                          !27-03-14
         module procedure echange4dr4_nonblock_itps                     !01-11-13
         module procedure echange3dr8_nonblock                          !04-07-10
         module procedure echange3dr8_nonblock_itps                     !04-07-10
         module procedure echange3dr4_nonblock_itps                     !01-11-13
         module procedure echange2dr8_nonblock                          !10-10-10
         module procedure echange2di4_nonblock                          !11-10-10
         module procedure echange3di4_nonblock                          !11-10-10
      end interface

      interface get_type_echange
         module procedure get_type_echange2d_r8                         !28-10-2013
         module procedure get_type_echange2d_r4                         !22-01-2015 !22-01-15
         module procedure get_type_echange2d_r8_itps                    !28-10-2013
         module procedure get_type_echange2d_r4_itps                    !28-10-2013
         module procedure get_type_echange3d_r8                         !28-10-2013
         module procedure get_type_echange3d_r8_itps                    !28-10-2013
         module procedure get_type_echange4d_r8                         !28-10-2013
         module procedure get_type_echange4d_r4                         !28-10-2013
         module procedure get_type_echange3d_r4_itps                    !28-10-2013
         module procedure get_type_echange3d_r4                         !28-10-2013
         module procedure get_type_echange3d_i4                         !28-10-2013
         module procedure get_type_echange3d_i2         
         module procedure get_type_echange2d_i4                         !28-10-2013
         module procedure get_type_echange2d_i4_itps
         module procedure get_type_echange2d_i2_itps
      end interface
!................................................................................



      
      interface echange_voisin
         module procedure echange_voisin2dr8
         module procedure echange_voisin3dr8
         module procedure echange_voisin4dr8
         module procedure echange_voisin4dr4
         module procedure echange_voisin3dr4
         module procedure echange_voisin2dr4
         module procedure echange_voisin3di4
         module procedure echange_voisin3di2
         module procedure echange_voisin2di4
      end interface

      interface echange_simple
         module procedure echange_r8_nonblock_2
         module procedure echange_r8_nonblock_1
         module procedure echange_r8_nonblock_p
      end interface

      interface par_gatherall_2d
        module procedure par_gatherall_2d_i1
        module procedure par_gatherall_2d_i4
        module procedure par_gatherall_2d_r4
        module procedure par_gatherall_2d_r8
      end interface

      
      contains
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine initialise_parallel(dimax,djmax,iglb,jglb,          &
           imax,jmax,kmaxp1,ipbc,jpbc)
!---------------------------------------------------------------------
!      INITIALISATION DES PAREMTRES DE PARALELLISATION :
!         RANK            :-> Numero du processus courant
!         NBDOM           :-> Nombre de Domaines
!         COMM2D          :-> Communicateur intra grille de processus
!         COORDS          :-> Coordonees de process dans la grille
!         TMECO           :-> Offset pour retrouver les indices globaux
!         TNECO           :-> Offset pour retrouver les indices globaux
!         TVOISIN         :-> Numero des processus voisins O,E,N,S
!---------------------------------------------------------------------
       !IMPLICIT NONE
       integer,intent(in)    :: dimax,djmax
       integer,intent(in)    :: iglb,jglb,kmaxp1,imax,jmax
       logical,intent(in)    :: ipbc,jpbc                              !15-10-10

       integer :: nbdom
       integer,dimension(ndims) :: dims
       integer,dimension(2)     :: coords,timax,tjmax,coordsvoisin
       logical,dimension(ndims) :: periods
       integer,dimension(dimax+1,djmax+1) :: tabdebi,tabdebi2,tabdebj,tabdebj2
!       integer,dimension(djmax+1) :: tabdebj,tabdebj2
       character :: c,c1,c2
       integer :: bcl
       character(len=3) dom_c1,dom_c2

       typinit="auto"
       nbdom=dimax*djmax
       periods(1)=ipbc   ! .true. periodic along Oi axis
       periods(2)=jpbc   ! .true. periodic along Oj axis

       ! Subcycling
       !allocate(glob_dte_lp(0:nbdom-1))
       !glob_dte_lp = 0.d00
       !allocate(glob_dte_lp_tmp(0:nbdom-1))
       !glob_dte_lp_tmp = 0.d00
       par%imax=imax
       par%jmax=jmax
       par%kmaxp1=kmaxp1

!       imax=(iglb-2)/dimax+2!+DMECO-1
!       jmax=(jglb-2)/djmax+2!+DNECO-1
       par%dimax=dimax
       par%djmax=djmax

!---------------------------------------------------------------
!     Initialisation de MPI: Communicateur, Rangs,...
!      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_symp, par%nbdom, ierr)
      call mpi_comm_rank(mpi_comm_symp, par%rank, ierr)
      if (par%nbdom /= nbdom) then
         write(0,*) " probleme : "
         write(0,*) "   nombre de domaines demandes :",nbdom
         write(0,*) "   nombre de processus dispos :",par%nbdom
      endif

!--------------------------------------------------------------
!    TAILLE GLOBALE
      par%iglb=iglb
      par%jglb=jglb

!---------------------------------------------------------------
! ALLOCATION Des tableaux GLOBAUX
      allocate(par%gtimax(0:nbdom-1,2))
      allocate(par%gtjmax(0:nbdom-1,2))
!---------------------------------------------------------------
! ALLOCATION Des tableaux GLOBAUX
      allocate(par%gtvoisin(0:nbdom-1,4))

!---------------------------------------------------------------
!     Creation de la typologie cartesienne
      dims(1) = dimax
      dims(2) = djmax

      call mpi_cart_create(mpi_comm_symp,ndims,                        &
           dims,periods,.false.,                                        &
           par%comm2d,ierr)
      call mpi_cart_coords(par%comm2d,par%rank,ndims,                   &
           coords,ierr)
      par%coords=coords
!--------------------------------------------------------
!     Calcul des bornes des dimensions des tableaux
!     Pour MECO
      if (coords(1) == 0) then
         timax(1) = 0
         timax(2) = imax
      else
         timax(1) = coords(1)*(imax-2)
         timax(2) = timax(1)+imax
      endif
!     Pour NECO
      if (coords(2) == 0) then
         tjmax(1) = 0
         tjmax(2) = jmax
      else
         tjmax(1) = coords(2)*(jmax-2)
         tjmax(2) = tjmax(1)+jmax
      endif
      par%timax=timax
      par%tjmax=tjmax
!--------------------------------------------------------
!     Calcul des bornes des dimensions des tableaux
!     Pour MECO avec possibilite de imax, jmax different sur les procs
      tabdebi = 0 
      tabdebj = 0
      tabdebi2 = 0 
      tabdebj2 = 0
      tabdebi(coords(1)+2,coords(2)+2) = imax-2
      tabdebj(coords(1)+2,coords(2)+2) = jmax-2
      call mpi_allreduce(tabdebi,tabdebi2,(dimax+1)*(djmax+1),mpi_integer,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(tabdebj,tabdebj2,(dimax+1)*(djmax+1),mpi_integer,mpi_sum,par%comm2d,ierr)
      timax(1) = sum(tabdebi2(1:coords(1)+1,coords(2)+2))
      timax(2) = timax(1)+imax
      tjmax(1) = sum(tabdebj2(coords(1)+2,1:coords(2)+1))
      tjmax(2) = tjmax(1)+jmax
      par%timax=timax
      par%tjmax=tjmax
      
      par%tvoisin(:) = mpi_proc_null
!--------------------------------------------------------
!     Les Voisins
!     Voisins EST OUEST
      call mpi_cart_shift(par%comm2d,0,1,                               &
         par%tvoisin(ouest),par%tvoisin(est),ierr)
!     Voisins NORD SUD
      call mpi_cart_shift(par%comm2d,1,-1,                              &
        par%tvoisin(nord),par%tvoisin(sud),ierr)
      ! voisin nord est
      if ( (par%tvoisin(est) /= mpi_proc_null).and.(par%tvoisin(nord) /= mpi_proc_null)) then
         coordsvoisin = coords + (/ 1,1 /)
         call MPI_CART_RANK (par%comm2d,coordsvoisin,par%tvoisin(nordest),ierr)
      endif
!     Les coins                                                                    !04-07-10
      ! voisin nord ouest
      if ( (par%tvoisin(ouest) /= mpi_proc_null).and.(par%tvoisin(nord) /= mpi_proc_null)) then
         coordsvoisin = coords + (/ -1,1 /)
         call MPI_CART_RANK (par%comm2d,coordsvoisin,par%tvoisin(nordouest),ierr)
      endif
      ! voisin sud est
      if ( (par%tvoisin(est) /= mpi_proc_null).and.(par%tvoisin(sud) /= mpi_proc_null)) then
         coordsvoisin = coords + (/ 1,-1 /)
         call MPI_CART_RANK (par%comm2d,coordsvoisin,par%tvoisin(sudest),ierr)
      endif
      ! voisin sud ouest
      if ( (par%tvoisin(ouest) /= mpi_proc_null).and.(par%tvoisin(sud) /= mpi_proc_null)) then
         coordsvoisin = coords + (/ -1,-1 /)
         call MPI_CART_RANK (par%comm2d,coordsvoisin,par%tvoisin(sudouest),ierr)
      endif

!     Le numero de domaine dans une chaine de caractere
      dom_c="000000"
      dom_c1="000"
      dom_c2="000"
      if(coords(1) < 10) then
         write(dom_c1(3:3),'(i1)')coords(1)
         if(coords(1) < 100) then
            write(dom_c1(2:3),'(i2)')coords(1)
         else
            if(coords(1) < 1000) then
               write(dom_c1(1:3),'(i3)')coords(1)
            endif
         endif
      endif
      dom_c(4:6)="000"
      if(coords(2) < 10) then
         write(dom_c2(3:3),'(i1)')coords(2)
      else
         if(coords(2) < 100) then
            write(dom_c2(2:3),'(i2)')coords(2)
         else
            if(coords(2) < 1000) then
               write(dom_c2(1:3),'(i3)')coords(2)
            endif
         endif
      endif
      dom_c=dom_c1//dom_c2
!     write(0,*) par%rank,coords,dom_c
!     if(nbdom == 1) dom_c="unique"
      if(nbdom == 1) dom_c="single"     !29-10-10


!     MISE A JOUR DE VALEURS LORSQUE MGLB n'est pas un multiple
!     du nombre de domaine.
      if(par%tvoisin(est) == mpi_proc_null) then
          par%timax(2)=iglb
      endif
      if(par%tvoisin(nord) == mpi_proc_null) then
          par%tjmax(2)=jglb
      endif
!     imax=par%timax(2)-par%timax(1)
!     jmax=par%tjmax(2)-par%tjmax(1)
      timax=par%timax
      tjmax=par%tjmax

!!! Cyril --- 2014-02-11
! Pour Le SubCycling
     par%tvoisin(ouest2) = par%tvoisin(ouest)
     par%tvoisin(est2)   = par%tvoisin(est)
     
      
! LES INDICES DEBUT ET FIN DE CHAQUE DOMAINES SONT CONNUS DE TOUS
      call mpi_allgather(par%timax(1),1,mpi_integer,                    &
           par%gtimax(:,1),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%timax(2),1,mpi_integer,                    &
           par%gtimax(:,2),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%tjmax(1),1,mpi_integer,                    &
           par%gtjmax(:,1),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%tjmax(2),1,mpi_integer,                    &
           par%gtjmax(:,2),1,mpi_integer,par%comm2d,ierr)
! LES VOSINS SERONT connus de TOUs
      call mpi_allgather(par%tvoisin(1),1,mpi_integer,                  &
           par%gtvoisin(:,1),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%tvoisin(2),1,mpi_integer,                  &
           par%gtvoisin(:,2),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%tvoisin(3),1,mpi_integer,                  &
           par%gtvoisin(:,3),1,mpi_integer,par%comm2d,ierr)
      call mpi_allgather(par%tvoisin(4),1,mpi_integer,                  &
           par%gtvoisin(:,4),1,mpi_integer,par%comm2d,ierr)

!---------------------------------------------------------
! LA POSITION DES INDICES DES ECHANGES
!    Envoi P1 vers P2 puis P3 vers P4
     call set_position_echange(imax,jmax)
 
!---------------------------------------------------------
! Pour Le SubCycling
! Commente le 14-02-14
!     mpi_neighbor_list=(/ ouest, est, nord, nordouest, nordest , sud, sudouest, sudest, ouest2, est2/) !14-02-14
      call init_borne_echange(imax,jmax)
!----------------------------------------------------------
!     Initailisation des differents type pour l echange des variables
!!      call par_type_gather2d
!      allocate(oldtype(0:par%nbdom-1))
!      allocate(lblock(0:par%nbdom-1))
!      allocate(displs(0:par%nbdom-1))
      !CALL SET_RACCORD_
      !call create_mpi_type()
!---------------------------------------------------------------------
!          FIN INIT AUTOMATIQUE
!---------------------------------------------------------------------
      end subroutine initialise_parallel
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine  initialise_parallel_manu(dimax,djmax,iglb,jglb,  & !16-11-09 nouvelle routine
           imax,jmax,kmaxp1,filepara)
!---------------------------------------------------------------------
!      INITIALISATION DES PAREMTRES DE PARALELLISATION :
!         RANK            :-> Numero du processus courant
!         NBDOM           :-> Nombre de Domaines
!         COMM2D          :-> Communicateur intra grille de processus
!         COORDS          :-> Coordonees de process dans la grille
!         TMECO           :-> Offset pour retrouver les indices globaux
!         TNECO           :-> Offset pour retrouver les indices globaux
!         TVOISIN         :-> Numero des processus voisins O,E,N,S
!---------------------------------------------------------------------
       !IMPLICIT NONE
       integer,intent(in)    :: dimax,djmax
       integer,intent(in)    :: iglb,jglb,kmaxp1
       integer               :: iglbin,jglbin,i_,j_
!      integer               :: imax,jmax,iglbin,jglbin
       integer,intent(inout)   :: imax,jmax
       character(len=*),intent(in) :: filepara
!
       integer :: nbdom,bcl,dnum,bcl2,bcl3,subcycle_
       integer,dimension(ndims) :: dims
       integer,dimension(2)     :: coords,timax,tjmax
       logical,dimension(ndims) :: periods
       character :: c,c1,c2
       character(len=6) :: ctmp

       typinit="manu"
       ! Lecture du fichier de description des domaines paralleles
       open(213,file=filepara)
       read(213,*) i_,j_,nbdom !  ; nbdom=i_*j_
       if(i_/=dimax)stop 'nbdom_imax notebook_grid /= carte mpi'
       if(j_/=djmax)stop 'nbdom_jmax notebook_grid /= carte mpi'
       read(213,*) iglbin,jglbin
!---------------------------------------------------------------
! ALLOCATION Des tableaux GLOBAUX
      allocate(par%gtimax(0:nbdom-1,2))
      allocate(par%gtjmax(0:nbdom-1,2))
      ! Subcycling
      !allocate(glob_dte_lp(0:nbdom-1))
      !glob_dte_lp = 0.d00
      !allocate(glob_dte_lp_tmp(0:nbdom-1))
      !glob_dte_lp_tmp = 0.d00
!---------------------------------------------------------------
! ALLOCATION Des tableaux GLOBAUX
      allocate(par%gtvoisin(0:nbdom-1,8))
      par%gtvoisin=mpi_proc_null
      par%tvoisin=mpi_proc_null
!       imax=(iglb-2)/dimax+2!+DMECO-1
!       jmax=(jglb-2)/djmax+2!+DNECO-1
       par%dimax=dimax
       par%djmax=djmax

       par%kmaxp1=kmaxp1

!---------------------------------------------------------------
!     Initialisation de MPI: Communicateur, Rangs,...
!      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_symp, par%nbdom, ierr)
      call mpi_comm_rank(mpi_comm_symp, par%rank, ierr)
      if (par%nbdom /= nbdom) then
         write(0,*) " probleme : "  
         write(0,*) "   nombre de domaines demandes :",nbdom
         write(0,*) "   nombre de processus dispos :",par%nbdom
      endif
      if ((iglb /= iglbin) .or. (jglb /= jglbin)) then
         write(0,*) " probleme : "
         write(0,*) "   iglb,jglb sont differents"
      endif

      par%comm2d = mpi_comm_symp

       ! Lecture du fichier de description des domaines paralleles suite
       do bcl=0,nbdom-1
          read(213,*)
          read(213,*)dnum
          read(213,*)coords(1),coords(2)
          read(213,*)subcycle_
          if(bcl==par%rank) then !---->
             par%coords(1)=coords(1)
             par%coords(2)=coords(2)
             par%subcycle=subcycle_
          endif                  !---->
!         print *,par%rank,"dom_c=",dom_c
          read(213,*) par%gtimax(dnum,1), par%gtimax(dnum,2)
          read(213,*) par%gtjmax(dnum,1), par%gtjmax(dnum,2)
          read(213,*) (par%gtvoisin(dnum,bcl3),bcl3=1,8) !Ouest, Est, Nord, Sud, 
          !                        sudouest,sudest,nordouest,nordest
!         write(300+par%rank,*) par%gtvoisin(dnum,:)
!         print *,"init nbvoisin=",nbvoisin
          do bcl2=1,8 
             if (par%gtvoisin(dnum,bcl2) < 0)                  &
                 par%gtvoisin(dnum,bcl2)=mpi_proc_null
          enddo
       enddo
       par%timax = par%gtimax(par%rank,:)
       par%tjmax = par%gtjmax(par%rank,:)
       par%tvoisin(1:4) = par%gtvoisin(par%rank,1:4) ! Ouest, Est, Nord, SudEst
       par%tvoisin(7:10) = par%gtvoisin(par%rank,5:8) !sudouest,sudest,nordouest,nordest
       coords(1) = par%coords(1)
       coords(2) = par%coords(2)
       imax = par%gtimax(par%rank,2)-par%gtimax(par%rank,1)
       jmax = par%gtjmax(par%rank,2)-par%gtjmax(par%rank,1)
       par%imax = imax
       par%jmax = jmax
       close(213)

!--------------------------------------------------------------
!     Le numero de domaine dans une chaine de caractere
      dom_c="000000"
      if(coords(1) < 10) then
         write(dom_c(3:3),'(i1)')coords(1)
      else
         if(coords(1) < 100) then
            write(dom_c(2:3),'(i2)')coords(1)
         else
            if(coords(1) < 1000) then
               write(dom_c(1:3),'(i3)')coords(1)
            endif
         endif
      endif
      if(coords(2) < 10) then
         write(dom_c(6:6),'(i1)')coords(2)
      else
         if(coords(2) < 100) then
            write(dom_c(5:6),'(i2)')coords(2)
         else
            if(coords(2) < 1000) then
               write(dom_c(4:6),'(i3)')coords(2)
            endif
         endif
      endif
      write(0,*) par%rank,coords,dom_c
!     if(nbdom == 1) dom_c="unique"
      if(nbdom == 1) dom_c="single"    !29-10-10

!--------------------------------------------------------------
!    TAILLE GLOBALE
      par%iglb=iglb
      par%jglb=jglb
      
!!! Cyril --- 2014-02-11
! Pour Le SubCycling
     par%tvoisin(ouest2) = par%tvoisin(ouest)
     par%tvoisin(est2)   = par%tvoisin(est)

      
! VERIFICATION des voisins
#ifdef bidon
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (par%gtimax(par%tvoisin(ouest),2) /= (par%timax(1)+2)) then
            print *,trim(filepara)
            print *,"domaine ",par%rank," et ",par%tvoisin(ouest),      &
                 " ont un probleme de connexion A",                     &
                 par%timax(1)+2,                                        &
                 par%gtimax(par%tvoisin(ouest),2)
            stop
         endif
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (par%gtimax(par%tvoisin(est),1) /= (par%timax(2)-2)) then
            print *,trim(filepara)
            print *,"domaine ",par%rank," et ",par%tvoisin(est),        &
                 " ont un probleme de connexion B",                     &
                 par%timax(2)-2,                                        &
                 par%gtimax(par%tvoisin(est),1)
            stop
         endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (par%gtjmax(par%tvoisin(nord),1) /= (par%tjmax(2)-2)) then
            print *,trim(filepara)
            print *,"domaine ",par%rank," et ",par%tvoisin(nord),       &
                 " ont un probleme de connexion C",                     &
                 par%tjmax(2)-2,                                        &
                 par%gtjmax(par%tvoisin(nord),1)
            stop
         endif
      endif
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (par%gtjmax(par%tvoisin(sud),2) /= (par%tjmax(1)+2)) then
            print *,trim(filepara)
            print *,"domaine ",par%rank," et ",par%tvoisin(sud),        &
                 " ont un probleme de connexion D",                     &
                 par%tjmax(1)+2,                                        &
                 par%gtjmax(par%tvoisin(sud),2)
         stop
         endif
      endif
#endif


!! IMPORTANT si decentre
!      imax = par%timax(2)+2
!      jmax = par%tjmax(2)+2
!      print *,par%rank,"imax=",imax,par%timax(2),par%timax(1),
!     &  par%timax(2)-par%timax(1)
!      print *,par%rank,"jmax=",jmax,par%tjmax(2),par%tjmax(1),
!     &  par%tjmax(2)-par%tjmax(1)
!      call mpi_barrier(par%comm2d,bcl)
!      stop

!      imax = par%timax(2)-par%timax(1)
!      jmax = par%tjmax(2)-par%tjmax(1)

!---------------------------------------------------------
! LA POSITION DES INDICES DES ECHANGES
      call set_position_echange(imax,jmax)

!---------------------------------------------------------
! Pour Le SubCycling
!     mpi_neighbor_list=(/ ouest, est, nord, nordouest, nordest , sud, sudouest, sudest, ouest2, est2/)
      call init_borne_echange(imax,jmax)


!---------------------------------------------------------------------
!          FIN INIT Manuel
!---------------------------------------------------------------------

      end subroutine initialise_parallel_manu
!---------------------------------------------------------------------


!---------------------------------------------------------------------
      subroutine set_position_echange(imax,jmax)
!---------------------------------------------------------------------
!     Initialise la position des echanges
!---------------------------------------------------------
      integer,intent(in) :: imax,jmax
!
! LA POSITION DES INDICES DES ECHANGES
!    Envoi P1 vers P2 puis P3 vers P4
!     TYPE X
      par%posechangex(ouestest,:) = (/ 3,imax+1,imax-1,1 /)
      par%posechangex(nordsud,:)  = (/ 2, jmax, jmax-1,1 /)
!     TYPE X2
      par%posechangex2(ouestest,:) = (/ 4,imax+2,imax-2,0 /)         !28-05-10
      par%posechangex2(nordsud,:)  = (/ 3, jmax+1, jmax-2,0 /)       !28-05-10
!     TYPE Y
      par%posechangey(ouestest,:) = (/ 2, imax, imax-1,1 /)
      par%posechangey(nordsud,:)  = (/ 3, jmax+1, jmax-1,1 /)
!     TYPE Y2
      par%posechangey2(ouestest,:) = (/ 3, imax+1, imax-2,0 /)       !28-05-10
      par%posechangey2(nordsud,:)  = (/ 4,jmax+2,jmax-2,0 /)         !28-05-10
!     TYPE Z0
      par%posechangez0(ouestest,:) = (/ 2,imax  , imax-1, 1/)
      par%posechangez0(nordsud,:)  = (/ 2,jmax  , jmax-1, 1/)
!     TYPE Z1
      par%posechangez1(ouestest,:) = (/ 3,imax+1, imax-2, 0/)
      par%posechangez1(nordsud,:)  = (/ 3,jmax+1, jmax-2, 0/)
!     TYPE Z2 Type Z ordre 2
      par%posechangez2(ouestest,:) = (/ 4,imax+2, imax-3,-1/)
      par%posechangez2(nordsud,:)  = (/ 4,jmax+2, jmax-3,-1/)
!     TYPE Z3
      par%posechangez3(ouestest,:) = (/ 5,imax+3, imax-4,-2/)
      par%posechangez3(nordsud,:)  = (/ 5,jmax+3, jmax-4,-2/)
!     TYPE P pour la pression
      par%posechangep(ouestest,:) = (/ 2, imax, imax-1,1 /)
      par%posechangep(nordsud,:)  = (/ 2, jmax, jmax-1,1 /)
!     TYPE R
      par%posechanger(ouestest,:) = (/ 3,imax+1,imax-1,1 /)
      par%posechanger(nordsud,:)  = (/ 3,jmax+1,jmax-1,1 /)
!     TYPE R1
      par%posechanger1(ouestest,:) = (/ 4,imax+2,imax-2,0 /) !2014-04-18
      par%posechanger1(nordsud,:)  = (/ 4,jmax+2,jmax-2,0 /) !2014-04-18
!     TYPE V1
      par%posechangev1(ouestest,:) = (/ 2, imax, imax-1,1 /)
      par%posechangev1(nordsud,:)  = (/ 3, jmax+1, jmax-1,1 /)
!     TYPE V2
      par%posechangev2(ouestest,:) = (/ 3,imax+1, imax-2,0 /)
      par%posechangev2(nordsud,:)  = (/ 4,jmax+2, jmax-2,0  /)
!     TYPE U1
      par%posechangeu1(ouestest,:) = (/ 3,imax+1, imax-1, 1 /)
      par%posechangeu1(nordsud,:)  = (/ 2,jmax, jmax-1, 1/)
!     TYPE U2
      par%posechangeu2(ouestest,:) = (/ 4,imax+2, imax-2, 0/)
      par%posechangeu2(nordsud,:)  = (/ 3,jmax+1, jmax-2, 0/)

!--------------------------------------------------------

!---------------------------------------------------------------------
!          FIN SET_POSITION_ECHANGE
!---------------------------------------------------------------------
      end subroutine set_position_echange
!---------------------------------------------------------------------


!---------------------------------------------------------------------
      subroutine finalise_paralelle
      end subroutine finalise_paralelle
!---------------------------------------------------------------------
!          FIN FINALISE_PARALELLE
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!-----------------                     ------------------------------
      function get_pos(type)
      implicit none
      character(len=2),intent(in) :: type
      integer,dimension(2,4)      :: get_pos
!      WRITE(0,*) PAR%RANK,"GET_POS :",TYPE,"|"
      select case(type)
         case('x ')
            get_pos=par%posechangex
         case('x2')
            get_pos=par%posechangex2
         case('y ')
            get_pos=par%posechangey
         case('y2')
            get_pos=par%posechangey2
         case('z0')
            get_pos=par%posechangez0
         case('z1')
            get_pos=par%posechangez1
         case('p')
           get_pos=par%posechangep
         case('z2')
            get_pos=par%posechangez2
         case('z3')
            get_pos=par%posechangez3
         case('r ')
            get_pos=par%posechanger
         case('r1')
            get_pos=par%posechanger1 ! 2014-04-18
         case('v1')
           get_pos=par%posechangev1
         case('v2')
           get_pos=par%posechangev2
         case('u1')
           get_pos=par%posechangeu1
         case('u2')
           get_pos=par%posechangeu2

         case('X ')
            get_pos=par%posechangex
         case('X2')
            get_pos=par%posechangex2
         case('Y ')
            get_pos=par%posechangey
         case('Y2')
            get_pos=par%posechangey2
         case('Z0')
            get_pos=par%posechangez0
         case('Z1')
            get_pos=par%posechangez1
         case('P')
           get_pos=par%posechangep
         case('Z2')
            get_pos=par%posechangez2
         case('Z3')
            get_pos=par%posechangez3
         case('R ')
            get_pos=par%posechanger
         case('R1')
            get_pos=par%posechanger1   ! 2014-04-18
         case('V1')
           get_pos=par%posechangev1
         case('V2')
           get_pos=par%posechangev2
         case('U1')
           get_pos=par%posechangeu1
         case('U2')
           get_pos=par%posechangeu2
         case default
            write(0,*) "type=",type,"|"
            stop 'type de donnee inconnu'
       end select
      !WRITE(0,*) "GET_POS=",GET_POS
      return
      end function get_pos
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------

! -------------------------------------------------------------------------------
! Gestion des type MPI d'echange de voisin a voisin
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_r8(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(2) :: lb,ub
        character(len=*) :: name
        double precision,dimension(lb(1):ub(1),lb(2)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_DOUBLE_PRECISION)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_r8
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_r4(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(2) :: lb,ub
        character(len=*) :: name
        real,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
!write(2000+par%rank,*) "nb_subtype_echange=",nb_subtype_echange
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_REAL)
!   !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_r4
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_i4(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(2) :: lb,ub
        character(len=*) :: name
        integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange))
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_INTEGER)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_i4
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_i4_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_INTEGER,itps)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_i4_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_i2_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_BYTE,itps)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_i2_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
!--------------------------------------------------------------------------------
! !-----------------                                 ------------------------------
!       subroutine get_type_echange2d_i4(type,name,tab,lb,ub,indice)
!         implicit none
!         character(len=*),intent(in) :: type
!         integer,dimension(2) :: lb,ub
!         character(len=*) :: name
!         integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
!         integer,intent(out) :: indice
!         !
!         integer :: bcl
!         logical :: trouve=.FALSE.
!         trouve=.FALSE.
!         bcl=1
!         do while((.not. trouve).and.(bcl <= nb_subtype_echange))
!            if (name_echange(bcl) == name) then
!               indice=bcl
!               trouve=.true.
!            else
!               bcl=bcl+1
!            endif
!         end do
!         if (.not.trouve) then
!            call create_sub_array2d_i4(type,tab,lb,ub)
!            indice=nb_subtype_echange
!            name_echange(nb_subtype_echange) = name
!         endif
!       end subroutine get_type_echange2d_i4
! !--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_r8_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
!write(2000+par%rank,*) "nb_subtype_echange=",nb_subtype_echange
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
!          write(2000+par%rank,*) name," | ",name_echange(bcl)
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
!          write(2000+par%rank,*) "create_subType_gen_itps :",type,lb,ub
           call create_subType_gen_itps(type,lb,ub,MPI_DOUBLE_PRECISION,itps)
!   !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_r8_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange2d_r4_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
!write(2000+par%rank,*) "nb_subtype_echange=",nb_subtype_echange
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_REAL,itps)
!   !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange2d_r4_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_r8_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(4) :: lb,ub
        character(len=*) :: name
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
!        write(200+par%rank,*) "Pass bien là..... 1"
!    	write(200+par%rank,*) "   ","lb=",lb,"ub=",ub
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_DOUBLE_PRECISION,itps)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_r8_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange4d_r8(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(4) :: lb,ub
        character(len=*) :: name
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_DOUBLE_PRECISION)
           !Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange4d_r8
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange4d_r4(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(4) :: lb,ub
        character(len=*) :: name
        real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_REAL)
           !Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange4d_r4
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_r4_itps(type,name,tab,lb,ub,itps,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(4) :: lb,ub
        character(len=*) :: name
        real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
        integer,intent(in) :: itps
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_REAL,itps)
           !Creation du nouveau type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_r4_itps
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_r8(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_DOUBLE_PRECISION)
 	   !Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_r8
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_r4(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        real,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_REAL)
!Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_r4
!--------------------------------------------------------------------------------
! !-----------------                                 ------------------------------
!       subroutine get_type_echange3d_r8(type,name,tab,lb,ub,indice)
!         implicit none 
!         character(len=*),intent(in) :: type
!         integer,dimension(3) :: lb,ub
!         character(len=*) :: name
!         double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
!         integer,intent(out) :: indice
!         !
!         integer :: bcl
!         logical :: trouve=.FALSE.
!         trouve=.FALSE.
!         bcl=1
!         do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
!            if (name_echange(bcl) == name) then
!               indice=bcl
!               trouve=.true.
!            else
!               bcl=bcl+1
!            endif
!         end do
!         if (.not. trouve) then
!            call create_sub_array3d_r8(type,tab,lb,ub)
!            indice=nb_subtype_echange
!            name_echange(nb_subtype_echange) = name
!         endif
!       end subroutine get_type_echange3d_r8
! !--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_i4(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_INTEGER)
 	   !Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_i4
!-----------------                                 ------------------------------
      subroutine get_type_echange3d_i2(type,name,tab,lb,ub,indice)
        implicit none 
        character(len=*),intent(in) :: type
        integer,dimension(3) :: lb,ub
        character(len=*) :: name
        integer(kind=1),dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        integer,intent(out) :: indice
        !
        integer :: bcl
        logical :: trouve=.FALSE.
        trouve=.FALSE.
        bcl=1
        do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
           if (name_echange(bcl) == name) then
              indice=bcl
              trouve=.true.
           else
              bcl=bcl+1
           endif
        end do
        if (.not.trouve) then
           call create_subType_gen_itps(type,lb,ub,MPI_BYTE)
 	   !Creation du nouvea type d'echange
           indice=nb_subtype_echange
           name_echange(nb_subtype_echange) = name
        endif
      end subroutine get_type_echange3d_i2
!-----------------                                 ------------------------------
! !--------------------------------------------------------------------------------
!       subroutine get_type_echange3d_i4(type,name,tab,lb,ub,indice)
!         implicit none 
!         character(len=*),intent(in) :: type
!         integer,dimension(3) :: lb,ub
!         character(len=*) :: name
!         integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
!         integer,intent(out) :: indice
!         !
!         integer :: bcl
!         logical :: trouve=.FALSE.
!         trouve=.FALSE.
!         bcl=1
!         do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
!            if (name_echange(bcl) == name) then
!               indice=bcl
!               trouve=.true.
!            else
!               bcl=bcl+1
!            endif
!         end do
!         if (.not. trouve) then
!            call create_sub_array3d_i4(type,tab,lb,ub)
!            indice=nb_subtype_echange
!            name_echange(nb_subtype_echange) = name
!         endif
!       end subroutine get_type_echange3d_i4
! !--------------------------------------------------------------------------------
! !-----------------                                 ------------------------------
!       subroutine get_type_echange3d_r8_itps_old(type,name,tab,lb,ub,itps,indice)
!         implicit none 
!         character(len=*),intent(in) :: type
!         integer,dimension(4) :: lb,ub
!         character(len=*) :: name
!         double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
!         integer,intent(in) :: itps
!         integer,intent(out) :: indice
!         !
!         integer :: bcl
!         logical :: trouve=.FALSE.
!         trouve=.FALSE.
!         bcl=1
!         do while((.not. trouve).and.(bcl <= nb_subtype_echange)) 
!            if (name_echange(bcl) == name) then
!               indice=bcl
!               trouve=.true.
!            else
!               bcl=bcl+1
!            endif
!         end do
!         if (.not.trouve) then
!            call create_sub_array3d_r8_itps(type,tab,lb,ub,itps)
!            indice=nb_subtype_echange
!            name_echange(nb_subtype_echange) = name
!         endif
!       end subroutine get_type_echange3d_r8_itps_old
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------

!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------
     subroutine get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,send,itps)
        implicit none
        character(len=2) :: type
        integer,intent(in) :: ndim,voisin,send,itps
        integer,dimension(:) :: lb,ub,profil_sous_tab,coord_debut
        !
        integer              :: imax, jmax, ksize, ierr
        integer              :: ghostisz, ghostjsz
        integer,dimension(2,4) :: pos
        integer              :: g0,h0,pox,poy,vox1,voy1,vox2,voy2,voxx,voyy
        integer              :: noi,noj,i1,j1
        !
        imax = ub(1)-1+lb(1)
        jmax = ub(2)-1+lb(2)
        if (ndim > 2) ksize = ub(3) - lb(3) + 1
        profil_sous_tab = 1
        coord_debut = 0
        ghostisz = 2-lb(1)
        ghostjsz = 2-lb(2)  
        if (ndim == 3) then
           if (itps == -10) then
              profil_sous_tab(3) = ksize
           else
              coord_debut(3) = itps-lb(3)
           endif
        endif
        if (ndim == 4) then
           profil_sous_tab(3) = ksize
           coord_debut(4) = itps-lb(4)
        endif
        vox1=0
        voy1=0
        vox2=0
        voy2=0
        voxx=0
        voyy=0
!       
!        print *,"get_info_subarray type=",type

        select case(type)
           case('z*')
              noi = imax-2
              noj = jmax-2
              j1=ghostjsz 
              i1=ghostisz
              if(par%tvoisin(ouest) == mpi_proc_null) then
                 noi = noi+ghostisz
                 i1=0
              endif
              if(par%tvoisin(est) == mpi_proc_null) then
                 noi = noi+ghostjsz
              endif
              if(par%tvoisin(sud) == mpi_proc_null) then
                 noj = noj+ghostjsz
                 j1=0
              endif
              if(par%tvoisin(nord) == mpi_proc_null) then
                 noj = noj+ghostjsz
              endif
              if (lb(1) == -2)  pos = get_pos('z3') !pm
              if (lb(1) == -1)  pos = get_pos('z2') 
              if (lb(1) ==  0)  pos = get_pos('z1') 
              select case(voisin) 
              case(ouest)
                 g0 = pos(ouestest,4)-lb(1)
                 profil_sous_tab(1:2) = (/ ghostisz, noj /)           
                 if(send == 0) coord_debut(1:2) = (/ g0, j1  /)
                 if(send == 1) coord_debut(1:2) = (/ ghostisz, j1 /)
              case(est)
                 profil_sous_tab(1:2) = (/ ghostisz, noj /)           
                 if(send == 0) coord_debut(1:2) = (/ imax-lb(1), j1 /)
                 if(send == 1) coord_debut(1:2) = (/ imax-lb(1)-ghostisz, j1 /)
              case(nord)
                 profil_sous_tab(1:2) = (/ noi, ghostjsz /)           
                 if(send == 0) coord_debut(1:2) = (/ i1, jmax-lb(2) /)
                 if(send == 1) coord_debut(1:2) = (/ i1, jmax-lb(2)-ghostjsz /)
              case(sud)
                 g0 = pos(nordsud,4)-lb(2)
                 profil_sous_tab(1:2) = (/ noi, ghostjsz /)           
                 if(send == 0) coord_debut(1:2) = (/ i1, g0 /)
                 if(send == 1) coord_debut(1:2) =  (/ i1, ghostjsz /)
              case(nordouest)
                 g0 = pos(ouestest,4)-lb(1)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, jmax-lb(2) /)
                 if(send == 1) coord_debut(1:2) = (/ ghostisz, jmax-lb(2)-ghostjsz /)
              case(nordest)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)         
                 if(send == 0) coord_debut(1:2) = (/ imax-lb(1), jmax-lb(2) /)
                 if(send == 1) coord_debut(1:2) = &
                      (/ imax-lb(1)-ghostisz, jmax-lb(2)-ghostjsz /)
              case(sudest)
                 g0 = pos(nordsud,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)         
                 if(send == 0) coord_debut(1:2) = (/ imax-lb(1), g0 /)
                 if(send == 1) coord_debut(1:2) = (/ imax-lb(1)-ghostisz, ghostjsz /)
              case(sudouest)
                 g0 = pos(nordsud,4)-lb(2)
                 h0 = pos(ouestest,4)-lb(1)
                 profil_sous_tab(1:2) =   (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ h0, g0 /)
                 if(send == 1) coord_debut(1:2) = (/ ghostisz, ghostjsz /)
              case default
                 print *, ' Type de voisin inconu'
              end select
           case ('z0' , 'z1' , 'z2' , 'z3' )
              pos = get_pos(type)
              ghostisz = 1
              ghostjsz = 1
              pox=0
              if (type(2:2) == '1') pox=1 
              if (type(2:2) == '2') pox=2 
              if (type(2:2) == '3') pox=3 
              select case(voisin) 
              case(ouest)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-2+2*pox /)           
                 g0 = pos(ouestest,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-pox-lb(2)  /)
                 g0 = pos(ouestest,1)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-pox-lb(2) /)
              case(est)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-2+2*pox /)           
                 g0 = pos(ouestest,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-pox-lb(2) /)
                 g0 = pos(ouestest,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-pox-lb(2) /)
              case(nord)
                 profil_sous_tab(1:2) = (/ imax-2+2*pox, ghostjsz /)           
                 g0 = pos(nordsud,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-pox-lb(1), g0 /)
                 g0 = pos(nordsud,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ 2-pox-lb(1), g0 /)
              case(sud)
                 profil_sous_tab(1:2) = (/ imax-2+2*pox, ghostjsz /)           
                 g0 = pos(nordsud,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-pox-lb(1), g0 /)
                 g0 = pos(nordsud,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ 2-pox-lb(1), g0 /)
              case(nordouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(nordest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case default
                 print *, ' Type de voisin inconu'
              end select
           case('x ', 'u1', 'u2')
              pos = get_pos(type)
              ghostisz = 1
              ghostjsz = 1  
              !if (voxx > 0) voxx=voxx-1
              if (par%tvoisin(sud)    == mpi_proc_null)  voy1=2
              if (par%tvoisin(nord)   == mpi_proc_null)  voy2=0
              if (par%tvoisin(ouest)  == mpi_proc_null)  vox1=2
              if (par%tvoisin(est)    == mpi_proc_null)  vox2=0
              select case(voisin) 
              case(ouest)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-2+voy1+voy2 /)           
                 g0 = pos(ouestest,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2  /)
                 g0 = pos(ouestest,1)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2 /)
              case(est)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-2+voy1+voy2 /)           
                 g0 = pos(ouestest,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2 /)
                 g0 = pos(ouestest,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2 /)
              case(nord)
                 profil_sous_tab(1:2) = (/ imax-1+vox1+vox2, ghostjsz /)           
                 g0 = pos(nordsud,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
                 g0 = pos(nordsud,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
              case(sud)
                 profil_sous_tab(1:2) = (/ imax-1+vox1+vox2, ghostjsz /)           
                 g0 = pos(nordsud,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
                 g0 = pos(nordsud,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
              case(nordouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(nordest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case default
                 print *, ' Type de voisin inconu'
              end select
           case('y ', 'v1', 'v2')
              pos = get_pos(type)
              ghostisz = 1
              ghostjsz = 1  
              !if (voyy > 0) voyy=voyy-1
              if (par%tvoisin(sud)    == mpi_proc_null)  voy1=2
              if (par%tvoisin(nord)   == mpi_proc_null)  voy2=0
              if (par%tvoisin(ouest)  == mpi_proc_null)  vox1=2
              if (par%tvoisin(est)    == mpi_proc_null)  vox2=0
              select case(voisin) 
              case(ouest)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-1+voy1+voy2 /)           
                 g0 = pos(ouestest,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2  /)
                 g0 = pos(ouestest,1)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-lb(2)-voy1-voy2 /)
              case(est)
                 profil_sous_tab(1:2) = (/ ghostisz, jmax-1+voy1+voy2 /)           
                 g0 = pos(ouestest,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ g0, 2-lb(2) -voy1-voy2 /)
                 g0 = pos(ouestest,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ g0, 2-lb(2) -voy1-voy2 /)
              case(nord)
                 profil_sous_tab(1:2) = (/ imax-2+vox1+vox2, ghostjsz /)           
                 g0 = pos(nordsud,2)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
                 g0 = pos(nordsud,3)-lb(1)
                 if(send == 1) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
              case(sud)
                 profil_sous_tab(1:2) = (/ imax-2+vox1+vox2, ghostjsz /)           
                 g0 = pos(nordsud,4)-lb(1)
                 if(send == 0) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
                 g0 = pos(nordsud,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ 2-lb(1)-vox1-vox2, g0 /)
              case(nordouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(nordest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,2)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,3)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudest)
                 g0 = pos(ouestest,2)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,3)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case(sudouest)
                 g0 = pos(ouestest,4)-lb(1)
                 h0 = pos(nordsud ,4)-lb(2)
                 profil_sous_tab(1:2) =  (/ ghostisz, ghostjsz /)        
                 if(send == 0) coord_debut(1:2) = (/ g0, h0 /)
                 g0 = pos(ouestest,1)-lb(1)
                 h0 = pos(nordsud ,1)-lb(2)
                 if(send == 1) coord_debut(1:2) = (/ g0, h0 /)
              case default
                 print *, ' Type de voisin inconu'
              end select
              
           case default
              print *, ' Type d echange inconu'
           end select
         end subroutine get_info_subarray
!--------------------------------------------------------------------------------
!-----------------                                 ------------------------------

!---------------------------------------------------------------------
! Creation des sub_type_array d'echange
!---------------------------------------------------------------------
!-----------------                      ------------------------------
      subroutine create_subType_gen_itps(type,lb,ub,typempi,itps)
        implicit none
        character(len=*),intent(in) :: type
        integer,dimension(:),intent(in) :: lb,ub
        integer,intent(in) :: typempi
        integer,intent(in),optional :: itps
        !
        integer :: voisin
        type (borneechange) :: thisborne
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2) 
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=typempi
        l_subtype_echange(nb_subtype_echange)%typechange=type

!write(200+par%rank,*) ,par%rank,"pass bien 2 type=",type
!       write(200+par%rank,*) "lb=",lb
!       write(200+par%rank,*) "ub=",ub
	
        select case(type)
         case('z0')
	      thisborne = all_borne_echange%Z0
         case('z1')
	      thisborne = all_borne_echange%Z1
         case('z2')
	      thisborne = all_borne_echange%Z2
         case('z3')
	      thisborne = all_borne_echange%Z3
         case('za')
	      thisborne = all_borne_echange%Z0Z1
         case('zb')
	      thisborne = all_borne_echange%Z1Z2
         case('zc')
	      thisborne = all_borne_echange%Z0Z1Z2
         case('u1')
	      thisborne = all_borne_echange%X
         case('v1')
	      thisborne = all_borne_echange%Y
         case('u2')
	      thisborne = all_borne_echange%X2
         case('v2')
	      thisborne = all_borne_echange%Y2
         case('r ')
	      thisborne = all_borne_echange%R
         case('r1')
	      thisborne = all_borne_echange%R1
	!call = get_borne(type)
          case default
            write(0,*) "type=",type,"|"
            stop 'type de donnee inconnu'
        end select
        do voisin=1,12
           if ( (voisin /= haut).and.(voisin /= bas)) then
            !write(1000+par%rank,*) voisin,":",est
            !write(1000+par%rank,*) thisborne%Recv(voisin)%imin,thisborne%Recv(voisin)%imax
            !write(1000+par%rank,*) thisborne%Send(voisin)%imin,thisborne%Send(voisin)%imax
	    call create_type_itps(thisborne%Recv(voisin), &
		    lb, ub, l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), &
		    voisin, typempi, itps) 
            !write(1000+par%rank,*) "RECV:",voisin,l_subtype_echange(nb_subtype_echange)%subtype_r(voisin)
	    call create_type_itps(thisborne%Send(voisin), &
		    lb,ub, l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), &
		    voisin, typempi, itps) 
            !write(1000+par%rank,*) "SEND:",voisin,l_subtype_echange(nb_subtype_echange)%subtype_s(voisin)
	  endif
	enddo
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !       case (Z0)         
!	call create_type(all_borne_echange%Z0%Send(bcl1), &
!	      lb,ub,varR_echTypeS(bcl1), bcl1, typempi) 
!            indice=nb_subtype_echange
!            name_echange(nb_subtype_echange) = name
! l_subtype_echange(nb_subtype_echange)%subtype_r(voisin)
! l_subtype_echange(nb_subtype_echange)%subtype_s(voisin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end subroutine create_subType_gen_itps
!-----------------                      ------------------------------
      subroutine create_sub_array2d_r8(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(2) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
        !
        integer,parameter              :: ndim=2
        integer,dimension(2) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
!	  print *,"create_sub_array2d_r8 type=",type
        itps = -10
        profil_tab = shape(tab)
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo

      end subroutine create_sub_array2d_r8
!---------------------------------------------------------------------
!-----------------                      ------------------------------
      subroutine create_sub_array2d_i4(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(2) :: lb,ub
        integer,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
        !
        integer,parameter              :: ndim=2
        integer,dimension(2) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
        itps = -10
        profil_tab = shape(tab)
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo

      end subroutine create_sub_array2d_i4
!---------------------------------------------------------------------
!-----------------                      ------------------------------
      subroutine create_sub_array2d_r8_itps(type,tab,lb,ub,itps) 
        implicit none 
        character(len=2) :: type
        integer,dimension(3) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        !
        integer,parameter              :: ndim=3
        integer,dimension(3) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !

!        print *,"create_sub_array2d_r8_itps type=",type
        profil_tab = shape(tab)
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
!	print *,"type=",type
!	print *,"profil_tab=",profil_tab
!	print *,"profil_sous_tab=",profil_sous_tab
!	print *,"coord_debut=",coord_debut
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
       
      end subroutine create_sub_array2d_r8_itps
!---------------------------------------------------------------------
!-----------------                      ------------------------------
     subroutine create_sub_array3d_r8(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(3) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        !
        integer,parameter              :: ndim=3
        integer,dimension(3) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
!	  print *,"create_sub_array3d_r8 type=",type
        itps = -10
        profil_tab = shape(tab)
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
       
      end subroutine create_sub_array3d_r8
!---------------------------------------------------------------------
!-----------------                      ------------------------------
     subroutine create_sub_array3d_i4(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(3) :: lb,ub
        integer,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        !
        integer,parameter              :: ndim=3
        integer,dimension(3) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
        itps = -10
        profil_tab = shape(tab)
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_integer, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
       
      end subroutine create_sub_array3d_i4
!---------------------------------------------------------------------
!-----------------                      ------------------------------
      subroutine create_sub_array3d_r8_itps(type,tab,lb,ub,itps)
        implicit none
        character(len=2) :: type
        integer,dimension(4) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
        !
        integer,parameter              :: ndim=4
        integer,dimension(4) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: nouveau_type, code, bcl,voisin,itps
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
! 	  print *,"create_sub_array3d_r8_itps type=",type
       profil_tab = shape(tab) 
        do voisin=1,4
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
        !les coins
         do voisin=7,10
           ! Type Recv
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,0,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_r(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(voisin),ierr)           
           ! Typ Send
           call  get_info_subarray(type,ndim,lb,ub,profil_sous_tab,coord_debut,voisin,1,itps)
           call mpi_type_create_subarray(ndim, profil_tab, profil_sous_tab, coord_debut,  &
                mpi_order_fortran, mpi_double_precision, &
                l_subtype_echange(nb_subtype_echange)%subtype_s(voisin), ierr)
           call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(voisin),ierr)           
        enddo
       
      end subroutine create_sub_array3d_r8_itps
!---------------------------------------------------------------------
!-----------------                      ------------------------------

!---------------------------------------------------------------------
! Creation des sub_type_array d'echange
!---------------------------------------------------------------------
      subroutine create_subtype_array2d_r8(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(3) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2)) :: tab
        !
        integer,parameter              :: nb_dims=2
        integer,dimension(2) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: ordre, ancien_tye
        integer              :: nouveau_type, code
        integer              :: imax, jmax, ksize, ierr
        integer              :: ghostisz, ghostjsz
        integer,dimension(2,4) :: pos
        integer              :: g0,h0,pox,poy
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
        ! Attention prevenir le cas des domaines de taille différente.................
        imax = (par%iglb-2)/par%dimax+2
        jmax = (par%jglb-2)/par%djmax+2
        ! Position des échanges
        pos = get_pos(type) 
        !
 !	  print *,"create_sub_array2d_r8 type=",type
        profil_tab = shape(tab)
        ghostisz = 2-lb(1)
        ghostjsz = 2-lb(2)  
        pox=0
        poy=0

        if ((type(1:1) == "x")) then
           ghostisz = ghostisz-1
           ghostjsz = ghostjsz-1
           pox=1
        endif
        if ((type(1:1) == "y")) then
           ghostisz = ghostisz-1
           ghostjsz = ghostjsz-1
           pox=0
           poy=1
        endif

         ! Frontière OUEST-EST.........................................................
        profil_sous_tab = (/ ghostisz, jmax-2+poy /)  
        !. Ouest..................................
        !Recv
        g0 = pos(ouestest,4)-lb(1)
        coord_debut = (/ g0, ghostjsz+poy+pox /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(ouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(ouest),ierr)
        !send
        g0 = pos(ouestest,1)-lb(1)
!        print *,"pox=",pox," .............................."
        coord_debut = (/ g0, ghostjsz+poy+pox   /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(ouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(ouest),ierr)
        !. Est..................................
        !Recv
        g0 = pos(ouestest,2)-lb(1)
        coord_debut = (/ g0, ghostjsz+poy+pox   /)  
!        print *," est coord_debut r =",coord_debut
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(est), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(est),ierr)
        !print *,par%rank," 3"
        !send
        g0 = pos(ouestest,3)-lb(1)
        coord_debut = (/ g0, ghostjsz+poy+pox   /) 
 !       print *," est coord_debut s =",coord_debut
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(est), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(est),ierr)
        !.................................................................................. 
        ! Frontière NORD-SUD..........................................................
        profil_sous_tab = (/ imax-2+pox, ghostjsz /)           
 !       print *,"  imax-2+pox=", imax-2+pox
        !. Nord..................................
        !Recv
        g0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ ghostisz+pox+poy, g0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nord), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nord),ierr)
        !print *,par%rank," 5"
        !send
        g0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ ghostisz+pox+poy, g0 /) !jmax-lb(2)-ghostjsz /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nord), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nord),ierr)
        !. Sud..................................
        !Recv
        g0 = pos(nordsud,4)-lb(2)
        coord_debut = (/ ghostisz+pox+poy, g0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sud), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sud),ierr)
        !print *,par%rank," 7"
        !send
        g0 = pos(nordsud,1)-lb(2)
        coord_debut = (/ ghostisz+pox+poy, g0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sud), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sud),ierr)
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        !---------------------------       Les coins               --------------------------------
        !------------------------------------------------------------------------------------------
        profil_sous_tab = (/ ghostisz, ghostjsz /)
        ! NordOuest.............................................................................
        g0 = pos(ouestest,4)-lb(1)
        h0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ g0, h0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nordouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nordouest),ierr)
        g0 = pos(ouestest,1)-lb(1)
        h0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ g0, h0 /)   !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nordouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nordouest),ierr)
        ! NordEst.............................................................................
        g0 = pos(ouestest,2)-lb(1)
        h0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ g0, h0 /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nordest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nordest),ierr)
        g0 = pos(ouestest,3)-lb(1)
        h0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ g0, h0 /)    !.. Send ............... 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nordest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nordest),ierr)
        ! SudEst.............................................................................
        g0 = pos(ouestest,2)-lb(1)
        h0 = pos(nordsud,4)-lb(1)
        coord_debut = (/ g0, h0 /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sudest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sudest),ierr)
        g0 = pos(ouestest,3)-lb(1)
        h0 = pos(nordsud,1)-lb(1)
        coord_debut = (/ g0, h0 /)   !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sudest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sudest),ierr)
        ! SudOuest.............................................................................
        g0 = pos(ouestest,4)-lb(1)
        h0 = pos(ouestest,4)-lb(1)
        coord_debut = (/ g0, h0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sudouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sudouest),ierr)
        g0 = pos(ouestest,1)-lb(1)
        h0 = pos(ouestest,1)-lb(1)
        coord_debut = (/ g0, h0 /)  !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sudouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sudouest),ierr)
      end subroutine create_subtype_array2d_r8
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      subroutine create_subtype_array3d_r8(type,tab,lb,ub)
        implicit none
        character(len=2) :: type
        integer,dimension(3) :: lb,ub
        double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
        !
        integer,parameter              :: nb_dims=3
        integer,dimension(3) :: profil_tab, profil_sous_tab, coord_debut
        integer              :: ordre, ancien_tye
        integer              :: nouveau_type, code
        integer              :: imax, jmax, ksize, ierr
        integer              :: ghostisz, ghostjsz
        integer,dimension(2,4) :: pos
        integer              :: g0,h0,pox,poy
        character(len=2) :: ltype
        !
        nb_subtype_echange=nb_subtype_echange+1
        l_subtype_echange(nb_subtype_echange)%imin=lb(1)
        l_subtype_echange(nb_subtype_echange)%jmin=ub(1)
        l_subtype_echange(nb_subtype_echange)%imax=lb(2)
        l_subtype_echange(nb_subtype_echange)%jmax=ub(2)
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%kmax=0
        l_subtype_echange(nb_subtype_echange)%typedata=kind(tab)
        l_subtype_echange(nb_subtype_echange)%typechange=type
        !
        ! Attention prevenir le cas des domaines de taille différente.................
        imax = (par%iglb-2)/par%dimax+2
        jmax = (par%jglb-2)/par%djmax+2
        ! Position des échanges
!  	  print *,"create_sub_array3d_r8 type=",type
       pos = get_pos(type) 
        !
        ksize =ub(3) - lb(3) + 1
        profil_tab = shape(tab)
        ghostisz = 2-lb(1)
        ghostjsz = 2-lb(2)        
        pox=0
        poy=0

        if ((type(1:1) == "x")) then
           ghostisz = ghostisz-1
           ghostjsz = ghostjsz-1
           pox=1
           print *,"TYPE x..............................."
        endif
        if ((type(1:1) == "y")) then
           ghostisz = ghostisz-1
           ghostjsz = ghostjsz-1
           pox=0
           poy=1
           print *,"TYPE y...............................",shape(tab),lbound(tab)
        endif
           ghostisz = 1
           ghostjsz = 1
        
!        print *," create_subtype_array3d_r8 0.5"
        ! Frontière OUEST-EST.........................................................
        profil_sous_tab = (/ ghostisz, jmax-2+poy, ksize /)  
        !. Ouest..................................
        !Recv
        g0 = pos(ouestest,4)-lb(1)
        coord_debut = (/ g0, 2+poy+pox, 0/)  
 !       print *,par%rank,"profil_sous_tab=",profil_sous_tab, profil_tab
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(ouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(ouest),ierr)
 !       print *,par%rank," create_subtype_array3d_r8 1", g0, ghostjsz+poy+pox
        !send
        g0 = pos(ouestest,1)-lb(1)
        coord_debut = (/ g0, 2+poy+pox, 0 /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(ouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(ouest),ierr)
        print *,par%rank," create_subtype_array3d_r8 2",g0, ghostjsz+poy+pox
        !. Est..................................
        !Recv
        g0 = pos(ouestest,2)-lb(1)
        coord_debut = (/ g0, 2+poy+pox, 0   /)  
!        print *," est coord_debut r =",coord_debut
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(est), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(est),ierr)
 !       print *," create_subtype_array3d_r8 3", g0, ghostjsz+poy+pox
        !print *,par%rank," 3"
        !send
        g0 = pos(ouestest,3)-lb(1)
        coord_debut = (/ g0, 2+poy+pox, 0  /) 
!        coord_debut = (/ imax-lb(1)-ghostisz, ghostjsz, 0 /) 
!        print *," est coord_debut s =",coord_debut
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(est), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(est),ierr)
!        print *," create_subtype_array3d_r8 4"
        !.................................................................................. 
        ! Frontière NORD-SUD..........................................................
        profil_sous_tab = (/ imax-2+pox, ghostjsz, ksize /)           
        !. Nord..................................
        !Recv
!        coord_debut = (/ ghostisz, jmax-lb(2), 0 /)  
        g0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ 2+pox+poy, g0, 0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nord), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nord),ierr)
 !       print *," create_subtype_array3d_r8 5"
       !print *,par%rank," 5"
        !send
!        coord_debut = (/ ghostisz, jmax-lb(2)-ghostjsz, 0 /)  
        g0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ 2+pox+poy, g0, 0 /) !jmax-lb(2)-ghostjsz /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nord), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nord),ierr)
!        print *," create_subtype_array3d_r8 6"
        !. Sud..................................
        !Recv
 !       g0 = pos(nordsud,4)-lb(1)
 !       coord_debut = (/ ghostisz, g0, 0 /)  
        g0 = pos(nordsud,4)-lb(2)
        coord_debut = (/ 2+pox+poy, g0, 0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sud), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sud),ierr)
!        print *," create_subtype_array3d_r8 7"
        !print *,par%rank," 7"
        !send
!        coord_debut = (/ ghostisz, ghostjsz, 0 /)  
        g0 = pos(nordsud,1)-lb(2)
        coord_debut = (/ 2+pox+poy, g0, 0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sud), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sud),ierr)
!        print *," create_subtype_array3d_r8 8"
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        !---------------------------       Les coins               --------------------------------
        !------------------------------------------------------------------------------------------
        profil_sous_tab = (/ ghostisz, ghostjsz, ksize /)
        ! NordOuest.............................................................................
!        g0 = pos(ouestest,4)-lb(1)
!        coord_debut = (/ g0, jmax-lb(2), 0 /)  
        g0 = pos(ouestest,4)-lb(1)
        h0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ g0, h0, 0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nordouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nordouest),ierr)
!        coord_debut = (/ ghostisz, jmax-lb(2)-ghostjsz, 0 /)   !.. Send ...............
        g0 = pos(ouestest,1)-lb(1)
        h0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ g0, h0, 0 /)   !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nordouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nordouest),ierr)
        ! NordEst.............................................................................
        !coord_debut = (/ imax-lb(1), jmax-lb(2), 0 /) 
        g0 = pos(ouestest,2)-lb(1)
        h0 = pos(nordsud,2)-lb(2)
        coord_debut = (/ g0, h0, 0 /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(nordest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(nordest),ierr)
        !coord_debut = (/ imax-lb(1)-ghostisz, jmax-lb(2)-ghostjsz, 0 /)    !.. Send ............... 
        g0 = pos(ouestest,3)-lb(1)
        h0 = pos(nordsud,3)-lb(2)
        coord_debut = (/ g0, h0, 0 /)    !.. Send ............... 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(nordest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(nordest),ierr)
        ! SudEst.............................................................................
!        g0 = pos(nordsud,4)-lb(1)
!        coord_debut = (/ imax-lb(1), g0, 0 /) 
        g0 = pos(ouestest,2)-lb(1)
        h0 = pos(nordsud,4)-lb(1)
        coord_debut = (/ g0, h0, 0 /) 
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sudest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sudest),ierr)
!        coord_debut = (/ imax-lb(1)-ghostisz, ghostjsz, 0 /)   !.. Send ...............
        g0 = pos(ouestest,3)-lb(1)
        h0 = pos(nordsud,1)-lb(1)
        coord_debut = (/ g0, h0, 0 /)   !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sudest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sudest),ierr)
        ! SudOuest.............................................................................
!        h0 = pos(ouestest,4)-lb(1)
!        coord_debut = (/ h0, g0, 0 /)  
        g0 = pos(ouestest,4)-lb(1)
        h0 = pos(ouestest,4)-lb(1)
        coord_debut = (/ g0, h0, 0 /)  
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_r(sudouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_r(sudouest),ierr)
!        coord_debut = (/ ghostisz, ghostjsz, 0 /)  !.. Send ...............
        g0 = pos(ouestest,1)-lb(1)
        h0 = pos(ouestest,1)-lb(1)
        coord_debut = (/ g0, h0, 0 /)  !.. Send ...............
        call mpi_type_create_subarray(nb_dims, profil_tab, profil_sous_tab, coord_debut,  &
             mpi_order_fortran, mpi_double_precision, &
             l_subtype_echange(nb_subtype_echange)%subtype_s(sudouest), ierr)
        call mpi_type_commit(l_subtype_echange(nb_subtype_echange)%subtype_s(sudouest),ierr)
      end subroutine create_subtype_array3d_r8
!---------------------------------------------------------------------




!---------------------------------------------------------------------
!              ECHANGES NON BLOQUANT
!---------------------------------------------------------------------
!_________________________________________________________________
!              --------     2D I4    ----------
      subroutine echange2di4_nonblock(type,tab,lb,ub)                   !04-07-10          
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      integer                        :: coin1,coin2,coin3,coin4
      integer                        :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)


      stop 'Stop echange2di4_nonblock. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1)
         call MPI_ISSEND( frsud,  sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2)
         call MPI_ISSEND( frnord,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:)
         call MPI_ISSEND(frouest,  sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:)
         call MPI_ISSEND(frest,  sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2di4_nonblock

! Mode dissocié OI, OJ....................................................................................
      subroutine echange2di4_nonblock_oi(type,tab,lb,ub,sens)                   !27-12-11         
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      integer                        :: coin1,coin2,coin3,coin4
      integer                        :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: ouestok=.true., estok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)

      ! gestion de la direction des echanges
      if (sens > 0) then
         ouestok=.false.
         estok=.false.
         if (sens == ouest) ouestok=.true.
         if (sens == est ) estok=.true.
         if ( (sens /= ouest).and.(sens /= est) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) then
         nbsend = nbsend+1
         frouest =  tab(p1,:)
         call MPI_ISSEND(frouest,  sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
         endif
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) then
         nbsend = nbsend+1
         frest =  tab(p2,:)
         call MPI_ISSEND(frest,  sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
         endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)     
      endif
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)   
      endif
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)           
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) tab(p4,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) tab(p3,:)=frest2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
        if (ouestok)  tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
        if (estok)  tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2di4_nonblock_oi

! (OJ)
      subroutine echange2di4_nonblock_oj(type,tab,lb,ub,sens)                   !27-12-11         
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      integer                        :: coin1,coin2,coin3,coin4
      integer                        :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: nordok=.true., sudok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! gestion de la direction des echanges
      if (sens > 0) then
         nordok=.false.
         sudok=.false.
         if (sens == nord) nordok=.true.
         if (sens == sud ) sudok=.true.
         if ( (sens /= nord).and.(sens /= sud) ) then
            print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) then
         nbsend = nbsend+1
         frsud =  tab(:,p1)
         call MPI_ISSEND( frsud,  sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)      
      endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) then
         nbsend = nbsend+1
         frnord =  tab(:,p2)
         call MPI_ISSEND( frnord,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr) 
      endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr) 
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) tab(:,p4)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) tab(:,p3)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
        if (sudok) tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2di4_nonblock_oj



!_________________________________________________________________
!              --------     2D     ----------
      subroutine echange2dr8_nonblock(type,tab,lb,ub)                   !04-07-10          
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      double precision                        :: coin1,coin2,coin3,coin4
      double precision                        :: rcoin1,rcoin2,rcoin3,rcoin4
      !

      pass = pass +1
      ilu = 10000*(pass-100)


      stop 'Stop echange2dr8_nonblock. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      pos = get_pos(type)
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2dr8_nonblock

! Mode dissocié OI, OJ....................................................................................
      subroutine echange2dr8_nonblock_oi(type,tab,lb,ub,sens)                   !27-12-11
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      integer, intent(in) :: sens
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      double precision                        :: coin1,coin2,coin3,coin4
      double precision                        :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: ouestok=.true., estok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)

      ! gestion de la direction des echanges
      if (sens > 0) then
         ouestok=.false.
         estok=.false.
         if (sens == ouest) ouestok=.true.
         if (sens == est ) estok=.true.
         if ( (sens /= ouest).and.(sens /= est) ) then
       print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) then
         nbsend = nbsend+1
         frouest =  tab(p1,:)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) then
         nbsend = nbsend+1
         frest =  tab(p2,:)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
          if (ouestok) tab(p4,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
          if (estok) tab(p3,:)=frest2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
          if (ouestok) tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
          if (ouestok) tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
          if (estok) tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
          if (estok) tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2dr8_nonblock_oi

! (OJ)
      subroutine echange2dr8_nonblock_oj(type,tab,lb,ub,sens)                   !27-12-211
      character(len=2),intent(in) :: type
      integer,dimension(2),intent(in) :: lb,ub
      integer, intent(in) :: sens
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2)),                          &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      double precision                        :: coin1,coin2,coin3,coin4
      double precision                        :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: nordok=.true., sudok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie2dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)


      ! gestion de la direction des echanges
      if (sens > 0) then
         nordok=.false.
         sudok=.false.
         if (sens == nord) nordok=.true.
         if (sens == sud ) sudok=.true.
         if ( (sens /= nord).and.(sens /= sud) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) then
         nbsend = nbsend+1
         frsud =  tab(:,p1)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) then
         nbsend = nbsend+1
         frnord =  tab(:,p2)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = 1
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) then
         coin1 = tab(p1,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) then
         coin2 = tab(p1,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      sz2= 1
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) then
         coin3 = tab(p2,p4)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) then
         coin4 = tab(p2,p3)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
          if (sudok) tab(:,p4)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) tab(:,p3)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
          if (sudok) tab(p2,p4) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) tab(p2,p3) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) tab(p1,p3) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
          if (sudok) tab(p1,p4) = rcoin4
      endif 

      end subroutine echange2dr8_nonblock_oj

!_________________________________________________________________
!              --------    FIN 2D     ----------
!_________________________________________________________________

!_________________________________________________________________
!              --------     3D     ----------
      subroutine echange3dr4_nonblock_itps(type,tab,lb,ub,itps)             !04-07-10
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer,intent(in) :: itps
      real*4,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      real*4,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      real*4,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      real*4                        :: coin1,coin2,coin3,coin4
      real*4                        :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange3dr4_nonblock_itps. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echange3dr4_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=ub(1)-lb(1)+1
      sz2=ub(2)-lb(2)+1
      !pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      !print *," 1 : POS=",pos," type=",type
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,itps)
         call MPI_ISSEND( frsud,  sz1, MPI_REAL, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_REAL, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,itps)
         call MPI_ISSEND( frnord,  sz1, MPI_REAL, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_REAL, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)
      endif
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,itps)
         call MPI_ISSEND(frouest,  sz2, MPI_REAL, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_REAL, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,itps)
         call MPI_ISSEND(frest,  sz2, MPI_REAL, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_REAL, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      ! Les coins
      pos = get_pos(type)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_REAL, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_REAL, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_REAL, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_REAL, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_REAL, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_REAL, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_REAL, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_REAL, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
      endif

      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)

!!!
      pos = get_pos(type)

      !print *," 2 : POS=",pos," type=",type

      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      !print *,"P3,P4=",P3,P4,ouestest
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,itps)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,itps)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,itps)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,itps)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,itps) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,itps) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,itps) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,itps) = rcoin4
      endif

      end subroutine echange3dr4_nonblock_itps
!_________________________________________________________________
!              --------     3D     ----------
      subroutine echange3dr8_nonblock(type,tab,lb,ub)                   !04-07-10          
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      double precision,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      double precision,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange3dr8_nonblock. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,:)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,:)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,:)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3dr8_nonblock
!Echange dissocies OI, OJ
      subroutine echange3dr8_nonblock_oi(type,tab,lb,ub,sens)                   !27-12-11
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      double precision,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      double precision,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: ouestok=.true., estok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)

      ! gestion de la direction des echanges
      if (sens > 0) then
         ouestok=.false.
         estok=.false.
         if (sens == ouest) ouestok=.true.
         if (sens == est ) estok=.true.
         if ( (sens /= ouest).and.(sens /= est) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr) 
      endif
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr) 
      endif
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif 
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif

      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) tab(p4,:,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) tab(p3,:,:)=frest2
      endif
      ! Les coins.........................................................
      !...................................................................
      pos = get_pos(type)
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3dr8_nonblock_oi
! (OJ)
      subroutine echange3dr8_nonblock_oj(type,tab,lb,ub,sens)                   !27-12-11  
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      double precision,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      double precision,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: nordok=.true., sudok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)


      ! gestion de la direction des echanges
      if (sens > 0) then
         nordok=.false.
         sudok=.false.
         if (sens == nord) nordok=.true.
         if (sens == sud ) sudok=.true.
         if ( (sens /= nord).and.(sens /= sud) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      ! Les coins
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) tab(:,p4,:)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
        if (nordok)  tab(:,p3,:)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
        if (sudok) tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3dr8_nonblock_oj

!_________________________________________________________________
!_________________________________________________________________
!              --------     3D     ----------
      subroutine echange3dr8_nonblock_itps(type,tab,lb,ub,itps)             !04-07-10         
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer,intent(in) :: itps
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      double precision                        :: coin1,coin2,coin3,coin4
      double precision                        :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange3dr8_nonblock_itps. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echange3dr8_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=ub(1)-lb(1)+1
      sz2=ub(2)-lb(2)+1
      !pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      !print *," 1 : POS=",pos," type=",type
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,itps)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,itps)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,itps)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,itps)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      pos = get_pos(type)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,itps) 
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      
      !print *," 2 : POS=",pos," type=",type

      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      !print *,"P3,P4=",P3,P4,ouestest
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,itps)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,itps)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,itps)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,itps)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,itps) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,itps) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,itps) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,itps) = rcoin4
      endif 

      end subroutine echange3dr8_nonblock_itps
!_________________________________________________________________
      subroutine echange3dr8_nonblock_itps_oij(type,tab,lb,ub,itps,sens)             !27-12-11
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,intent(in) :: itps
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      
!      print *, par%rank, "echange3dr8_nonblock_itps_oi Entree"
!      print *, par%rank, "echange3dr8_nonblock_itps_oi Sortie"
      call mpi_finalize(ierr)
      stop 'STOP echange3dr8_nonblock_itps_oij'
 
      end subroutine echange3dr8_nonblock_itps_oij
!---------------------------------------------------------------------------------------------------
      subroutine echange3dr8_nonblock_itps_oi(type,tab,lb,ub,itps,sens)             !27-12-11
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,intent(in) :: itps
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
    end subroutine echange3dr8_nonblock_itps_oi
!---------------------------------------------------------------------------------------------------
      subroutine echange3dr8_nonblock_itps_oj(type,tab,lb,ub,itps,sens)             !27-12-11
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,intent(in) :: itps
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1)) :: frnord,frsud,frnord2,frsud2
      double precision                        :: coin1,coin2,coin3,coin4
      double precision                        :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: nordok=.true., sudok=.true.
      !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echange3dr8_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=ub(1)-lb(1)+1
      sz2=ub(2)-lb(2)+1
      !pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      !print *," 1 : POS=",pos," type=",type
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! gestion de la direction des echanges
      if (sens > 0) then
         nordok=.false.
         sudok=.false.
         if (sens == nord) nordok=.true.
         if (sens == sud ) sudok=.true.
         if ( (sens /= nord).and.(sens /= sud) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,itps)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,itps)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      ! Les coins
      pos = get_pos(type)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) then
         coin1 = tab(p1,p3,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) then
         coin2 = tab(p1,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) then
         coin3 = tab(p2,p4,itps)
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) then
         coin4 = tab(p2,p3,itps) 
         sz2 = 1
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      pos = get_pos(type)
      
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) tab(:,p4,itps)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) tab(:,p3,itps)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) tab(p2,p4,itps) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) tab(p2,p3,itps) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) tab(p1,p3,itps) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) tab(p1,p4,itps) = rcoin4
      endif 

      end subroutine echange3dr8_nonblock_itps_oj
!_________________________________________________________________
!_________________________________________________________________


!_________________________________________________________________
!              --------     3D     ----------
      subroutine echange3di4_nonblock(type,tab,lb,ub)                   !04-07-10          
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      integer,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      integer,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)


      stop 'Stop echange3di4_nonblock. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:)
         call MPI_ISSEND( frsud,  sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:)
         call MPI_ISSEND( frnord,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:)
         call MPI_ISSEND(frouest,  sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:)
         call MPI_ISSEND(frest,  sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,:)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,:)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,:)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3di4_nonblock
      subroutine echange3di4_nonblock_oi(type,tab,lb,ub,sens)                   !27-12-11 
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      integer,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      integer,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
      logical                        :: ouestok=.true., estok=.true.
     !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      p1=pos(ouestest,1)
      p2=pos(ouestest,3)

      ! gestion de la direction des echanges
      if (sens > 0) then
         ouestok=.false.
         estok=.false.
         if (sens == ouest) ouestok=.true.
         if (sens == est ) estok=.true.
         if ( (sens /= ouest).and.(sens /= est) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
       endif

      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         if (ouestok) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:)
         call MPI_ISSEND(frouest,  sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_INTEGER, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)   
      endif
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:)
         call MPI_ISSEND(frest,  sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_INTEGER, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)      
      endif
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr) 
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (ouestok) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)   
      endif
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)       
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
          if (ouestok) tab(p4,:,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         if (estok) tab(p3,:,:)=frest2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
          if (ouestok) tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
          if (ouestok) tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (estok) tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3di4_nonblock_oi
      subroutine echange3di4_nonblock_oj(type,tab,lb,ub,sens)                   !27-12-11  
      character(len=2),intent(in) :: type
      integer,dimension(3),intent(in) :: lb,ub
      integer, intent(in) :: sens
      integer,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),              &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      integer,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      integer,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      integer,dimension(lb(3):ub(3))            :: coin1,coin2,coin3,coin4
      integer,dimension(lb(3):ub(3))         :: rcoin1,rcoin2,rcoin3,rcoin4
       logical                        :: nordok=.true., sudok=.true.
     !
      pass = pass +1
      ilu = 10000*(pass-100)
      !print *,par%rank,"Entree echangie3dr8_nonblock_3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      pos = get_pos(type) 
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz1 = size(frsud)
      sz2 = size(frouest)

      ! gestion de la direction des echanges
      if (sens > 0) then
         nordok=.false.
         sudok=.false.
         if (sens == nord) nordok=.true.
         if (sens == sud ) sudok=.true.
         if ( (sens /= nord).and.(sens /= sud) ) then
         print *," erreur echange_oi mauvaise direction d'echange rank=",par%rank
            call mpi_finalize(ierr)
            stop 
         endif
      endif

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         if (sudok) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:)
         call MPI_ISSEND( frsud,  sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_INTEGER, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:)
         call MPI_ISSEND( frnord,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_INTEGER, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2 = size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (sudok) then
         coin1 = tab(p1,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_INTEGER, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) then
         coin2 = tab(p1,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_INTEGER, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) then
         coin3 = tab(p2,p4,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_INTEGER, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (sudok) then
         coin4 = tab(p2,p3,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_INTEGER, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
        if (sudok) tab(:,p4,:)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         if (nordok) tab(:,p3,:)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
        if (sudok) tab(p2,p4,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         if (nordok) tab(p2,p3,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         if (nordok) tab(p1,p3,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
        if (sudok) tab(p1,p4,:) = rcoin4
      endif 

      end subroutine echange3di4_nonblock_oj

!_________________________________________________________________
!              --------    FIN 3D     ----------
!_________________________________________________________________

!_________________________________________________________________
!              --------     4D     ----------
!_________________________________________________________________
!_________________________________________________________________
!              --------     4D non block     ----------
      subroutine echange4dr8_nonblock_itps(type,tab,lb,ub,itps)                !04-07-10    
      character(len=2),intent(in) :: type
      integer,dimension(4),intent(in) :: lb,ub
      integer,intent(in) :: itps
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)),  &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      double precision,dimension(lb(3):ub(3))             :: coin1,coin2,coin3,coin4
      double precision,dimension(lb(3):ub(3))             :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange4dr8_nonblock_itps. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echange4dr8_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=size(frsud)
      sz2=size(frouest)
      pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:,itps)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:,itps)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:,itps)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:,itps)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2=size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,:,itps)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,:,itps)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,:,itps)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,:,itps)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,:,itps) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,:,itps) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,:,itps) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,:,itps) = rcoin4
      endif 

      end subroutine echange4dr8_nonblock_itps
!_________________________________________________________________
!_________________________________________________________________
!              --------     4D non block     ----------
      subroutine echange4dr8_nonblock(type,tab,lb,ub)                !04-07-10    
      character(len=2),intent(in) :: type
      integer,dimension(4),intent(in) :: lb,ub
      double precision,                                                 &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)),  &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      double precision,dimension(lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: frest,frouest,frest2,frouest2
      double precision,dimension(lb(1):ub(1),lb(3):ub(3),lb(4):ub(4)) :: frnord,frsud,frnord2,frsud2
      double precision,dimension(lb(3):ub(3),lb(4):ub(4))             :: coin1,coin2,coin3,coin4
      double precision,dimension(lb(3):ub(3),lb(4):ub(4))             :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange4dr8_nonblock. Adopter nouvelle procedure!' !02-07-14

      !print *,par%rank,"Entree echange4dr8_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=size(frsud)
      sz2=size(frouest)
      pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:,:)
         call MPI_ISSEND( frsud,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:,:)
         call MPI_ISSEND( frnord,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_DOUBLE_PRECISION, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)         
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:,:)
         call MPI_ISSEND(frouest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:,:)
         call MPI_ISSEND(frest,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)         
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2=size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,:,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,:,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,:,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)                
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,:,:)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_DOUBLE_PRECISION, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)                
      endif
      
      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)
      
!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,:,:)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,:,:)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,:,:)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,:,:)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,:,:) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,:,:) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,:,:) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,:,:) = rcoin4
      endif 

      end subroutine echange4dr8_nonblock
!_________________________________________________________________
!_________________________________________________________________

      subroutine echange4dr4_nonblock_itps(type,tab,lb,ub,itps)        !01-12-13
      character(len=2),intent(in) :: type
      integer,dimension(4),intent(in) :: lb,ub
      integer,intent(in) :: itps
      real*4,                                                          &
           dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)), &
           intent(inout) :: tab
      !
      integer,dimension(2,4) :: pos
      integer :: sz1,direction,status,irc,p1,p2,p3,p4,i,j,ii,nbsend,sz2,FLAG
      logical :: postage1,postage2
      integer,parameter :: TAGOUEST=5000, TAGEST=5010, TAGSUD=6000, TAGNORD=6010
      integer,parameter :: TAGSUDOUEST=15000, TAGSUDEST=15010, TAGNORDOUEST=16000, TAGNORDEST=16010
      integer,dimension(18) :: tabreq
      integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
      integer,save :: pass=100,ilu
      !
      !
      real*4,dimension(lb(2):ub(2),lb(3):ub(3)) :: frest,frouest,frest2,frouest2
      real*4,dimension(lb(1):ub(1),lb(3):ub(3)) :: frnord,frsud,frnord2,frsud2
      real*4,dimension(lb(3):ub(3))             :: coin1,coin2,coin3,coin4
      real*4,dimension(lb(3):ub(3))             :: rcoin1,rcoin2,rcoin3,rcoin4
      !
      pass = pass +1
      ilu = 10000*(pass-100)

      stop 'Stop echange4dr4_nonblock_itps. Adopter nouvelle procedure!'!02-07-14

      !print *,par%rank,"Entree echange4dr4_nonblock_itps3"
      !RECUPERATION DES POSITIONS DES ECHANGES
      sz1=size(frsud)
      sz2=size(frouest)
      pos = get_pos(type)
      nbsend = 0
      pos = get_pos(type)
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)

      ! ECHANGE NORD-SUD
      p1=pos(nordsud,1)
      p2=pos(nordsud,3)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         nbsend = nbsend+1
         frsud =  tab(:,p1,:,itps)
         call MPI_ISSEND( frsud,  sz1, MPI_REAL, PAR%tvoisin(sud), &
              TAGSUD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frsud2, sz1, MPI_REAL, PAR%tvoisin(sud), &
              TAGNORD, par%comm2d, tabreq(nbsend), ierr)
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         nbsend = nbsend+1
         frnord =  tab(:,p2,:,itps)
         call MPI_ISSEND( frnord,  sz1, MPI_REAL, PAR%tvoisin(nord), &
              TAGNORD,   par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frnord2,  sz1, MPI_REAL, PAR%tvoisin(nord), &
              TAGSUD, par%comm2d, tabreq(nbsend), ierr)
      endif
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      ! ECHANGE OUEST-EST
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         nbsend = nbsend+1
         frouest =  tab(p1,:,:,itps)
         call MPI_ISSEND(frouest,  sz2, MPI_REAL, PAR%tvoisin(ouest), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frouest2, sz2, MPI_REAL, PAR%tvoisin(ouest), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         nbsend = nbsend+1
         frest =  tab(p2,:,:,itps)
         call MPI_ISSEND(frest,  sz2, MPI_REAL, PAR%tvoisin(est), &
              TAGEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(frest2, sz2, MPI_REAL, PAR%tvoisin(est), &
              TAGOUEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      ! Les coins
      p1=pos(ouestest,1)
      p2=pos(ouestest,3)
      p3=pos(nordsud,1)
      p4=pos(nordsud,3)
      sz2=size(coin1)
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin1 = tab(p1,p3,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin1,  sz2, MPI_REAL, PAR%tvoisin(sudouest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin1, sz2, MPI_REAL, PAR%tvoisin(sudouest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         coin2 = tab(p1,p4,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin2,  sz2, MPI_REAL, PAR%tvoisin(nordouest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin2, sz2, MPI_REAL, PAR%tvoisin(nordouest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
      endif
      sz2=size(coin3)
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin3 = tab(p2,p4,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin3,  sz2, MPI_REAL, PAR%tvoisin(nordest), &
              TAGNORDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin3, sz2, MPI_REAL, PAR%tvoisin(nordest), &
              TAGSUDOUEST, par%comm2d, tabreq(nbsend), ierr)
       endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         coin4 = tab(p2,p3,:,itps)
         nbsend = nbsend+1
         call MPI_ISSEND(coin4,  sz2, MPI_REAL, PAR%tvoisin(sudest), &
              TAGSUDEST, par%comm2d, tabreq(nbsend), ierr)
         nbsend = nbsend+1
         call MPI_IRECV(rcoin4, sz2, MPI_REAL, PAR%tvoisin(sudest), &
              TAGNORDOUEST, par%comm2d, tabreq(nbsend), ierr)
      endif

      ! La salle d'attente
      CALL MPI_WAITALL(nbsend, tabreq(1:nbsend), tstatus(:,1:nbsend), IERR)

!!!
      p3=pos(ouestest,2)
      p4=pos(ouestest,4)
      if (par%tvoisin(ouest) /= mpi_proc_null) then
         tab(p4,:,:,itps)=frouest2
      endif
      if (par%tvoisin(est) /= mpi_proc_null) then
         tab(p3,:,:,itps)=frest2
      endif
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)
      if (par%tvoisin(sud) /= mpi_proc_null) then
         tab(:,p4,:,itps)=frsud2
      endif
      if (par%tvoisin(nord) /= mpi_proc_null) then
         tab(:,p3,:,itps)=frnord2
      endif
      ! Les coins.........................................................
      !...................................................................
      p1=pos(ouestest,2)
      p2=pos(ouestest,4)
      p3=pos(nordsud,2)
      p4=pos(nordsud,4)

      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p4,:,itps) = rcoin1
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(ouest) /= mpi_proc_null)) then
         tab(p2,p3,:,itps) = rcoin2
      endif
      if ( (par%tvoisin(nord) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p3,:,itps) = rcoin3
      endif
      if ( (par%tvoisin(sud) /= mpi_proc_null).and. (par%tvoisin(est) /= mpi_proc_null)) then
         tab(p1,p4,:,itps) = rcoin4
      endif

      end subroutine echange4dr4_nonblock_itps

!_________________________________________________________________
!_________________________________________________________________


!_________________________________________________________________
!              --------    Echanges SIMPLE     ----------
!_________________________________________________________________
      subroutine echange_r8_nonblock_2(tabin,tabout,direction)
        implicit none
        double precision, dimension(:,:),intent(in) :: tabin
        double precision, dimension(:,:),intent(out) :: tabout
        integer,intent(in) :: direction
        !
        integer                               :: tabreqsend,tabreqrecv
        integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
        integer,save :: pass=100,ilu
        integer      :: sz
        integer,dimension(10),parameter :: tagsend = (/  1, 2, 3, 4, -10, -10, 7, 8, 9,10  /)
        integer,dimension(10),parameter :: tagrecv = (/  2, 1, 4, 3, -10, -10,10, 9, 8, 7  /)
        integer :: tags,tagr
        integer :: voisin
        integer :: bcl,nb
        !
!        write(500+par%rank,*) 'Envoie :',tabin
        voisin = par%tvoisin(direction)
        if (voisin /= mpi_proc_null) then
!           write(10000+par%rank,*) 'rank=',par%rank," voisin=",voisin
           tags = tagsend(direction) + 20*l_echange%nb *0
           tagr = tagrecv(direction) + 20*l_echange%nb *0
           sz = size(tabin) 
           call MPI_ISSEND( tabin,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tags,   par%comm2d, tabreqsend, ierr)
!           write(10000+par%rank,*) tabreqsend,'tagsend=',tags
           call MPI_IRECV( tabout,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tagr,   par%comm2d, tabreqrecv, ierr)
!           write(10000+par%rank,*) tabreqrecv,'tagrecv=',tagr
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqsend
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqrecv
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr 
            
           !tabout = -23000-par%rank
           nb=l_echange%nb
       do bcl=1,nb
       enddo
          
        endif
      end subroutine echange_r8_nonblock_2
!_________________________________________________________________
      subroutine echange_r8_nonblock_1(tabin,tabout,direction)
        implicit none
        double precision, dimension(:),intent(in) :: tabin
        double precision, dimension(:),intent(out) :: tabout
        integer,intent(in) :: direction
        !
        integer                               :: tabreqsend,tabreqrecv
        integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
        integer,save :: pass=100,ilu
        integer      :: sz
        integer,dimension(10),parameter :: tagsend = (/  1, 2, 3, 4, -10, -10, 7, 8, 9,10  /)
        integer,dimension(10),parameter :: tagrecv = (/  2, 1, 4, 3, -10, -10,10, 9, 8, 7  /)
        integer :: tags,tagr
        integer :: voisin
        !
        voisin = par%tvoisin(direction)
        if (voisin /= mpi_proc_null) then
           sz = size(tabin) 
           tags = tagsend(direction) + 20*l_echange%nb *0
           tagr = tagrecv(direction) + 20*l_echange%nb *0
           call MPI_ISSEND( tabin,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tags,   par%comm2d, tabreqsend, ierr)
           call MPI_IRECV( tabout,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tagr,   par%comm2d, tabreqrecv, ierr)
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqsend
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqrecv
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr          
        endif
      end subroutine echange_r8_nonblock_1
!_________________________________________________________________
      subroutine echange_r8_nonblock_p(tabin,tabout,direction)
        implicit none
        double precision, intent(in) :: tabin
        double precision, intent(out) :: tabout
        integer,intent(in) :: direction
        !
        integer                               :: tabreqsend,tabreqrecv
        integer,dimension(MPI_STATUS_SIZE,18) :: tstatus
        integer,save :: pass=100,ilu
        integer      :: sz
        integer,dimension(10),parameter :: tagsend = (/  1, 2, 3, 4, -10, -10, 7, 8, 9,10  /)
        integer,dimension(10),parameter :: tagrecv = (/  2, 1, 4, 3, -10, -10,10, 9, 8, 7  /)
        integer :: tags,tagr
        integer :: voisin
        !
        voisin = par%tvoisin(direction)
        if (voisin /= mpi_proc_null) then
           sz = 1 
           tags = tagsend(direction) + 20*l_echange%nb *0
           tagr = tagrecv(direction) + 20*l_echange%nb *0
           call MPI_ISSEND( tabin,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tagr,   par%comm2d, tabreqsend, ierr)
           call MPI_IRECV( tabout,  sz, MPI_DOUBLE_PRECISION, voisin, &
                tags,   par%comm2d, tabreqrecv, ierr)
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqsend
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr
           l_echange%nb = l_echange%nb+1 
           l_echange%req(l_echange%nb) = tabreqrecv
           l_echange%sender(l_echange%nb) = par%rank
           l_echange%recv(l_echange%nb) = voisin
           l_echange%tagsend(l_echange%nb) = tags
           l_echange%tagrecv(l_echange%nb) = tagr          
        endif
     end subroutine echange_r8_nonblock_p
!------------------------------------------------------------------      
     subroutine end_echange()
       implicit none
       integer :: nb,bcl
       !
       nb = l_echange%nb  
       print *,"       ",par%rank," voisin : ",par%tvoisin
       print *,"       ",par%rank," nb echange : ",nb
       CALL MPI_WAITALL(nb, l_echange%req(1:nb), l_echange%tstatus(:,1:nb), IERR)
!       print *,par%rank,IERR
       l_echange%nb = 0
     end subroutine end_echange
!_________________________________________________________________
!              --------    FIN Echanges SIMPLE     ----------
!_________________________________________________________________


!_________________________________________________________________
!              --------    GATHERs     ----------
!_________________________________________________________________
      subroutine par_gatherall_2d_zi4(tab_out,lbout,ubout)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
      integer(kind=ki8),dimension(0:(ubout(1)-lbout(1)),                &
                                  0:(ubout(2)-lbout(2)))                &
           ,target :: tab_out
      integer,dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2

      if (par%nbdom > 1) then                                 ! 21-10-09
      do ii=0,par%nbdom-1
         m1=par%gtimax(ii,1)
         m2=par%gtimax(ii,2)+1
         n1=par%gtjmax(ii,1)
         n2=par%gtjmax(ii,2)+1
         xp => tab_out(par%gtimax(ii,1):par%gtimax(ii,2)+1,             &
                       par%gtjmax(ii,1):par%gtjmax(ii,2)+1)
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_integer                                   & ! type d'info
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif
      end subroutine par_gatherall_2d_zi4

! . .  .  .
      subroutine par_gatherall_2d_zr8(tab_out,lbout,ubout)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
      double precision,dimension(0:(ubout(1)-lbout(1)),                 &
                                 0:(ubout(2)-lbout(2)))                 &
           ,target :: tab_out
      double precision,dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2

      if (par%nbdom > 1) then                                ! 21-10-09
      do ii=0,par%nbdom-1
         m1=par%gtimax(ii,1)
         m2=par%gtimax(ii,2)+1
         n1=par%gtjmax(ii,1)
         n2=par%gtjmax(ii,2)+1
         xp => tab_out(par%gtimax(ii,1):par%gtimax(ii,2)+1,             &
                       par%gtjmax(ii,1):par%gtjmax(ii,2)+1)
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_double_precision                          & ! type d'info
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif
      end subroutine par_gatherall_2d_zr8


! . .  .  .
      subroutine par_gatherall_2d_r8(tab_out,lbout,ubout,dims,nbdom)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
      double precision,dimension(lbout(1):ubout(1),                     &
                                 lbout(2):ubout(2))                     &
           ,target :: tab_out
      integer,dimension(4),intent(in) :: dims
      integer,intent(in) :: nbdom
      double precision,dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2
      integer,dimension(4,nbdom) :: globdim
      !
      if (par%nbdom > 1) then                               ! 21-10-09
      do ii=0,par%nbdom-1
         globdim(:,ii+1) = dims(:)                          ! 16-11-09 Debug
         call mpi_bcast(globdim(:,ii+1),4,                              &
              mpi_integer,ii,par%comm2d,ierr)
      enddo
      do ii=0, par%nbdom-1
         xp => tab_out(globdim(1,ii+1):globdim(2,ii+1),      & ! 16-11-09 Debug
              globdim(3,ii+1):globdim(4,ii+1))
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_double_precision                          & ! type d'info
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif
      end subroutine par_gatherall_2d_r8
! . .  .  .
! . .  .  .
      subroutine par_gatherall_2d_r4(tab_out,lbout,ubout,dims,nbdom)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
      real,dimension(lbout(1):ubout(1),                                 &
                                 lbout(2):ubout(2))                     &
           ,target :: tab_out
      integer,dimension(4),intent(in) :: dims
      integer,intent(in) :: nbdom
      real,dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2
      integer,dimension(4,nbdom) :: globdim
      !
      if (par%nbdom > 1) then                               ! 21-10-09
      do ii=0,par%nbdom-1
         globdim(:,ii+1) = dims(:)                          ! 16-11-09 Debug
         call mpi_bcast(globdim(:,ii+1),4,                              &
              mpi_integer,ii,par%comm2d,ierr)
      enddo
      do ii=0, par%nbdom-1
         xp => tab_out(globdim(1,ii+1):globdim(2,ii+1),      & ! 16-11-09 Debug
              globdim(3,ii+1):globdim(4,ii+1))
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_real                                      & ! type d'info  !08-11-09
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif
      end subroutine par_gatherall_2d_r4
! . .  .  .

      subroutine par_gatherall_2d_i1(tab_out,lbout,ubout,dims,nbdom)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
      integer(kind=1),dimension(lbout(1):ubout(1),                     &
                                lbout(2):ubout(2))                     &
                     ,target :: tab_out
      integer,dimension(4),intent(in) :: dims
      integer,intent(in) :: nbdom
      integer(kind=1),dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2
      integer,dimension(4,nbdom) :: globdim
      !
      if (par%nbdom > 1) then                               ! 21-10-09
      do ii=0,par%nbdom-1
         globdim(:,ii+1) = dims(:)                          ! 16-11-09 Debug
         call mpi_bcast(globdim(:,ii+1),4,                              &
              mpi_integer,ii,par%comm2d,ierr)
      enddo
      do ii=0, par%nbdom-1
         xp => tab_out(globdim(1,ii+1):globdim(2,ii+1),      & ! 16-11-09 Debug
              globdim(3,ii+1):globdim(4,ii+1))
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_integer1                                    & ! type d'info  !08-11-09
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif
      end subroutine par_gatherall_2d_i1
! . .  .  .
! . .  .  .
      subroutine par_gatherall_2d_i4(tab_out,                           &
           lbout,ubout,dims,nbdom)
      implicit none
      integer,dimension(:),intent(in) :: lbout,ubout
!      integer(kind=ki8),dimension(0:(ubout(1)-lbout(1)),
!     &                          0:(ubout(2)-lbout(2)))
!     &     ,target :: tab_out
      integer(kind=ki8),dimension(lbout(1):ubout(1),                    &
                                lbout(2):ubout(2))                      &
           ,target :: tab_out
      integer,dimension(4),intent(in) :: dims
      integer,intent(in) :: nbdom
       integer(kind=ki8),dimension(:,:),pointer :: xp
      integer :: ii,m1,m2,n1,n2
      integer,dimension(4,nbdom) :: globdim
      !
      if (par%nbdom > 1) then                                ! 21-10-09
      do ii=0,par%nbdom-1
         globdim(:,ii+1) = dims(:)                           ! 16-11-09 DeBug
         call mpi_bcast(globdim(:,ii+1),4,                              &
              mpi_integer,ii,par%comm2d,ierr)
      enddo
!      do ii=0,par%nbdom-1
!         write(100+par%rank,*) ii," => ",globdim(ii+1,:)
!      enddo
      do ii=0,par%nbdom-1
         xp => tab_out(globdim(1,ii+1):globdim(2,ii+1),      & ! 16-11-09 Debug
              globdim(3,ii+1):globdim(4,ii+1))
         call mpi_bcast(xp                                   & !tableau et debut
              ,size(xp)                                      & ! nbre d'elements
              ,mpi_integer                                   & ! type d'info
              ,ii                                            & ! expediteur. bcast envoie à tous les proc
                           ,par%comm2d,ierr)                 ! communicateur, status d'erreur
      enddo
      endif

      end subroutine par_gatherall_2d_i4

!---------------------------------------------------------------------
      subroutine init_borne_echange(imaxin,jmaxin)
    integer,intent(in) :: imaxin,jmaxin
    integer :: imax,jmax
    imax=imaxin
    jmax=jmaxin
    !!! Si pas de voisins alors extension a faire a la fin ???? 
    !!!! A verifer si les echanges sont necessaires
!------------------------------------------------------------------------------------------------
! Type Z0 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z0%Send(nord), (/ 2, imax-1/), (/ jmax-1, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%Recv(nord), (/ 2, imax-1/), (/ jmax  , jmax/) )
     !      Nord-Ouest
          call init_subarray( &
	   all_borne_echange%Z0%Send(nordouest), (/ 2, 2/) , (/ jmax-1,  jmax-1/) )
          call init_subarray( &
           all_borne_echange%Z0%Recv(nordouest), (/ 1, 1/) , (/ jmax  ,  jmax/) )
!     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z0%Send(ouest), (/ 2, 2/) , (/ 2, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(ouest), (/ 1, 1/) , (/ 2, jmax-1/) )
!     !      Ouest2 --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0%Send(ouest2), (/ 2, 2/) , (/ 1, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(ouest2), (/ 1, 1/) , (/ 1, jmax-1/) )
!     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%Z0%Send(sudouest), (/ 2, 2/) , (/ 2, 2/) )
           call init_subarray( &
 	    all_borne_echange%Z0%RECV(sudouest), (/ 1, 1/) , (/ 1, 1/) )
!     !      Sud
          call init_subarray( &
	    all_borne_echange%Z0%Send(sud), (/ 2, imax-1 /), (/  2, 2/) )
          call init_subarray( &
	    all_borne_echange%Z0%Recv(sud), (/ 2, imax-1 /), (/  1, 1/) )
    !      Sud-Est
          call init_subarray( &
	    all_borne_echange%Z0%Send(sudest), (/ imax-1, imax-1/) , (/ 2, 2/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(sudest), (/ imax, imax/) , (/ 1, 1/) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z0%Send(est), (/ imax-1, imax-1/) , (/ 2, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(est), (/ imax, imax/) , (/ 2, jmax-1/) )
      !      Est2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0%Send(est2), (/ imax-1, imax-1/) , (/ 1, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(est2), (/ imax, imax/) , (/ 1, jmax-1/) )
!     !      Nord-Est
          call init_subarray( &
	    all_borne_echange%Z0%Send(nordest), (/ imax-1, imax-1/) , (/ jmax-1, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z0%RECV(nordest), (/ imax, imax/) , (/ jmax, jmax/) )
!---------------  FIN TYPE Z0------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z1 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z1%Send(nord), (/ 2, imax-1 /), (/ jmax-2, jmax-2 /) )
          call init_subarray( &
	    all_borne_echange%Z1%Recv(nord), (/ 2, imax-1 /), (/ jmax+1, jmax+1 /) )
     !      Nord-Ouest
          call init_subarray( &
	    all_borne_echange%Z1%Send(nordouest), (/ 2, 2, 3, 3 /) , (/ jmax-2, jmax-2, jmax-2, jmax-1 /) )
          call init_subarray( &
            all_borne_echange%Z1%Recv(nordouest), (/ 0, 0, 1, 1 /) , (/ jmax, jmax+1, jmax+1, jmax+1 /) )
! !     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z1%Send(ouest), (/ 3, 3/) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(ouest), (/ 0, 0/) , (/ 2, jmax-1 /) )
! !     !      Ouest2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z1%Send(ouest2), (/ 3, 3/) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(ouest2), (/ 0, 0/) , (/ 1, jmax-1 /) )
! !     !      Sud-Ouest
          call init_subarray( &
           all_borne_echange%Z1%Send(sudouest), (/ 2, 2, 3, 3 /) , (/ 3, 3, 2, 3 /) )
          call init_subarray( &
           all_borne_echange%Z1%RECV(sudouest), (/ 0, 0, 1, 1/) , (/ 0, 1, 0, 0 /) )
!     !      Sud
          call init_subarray( &
	    all_borne_echange%Z1%Send(sud), (/ 2, imax-1/), (/  3, 3 /) )
          call init_subarray( &
	    all_borne_echange%Z1%Recv(sud), (/ 2, imax-1/), (/  0, 0 /) )
!     !      Sud-Est
          call init_subarray( &
	    all_borne_echange%Z1%Send(sudest), (/ imax-2, imax-2, imax-1, imax-1/) , (/ 2, 3, 3, 3 /) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(sudest), (/ imax, imax, imax+1, imax+1 /) , (/ 0, 0, 0, 1 /) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z1%Send(est), (/ imax-2, imax-2/) , (/ 2, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(est), (/ imax+1, imax+1/) , (/ 2, jmax-1/) )
      !      Est2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z1%Send(est2), (/ imax-2, imax-2/) , (/ 1, jmax-1/) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(est2), (/ imax+1, imax+1/) , (/ 1, jmax-1/) )
!     !      Nord-Est
          call init_subarray( &
 	    all_borne_echange%Z1%Send(nordest), (/ imax-2, imax-2, imax-1, imax-1 /) , (/ jmax-2, jmax-1, jmax-2, jmax-2 /) )
          call init_subarray( &
	    all_borne_echange%Z1%RECV(nordest), (/ imax, imax, imax+1, imax+1 /) , (/ jmax+1, jmax+1, jmax, jmax+1 /) )
! !---------------  FIN TYPE Z1------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z2 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z2%Send(nord), (/ 2, imax-1 /), (/ jmax-3, jmax-3 /) )
          call init_subarray( &
	    all_borne_echange%Z2%Recv(nord), (/ 2, imax-1 /), (/  jmax+2 , jmax+2 /) )
     !      Nord-Ouest
          call init_subarray( &
	    all_borne_echange%Z2%Send(nordouest), (/  2,  2, 3, 3, 4, 4/) , (/ jmax-3, jmax-3, jmax-3, jmax-3, jmax-3, jmax-1 /) )
          call init_subarray( &
            all_borne_echange%Z2%Recv(nordouest), (/ -1, -1, 0, 0, 1, 1 /) , (/ jmax, jmax+2, jmax+2, jmax+2, jmax+2, jmax+2  /) )
!     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z2%Send(ouest), (/  4,  4/) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(ouest), (/ -1, -1 /) , (/ 2, jmax-1 /) )
!     !      Ouest2 --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z2%Send(ouest2), (/  4,  4/) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(ouest2), (/ -1, -1 /) , (/ 1, jmax-1 /) )
!     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%Z2%Send(sudouest), (/ 2, 2, 3, 3, 4, 4 /) , (/ 4, 4, 4, 4, 2, 4 /) )
           call init_subarray( &
 	    all_borne_echange%Z2%RECV(sudouest), (/ -1, -1, 0, 0, 1, 1 /) , (/ -1, 1, -1, -1, -1, -1 /) )
!     !      Sud
          call init_subarray( &
	    all_borne_echange%Z2%Send(sud), (/ 2, imax-1 /), (/  4, 4 /) )
          call init_subarray( &
	    all_borne_echange%Z2%Recv(sud), (/ 2, imax-1 /), (/  -1, -1 /) )
!     !      Sud-Est
          call init_subarray( &
	    all_borne_echange%Z2%Send(sudest), (/ imax-3, imax-3, imax-2, imax-2, imax-1, imax-1 /) , (/ 2, 4, 4, 4, 4, 4/) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(sudest), (/ imax, imax, imax+1, imax+1, imax+2, imax+2 /) , (/ -1, -1, -1, -1, -1, 1 /) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z2%Send(est), (/ imax-3, imax-3 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(est), (/ imax+2, imax+2 /) , (/ 2, jmax-1 /) )
      !      Est2 --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z2%Send(est2), (/ imax-3, imax-3 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(est2), (/ imax+2, imax+2 /) , (/ 1, jmax-1 /) )
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%Z2%Send(nordest), (/ imax-3, imax-3, imax-2, imax-2, imax-1, imax-1 /) , (/ jmax-3, jmax-1, jmax-3, jmax-3, jmax-3, jmax-3 /) )
          call init_subarray( &
	    all_borne_echange%Z2%RECV(nordest), (/ imax, imax, imax+1, imax+1, imax+2, imax+2/), (/jmax+2, jmax+2, jmax+2, jmax+2, jmax, jmax+2 /) )
!---------------  FIN TYPE Z2------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z3 
     !      Nord Send, Recv
      call init_subarray(                                       &
       all_borne_echange%Z3%Send(nord), (/ 2     , imax-1 /)    & ! imin,imax
                                      , (/ jmax-4, jmax-4 /) )    ! jmin,jmax
      call init_subarray(                                       &
       all_borne_echange%Z3%Recv(nord), (/ 2     , imax-1 /)    & ! imin,imax
                                      , (/jmax+3 , jmax+3 /) )    ! jmin,jmax
     !      Nord-Ouest
      call init_subarray(                          &
       all_borne_echange%Z3%Send(nordouest)        &
            , (/  2    , 5     , 5     , 5     /)  & ! imin,imax 1  imin,imax 2
            , (/ jmax-4, jmax-4, jmax-3, jmax-1/) )  ! jmin,jmax 1  jmin,jmax 2
          call init_subarray(                      &
           all_borne_echange%Z3%Recv(nordouest)    &
           , (/ -2    , 1     , -2    , -2    /)   & ! imin,imax 1  imin,imax 2
           , (/ jmax+3, jmax+3, jmax  , jmax+2/) )   ! jmin,jmax 1  jmin,jmax 2

!     !      Ouest
      call init_subarray( &
       all_borne_echange%Z3%Send(ouest), (/ 5,  5    /) & ! imin,imax
                                       , (/ 2, jmax-1/) ) ! jmin,jmax
      call init_subarray( &
       all_borne_echange%Z3%RECV(ouest), (/-2, -2    /) & ! imin,imax
                                       , (/ 2, jmax-1/) ) ! jmin,jmax

!     !      Sud-Ouest
      call init_subarray(                     &
       all_borne_echange%Z3%Send(sudouest)    & 
       , (/ 2, 5, 5, 5 /)                     & ! imin,imax 1   imin,imax 2
       , (/ 5, 5, 2, 4 /) )                     ! jmin,jmax 1   jmin,jmax 2
      call init_subarray(                     &
       all_borne_echange%Z3%RECV(sudouest)    &
       , (/ -2, 1, -2, -2 /)                  & ! imin,imax 1   imin,imax 2
       , (/ -2,-2, -1,  1 /) )                  ! jmin,jmax 1   jmin,jmax 2 !17-04-15

!     !      Sud
      call init_subarray( &
       all_borne_echange%Z3%Send(sud), (/ 2, imax-1 /)  & ! imin,imax
                                     , (/ 5, 5      /) )  ! jmin,jmax
      call init_subarray( &
       all_borne_echange%Z3%Recv(sud), (/ 2, imax-1 /)  & ! imin,imax
                                     , (/-2, -2     /) )  ! jmin,jmax

!     !      Sud-Est
      call init_subarray(                      &
       all_borne_echange%Z3%Send(sudest)       &
       , (/ imax-4, imax-1, imax-4, imax-4 /)  & ! imin,imax 1    imin,imax 2
       , (/ 5     , 5     , 2     , 4      /) )  ! jmin,jmax 1    jmin,jmax 2
      call init_subarray(                      &
       all_borne_echange%Z3%RECV(sudest)       & 
       , (/ imax, imax+3, imax+3, imax+3 /)    & ! imin,imax 1    imin,imax 2
       , (/ -2  , -2    , -1    , 1      /) )    ! jmin,jmax 1    jmin,jmax 2

      !      Est
      call init_subarray( &
       all_borne_echange%Z3%Send(est), (/ imax-4, imax-4 /) & ! imin,imax
                                     , (/ 2     , jmax-1 /) ) ! jmin,jmax
      call init_subarray( &
       all_borne_echange%Z3%RECV(est), (/ imax+3, imax+3 /) & ! imin,imax
                                     , (/ 2     , jmax-1 /) ) ! jmin,jmax

!     !      Nord-Est
       call init_subarray(                          & 
        all_borne_echange%Z3%Send(nordest)          &
            , (/ imax-4, imax-1, imax-4, imax-4 /)  & ! imin,imax 1   imin,imax 2
            , (/ jmax-4, jmax-4, jmax-3, jmax-1 /) )  ! jmin,jmax 1   jmin,jmax 2
          call init_subarray(                       &
           all_borne_echange%Z3%RECV(nordest)       & 
           , (/ imax  , imax+3,  imax+3, imax+3/)   & ! imin,imax 1   imin,imax 2
           , (/ jmax+3, jmax+3,  jmax  , jmax+2/) )   ! jmin,jmax 1   jmin,jmax 2

!---------------  FIN TYPE Z3------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z0Z1Z2 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(nord), (/ 2, imax-1 /), (/ jmax-3, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Recv(nord), (/ 2, imax-1 /), (/ jmax, jmax+2  /) )
     !      Nord-Ouest
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(nordouest), (/ 2, 4 /) , (/ jmax-3, jmax-1  /) )
          call init_subarray( &
           all_borne_echange%Z0Z1Z2%Recv(nordouest), (/ -1, 1 /) , (/ jmax, jmax+2  /) )
!     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(ouest), (/  2, 4 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(ouest), (/ -1, 1 /) , (/ 2, jmax-1 /) )
!     !      Ouest2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(ouest2), (/  2, 4 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(ouest2), (/ -1, 1 /) , (/ 1, jmax-1 /) )
!     !      Sud-Ouest
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(sudouest), (/  2, 4 /) , (/ 2, 4 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(sudouest), (/ -1, 1 /) , (/ -1, 1 /) )
!     !      Sud
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(sud), (/ 2, imax-1 /), (/  2, 4 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Recv(sud), (/ 2, imax-1 /), (/  -1, 1 /) )
!     !      Sud-Est 
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(sudest), (/ imax-3 ,  imax-1 /) , (/ 2, 4 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(sudest), (/  imax   , imax+2 /) , (/ -1, 1 /) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(est), (/ imax-3 ,  imax-1 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(est), (/ imax   ,  imax+2 /) , (/ 2, jmax-1 /) )
      !      Est2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(est2), (/ imax-3 ,  imax-1 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(est2), (/ imax   ,  imax+2 /) , (/ 1, jmax-1 /) )
!     !      Nord-Est
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%Send(nordest), (/ imax-3 ,  imax-1 /) , (/ jmax-3, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1Z2%RECV(nordest), (/ imax   ,  imax+2 /) , (/ jmax, jmax+2 /) )
!---------------  FIN TYPE Z0Z1Z2------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z1Z2 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(nord), (/ 2, imax-1 /), (/ jmax-3, jmax-2  /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%Recv(nord), (/ 2, imax-1 /), (/ jmax+1, jmax+2  /) )
     !      Nord-Ouest
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(nordouest), (/ 2, 4, 3, 4 /) , (/ jmax-3, jmax-2, jmax-1, jmax-1 /) )
          call init_subarray( &
            all_borne_echange%Z1Z2%Recv(nordouest), (/ -1, 1, -1, 0 /) , (/ jmax+1,  jmax+2, jmax, jmax /) )
!     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(ouest), (/  3, 4 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(ouest), (/ -1, 0 /) , (/ 2, jmax-1 /) )
!     !      Ouest2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(ouest2), (/  3, 4 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(ouest2), (/ -1, 0 /) , (/ 1, jmax-1 /) )
!     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%Z1Z2%Send(sudouest), (/ 3, 4,  2, 4 /) , (/ 2, 2, 3, 4/) )
           call init_subarray( &
 	    all_borne_echange%Z1Z2%RECV(sudouest), (/ -1, 1, -1, 0 /) , (/ -1, 0, 1, 1/) )
!     !      Sud
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(sud), (/ 2, imax-1 /), (/  3, 4 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%Recv(sud), (/ 2, imax-1 /), (/  -1, 0 /) )
!     !      Sud-Est
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(sudest), (/ imax-3, imax-1, imax-3, imax-2 /) , (/ 3, 4, 2, 2 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(sudest), (/ imax, imax+2, imax+1, imax+2 /) , (/ -1, 0,  1, 1/) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(est), (/ imax-3, imax-2 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(est), (/ imax+1, imax+2 /) , (/ 2, jmax-1 /) )
      !      Est2  --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(est2), (/ imax-3, imax-2 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(est2), (/ imax+1, imax+2 /) , (/ 1, jmax-1 /) )
!     !      Nord-Est
          call init_subarray( &
	    all_borne_echange%Z1Z2%Send(nordest), (/ imax-3, imax-1, imax-3, imax-2 /) , (/ jmax-3, jmax-2, jmax-1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z1Z2%RECV(nordest), (/ imax+1, imax+2, imax, imax+2/) , (/ jmax, jmax, jmax+1, jmax+2 /) )
!---------------  FIN TYPE Z1Z2------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Z0Z1 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(nord), (/ 2, imax-1 /), (/ jmax-2, jmax-1  /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%Recv(nord), (/ 2, imax-1 /), (/ jmax, jmax+1  /)  )
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%Z0Z1%Send(nordouest), (/ 2, 3 /) , (/ jmax-2, jmax-1 /) )
           call init_subarray( &
            all_borne_echange%Z0Z1%Recv(nordouest), (/ 0, 1 /) , (/ jmax,  jmax+1/) )
!     !      Ouest
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(ouest), (/ 2, 3 /) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%RECV(ouest), (/ 0, 1 /) , (/ 2, jmax-1 /) )
!     !      Ouest2 --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(ouest2), (/ 2, 3 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%RECV(ouest2), (/ 0, 1 /) , (/ 1, jmax-1 /) )
     !      Sud-Ouest
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(sudouest), (/ 2, 3 /) , (/ 2, 3 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%RECV(sudouest), (/ 0, 1 /) , (/ 0, 1 /) )
    !      Sud
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(sud), (/ 2, imax-1 /), (/  2, 3 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%Recv(sud), (/ 2, imax-1 /), (/  0, 1 /) )
!     !      Sud-Est
           call init_subarray( &
 	    all_borne_echange%Z0Z1%Send(sudest), (/ imax-2, imax-1 /) , (/ 2, 3/) )
           call init_subarray( &
 	    all_borne_echange%Z0Z1%RECV(sudest), (/ imax, imax+1 /) , (/ 0, 1/) )
      !      Est
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(est), (/ imax-2, imax-1/) , (/ 2, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%RECV(est), (/ imax, imax+1/) , (/ 2, jmax-1 /) )
      !      Est2 --- Cyril 2014-02-11
          call init_subarray( &
	    all_borne_echange%Z0Z1%Send(est2), (/ imax-2, imax-1/) , (/ 1, jmax-1 /) )
          call init_subarray( &
	    all_borne_echange%Z0Z1%RECV(est2), (/ imax, imax+1/) , (/ 1, jmax-1 /) )
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%Z0Z1%Send(nordest), (/ imax-2, imax-1 /) , (/ jmax-2, jmax-1 /) )
           call init_subarray( &
 	    all_borne_echange%Z0Z1%RECV(nordest), (/ imax, imax+1/) , (/ jmax, jmax+1 /) )
!---------------  FIN TYPE Z0Z1------------------------------------------------------------------------!------------------------------------------------------------------------------------------------
! Type X ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
           call init_subarray( &
	    all_borne_echange%X%Send(nord), (/ 2, imax /), (/ jmax-1, jmax-1  /) )
           call init_subarray( &
	    all_borne_echange%X%Recv(nord), (/ 2, imax /), (/ jmax  , jmax    /) ) !....
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%X%Send(nordouest), (/ 3, 3 /) , (/ jmax-1, jmax-1  /) ) !.....
           call init_subarray( &
	    all_borne_echange%X%Recv(nordouest), (/ 1, 1 /) , (/ jmax  , jmax    /) ) !.....
!     !      Ouest
           call init_subarray( &
 	    all_borne_echange%X%Send(ouest), (/ 3, 3 /) , (/ 2, jmax-1/) ) !....
           call init_subarray( &
 	    all_borne_echange%X%RECV(ouest), (/ 1, 1 /) , (/ 2, jmax-1/) ) !....
!     !      Ouest2 --- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%X%Send(ouest2), (/ 3, 3 /) , (/ 1, jmax-1/) ) !....
           call init_subarray( &
 	    all_borne_echange%X%RECV(ouest2), (/ 1, 1 /) , (/ 1, jmax-1/) ) !....
!     !      Sud-Ouest
          call init_subarray( &
	    all_borne_echange%X%Send(sudouest), (/ 3, 3 /) , (/ 2, 2/) ) !....
          call init_subarray( &
	    all_borne_echange%X%RECV(sudouest), (/ 1, 1 /) , (/ 1, 1/) ) !....
!     !      Sud
           call init_subarray( &
	    all_borne_echange%X%Send(sud), (/ 2, imax /), (/  2,  2  /) ) !....
           call init_subarray( &
	    all_borne_echange%X%Recv(sud), (/ 2, imax /), (/  1,  1  /) ) !....
!     !      Sud-Est
          call init_subarray( &
	    all_borne_echange%X%Send(sudest), (/ imax-1, imax-1 /) , (/ 2, 2 /) ) !....
          call init_subarray( &
	    all_borne_echange%X%RECV(sudest), (/ imax+1, imax+1 /) , (/ 1, 1 /) ) !....
      !      Est
           call init_subarray( &
 	    all_borne_echange%X%Send(est), (/ imax-1, imax-1 /) , (/ 2, jmax-1 /) ) !....
           call init_subarray( &
 	    all_borne_echange%X%RECV(est), (/ imax+1, imax+1 /) , (/ 2, jmax-1 /) ) !....
      !      Est2 --- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%X%Send(est2), (/ imax-1, imax-1 /) , (/ 1, jmax-1 /) ) !....
           call init_subarray( &
 	    all_borne_echange%X%RECV(est2), (/ imax+1, imax+1 /) , (/ 1, jmax-1 /) ) !....
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%X%Send(nordest), (/ imax-1, imax-1 /) , (/ jmax-1, jmax-1 /) ) !....
           call init_subarray( &
 	    all_borne_echange%X%RECV(nordest), (/ imax+1, imax+1 /) , (/ jmax, jmax /) ) !....
!---------------  FIN TYPE X-------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type X2 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
           call init_subarray( &
	    all_borne_echange%X2%Send(nord), (/ 2, imax /), (/ jmax-2, jmax-2 /) )
           call init_subarray( &
	    all_borne_echange%X2%Recv(nord), (/ 2, imax /), (/ jmax+1, jmax+1 /) )
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%X2%Send(nordouest), (/ 3, 3, 4, 4 /) , (/ jmax-2, jmax-2, jmax-2, jmax-1 /) )
           call init_subarray( &
            all_borne_echange%X2%Recv(nordouest), (/ 0, 0, 1, 1/) , (/ jmax, jmax+1, jmax+1, jmax+1 /) )
!     !      Ouest
           call init_subarray( &
 	    all_borne_echange%X2%Send(ouest), (/ 4, 4 /) , (/ 2, jmax-1 /) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(ouest), (/ 0, 0 /) , (/ 2, jmax-1/) )
!!     !      Ouest2 --- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%X2%Send(ouest2), (/ 4, 4 /) , (/ 1, jmax-1 /) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(ouest2), (/ 0, 0 /) , (/ 1, jmax-1/) )
!     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%X2%Send(sudouest), (/ 3, 3, 4, 4/) , (/ 3, 3, 2, 3 /) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(sudouest), (/ 0, 0, 1, 1/) , (/ 0, 1, 0, 0 /) )
!     !      Sud
           call init_subarray( &
	    all_borne_echange%X2%Send(sud), (/ 2, imax /), (/  3, 3/) )
           call init_subarray( &
	    all_borne_echange%X2%Recv(sud), (/ 2, imax /), (/  0, 0 /) )
!     !      Sud-Est
           call init_subarray( &
 	    all_borne_echange%X2%Send(sudest), (/ imax-2, imax-2, imax-1, imax-1 /) , (/ 2, 3, 3, 3/) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(sudest), (/ imax+1, imax+1, imax+2, imax+2 /) , (/ 0, 0, 0, 1/) )
      !      Est
           call init_subarray( &
 	    all_borne_echange%X2%Send(est), (/ imax-2, imax-2 /) , (/ 2, jmax-1 /) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(est), (/ imax+2, imax+2 /) , (/ 2, jmax-1 /) )
       !      Est2 --- Cyril 2014-02-11
          call init_subarray( &
 	    all_borne_echange%X2%Send(est2), (/ imax-2, imax-2 /) , (/ 1, jmax-1 /) )
          call init_subarray( &
 	    all_borne_echange%X2%RECV(est2), (/ imax+2, imax+2 /) , (/ 1, jmax-1 /) )
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%X2%Send(nordest), (/ imax-2, imax-2, imax-1, imax-1 /) , (/ jmax-2, jmax-1, jmax-2, jmax-2 /) )
           call init_subarray( &
 	    all_borne_echange%X2%RECV(nordest), (/ imax+1, imax+1, imax+2, imax+2 /) , (/ jmax+1, jmax+1, jmax, jmax+1 /) )
!---------------  FIN TYPE X2------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Y ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
          call init_subarray( &
	    all_borne_echange%Y%Send(nord), (/ 2, imax-1 /), (/ jmax-1, jmax-1/) )
	  call init_subarray( &
	    all_borne_echange%Y%Recv(nord), (/ 2, imax-1 /), (/ jmax+1, jmax+1 /) )
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%Y%Send(nordouest), (/ 2, 2/) , (/ jmax-1,  jmax-1/) )
           call init_subarray( &
            all_borne_echange%Y%Recv(nordouest), (/ 1, 1/) , (/ jmax+1,  jmax+1 /) )
!     !      Ouest
           call init_subarray( &
 	    all_borne_echange%Y%Send(ouest), (/ 2, 2/) , (/ 2, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(ouest), (/ 1, 1 /) , (/ 2, jmax /) )
!     !      Ouest2 --- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%Y%Send(ouest2), (/ 2, 2 /) , (/ 1, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(ouest2), (/ 1, 1 /) , (/ 1, jmax /) )
!     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%Y%Send(sudouest), (/ 2, 2/) , (/ 3, 3 /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(sudouest), (/ 1, 1/) , (/ 1, 1 /) )
!     !      Sud
           call init_subarray( &
	    all_borne_echange%Y%Send(sud), (/ 2, imax-1 /), (/  3, 3 /) )
           call init_subarray( &
	    all_borne_echange%Y%Recv(sud), (/ 2, imax-1 /), (/  1, 1 /) )
!     !      Sud-Est
           call init_subarray( &
 	    all_borne_echange%Y%Send(sudest), (/ imax-1, imax-1 /) , (/ 3 , 3 /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(sudest), (/ imax, imax/) , (/ 1, 1/) )
      !      Est
           call init_subarray( &
 	    all_borne_echange%Y%Send(est), (/ imax-1, imax-1 /) , (/ 2, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(est), (/ imax, imax /) , (/ 2, jmax /) )
      !      Est2 --- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%Y%Send(est2), (/ imax-1, imax-1 /) , (/ 1, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(est2), (/ imax, imax /) , (/ 1, jmax /) )
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%Y%Send(nordest), (/ imax-1, imax-1 /) , (/ jmax-1, jmax-1/) )
           call init_subarray( &
 	    all_borne_echange%Y%RECV(nordest), (/ imax, imax/) , (/ jmax+1, jmax+1 /) )
!---------------  FIN TYPE Y------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Type Y2 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
           call init_subarray( &
	    all_borne_echange%Y2%Send(nord), (/ 2, imax-1 /), (/ jmax-2, jmax-2 /) )
           call init_subarray( &
	    all_borne_echange%Y2%Recv(nord), (/  2, imax-1 /), (/ jmax+2, jmax+2 /) )
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%Y2%Send(nordouest), (/ 2, 2, 3, 3 /) , (/ jmax-2, jmax-2, jmax-2, jmax-1 /) )
           call init_subarray( &
            all_borne_echange%Y2%Recv(nordouest), (/ 0, 0, 1, 1 /) , (/ jmax+1,  jmax+2, jmax+2, jmax+2 /) )
!     !      Ouest
           call init_subarray( & 
 	    all_borne_echange%Y2%Send(ouest), (/ 3, 3/) , (/ 2, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(ouest), (/ 0, 0 /) , (/ 2, jmax /) )
!     !      Ouest2--- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%Y2%Send(ouest2), (/ 3, 3/) , (/ 1, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(ouest2), (/ 0, 0 /) , (/ 1, jmax /) )
!     !      Sud-Ouest
          call init_subarray( &
	    all_borne_echange%Y2%Send(sudouest), (/ 2, 2, 3, 3 /) , (/ 4, 4, 3, 4 /) )
          call init_subarray( &
	    all_borne_echange%Y2%RECV(sudouest), (/ 0, 0, 1, 1/) , (/ 0, 1, 0, 0 /) )
!     !      Sud
           call init_subarray( &
	    all_borne_echange%Y2%Send(sud), (/ 2, imax-1 /), (/  4, 4 /) )
           call init_subarray( &
	    all_borne_echange%Y2%Recv(sud), (/ 2, imax-1 /), (/  0, 0 /) )
!     !      Sud-Est
           call init_subarray( &
 	    all_borne_echange%Y2%Send(sudest), (/ imax-2, imax-2, imax-1, imax-1 /) , (/ 3, 4, 4, 4 /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(sudest), (/ imax, imax, imax+1, imax+1 /) , (/ 0, 0, 0, 1 /) )
!      !      Est
           call init_subarray( &
 	    all_borne_echange%Y2%Send(est), (/ imax-2, imax-2 /) , (/ 2, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(est), (/ imax+1, imax+1 /) , (/ 2, jmax /) )
!      !      Est2--- Cyril 2014-02-11
           call init_subarray( &
 	    all_borne_echange%Y2%Send(est2), (/ imax-2, imax-2 /) , (/ 1, jmax /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(est2), (/ imax+1, imax+1 /) , (/ 1, jmax /) )
!     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%Y2%Send(nordest), (/ imax-2, imax-2, imax-1, imax-1 /) , (/ jmax-2, jmax-1, jmax-2, jmax-2 /) )
           call init_subarray( &
 	    all_borne_echange%Y2%RECV(nordest), (/ imax, imax, imax+1, imax+1 /) , (/ jmax+2, jmax+2, jmax+1, jmax+2 /) )
!---------------  FIN TYPE Y2------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
! Type R ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
! Internal loop  2:imax & 2:jmax
! v(2:imax,jmax+1)=v_N (2:imax,3     )
! v(1     ,jmax+1)=v_NW(imax-1,3     )
! v(1     ,2:jmax)=v_W (imax-1,2:jmax)
! v(1     ,1     )=v_SW(imax-1,jmax-1)
! v(2:imax,1     )=v_S (2:imax,jmax-1)
! v(imax+1,1     )=v_SE(3     ,jmax-1)
! v(imax+1,2:jmax)=v_E (3     ,2:jmax)
! v(imax+1,jmax+1)=v_NE(3     ,3     )

     !      Nord Send, Recv
           call init_subarray( &
      all_borne_echange%R%Send(nord), (/ 2, imax /), (/ jmax-1, jmax-1 /) ) !07-02-15
           call init_subarray( &
      all_borne_echange%R%Recv(nord), (/ 2, imax /), (/ jmax+1, jmax+1 /) ) !06-02-15! Bug fix Cyril
     !      Nord-Ouest
           call init_subarray( &
      all_borne_echange%R%Send(nordouest), (/ 3, 3 /) , (/ jmax-1,  jmax-1 /) )
           call init_subarray( &
      all_borne_echange%R%Recv(nordouest), (/ 1, 1 /) , (/ jmax+1,  jmax+1 /) )
!     !      Ouest
           call init_subarray( &
      all_borne_echange%R%Send(ouest), (/ 3, 3 /) , (/ 2, jmax /) )
           call init_subarray( &
      all_borne_echange%R%RECV(ouest), (/ 1, 1 /) , (/ 2, jmax /) )
!     !      Ouest2--- Cyril 2014-02-11
           call init_subarray( &
      all_borne_echange%R%Send(ouest2), (/ 3, 3/) , (/ 1, jmax-1 /) )
           call init_subarray( &
      all_borne_echange%R%RECV(ouest2), (/ 1, 1 /) , (/ 1, jmax-1 /) )
!     !      Sud-Ouest
           call init_subarray( &
      all_borne_echange%R%Send(sudouest), (/ 3, 3/) , (/ 3, 3 /) )
           call init_subarray( &
      all_borne_echange%R%RECV(sudouest), (/ 1, 1 /) , (/ 1, 1 /) )
!     !      Sud
           call init_subarray( &
      all_borne_echange%R%Send(sud), (/ 2, imax /), (/  3, 3 /) )
           call init_subarray( &
      all_borne_echange%R%Recv(sud), (/ 2, imax /), (/  1, 1 /) )
!     !      Sud-Est
           call init_subarray( &
      all_borne_echange%R%Send(sudest), (/ imax-1, imax-1 /) , (/ 3, 3 /) )
           call init_subarray( &
      all_borne_echange%R%RECV(sudest), (/ imax+1, imax+1 /) , (/ 1, 1 /) )
      !      Est
           call init_subarray( &
      all_borne_echange%R%Send(est), (/ imax-1, imax-1 /) , (/ 2, jmax /) )
           call init_subarray( &
      all_borne_echange%R%RECV(est), (/ imax+1, imax+1 /) , (/ 2, jmax /) )
      !      Est2--- Cyril 2014-02-11
           call init_subarray( &
      all_borne_echange%R%Send(est2), (/ imax-1, imax-1 /) , (/ 1, jmax-1 /) )
           call init_subarray( &
      all_borne_echange%R%RECV(est2), (/ imax+1, imax+1 /) , (/ 1, jmax-1 /) )
!     !      Nord-Est
           call init_subarray( &
      all_borne_echange%R%Send(nordest), (/ imax-1, imax-1 /) , (/ jmax-1, jmax-1 /) )
           call init_subarray( &
      all_borne_echange%R%RECV(nordest), (/ imax+1, imax+1/) , (/ jmax+1, jmax+1 /) )
!---------------  FIN TYPE R------------------------------------------------------------------------
!!!!
!!!!   2014-04-18     EchangeR1 .....
!------------------------------------------------------------------------------------------------
! Type R1 ! (/ imin, imax /)*N , (/ jmin, jmax /)*N ............................................
     !      Nord Send, Recv
           call init_subarray( &
	    all_borne_echange%R1%Send(nord), (/ 1, imax+1 /), (/ jmax-2, jmax-2 /) ) !...
           call init_subarray( &
	    all_borne_echange%R1%Recv(nord), (/ 1, imax+1 /), (/ jmax+2, jmax+2 /) ) !...
     !      Nord-Ouest
           call init_subarray( &
	    all_borne_echange%R1%Send(nordouest), (/ 3, 3, 4, 4  /) , (/ jmax-2,  jmax-2, jmax-2, jmax-1 /) ) !...
           call init_subarray( &
            all_borne_echange%R1%Recv(nordouest), (/ 0, 0, 1, 1  /) , (/jmax+1, jmax+2, jmax+2, jmax+2 /) )  !...
!     !      Ouest
           call init_subarray( &
 	    all_borne_echange%R1%Send(ouest), (/ 4, 4/) , (/ 1, jmax+1 /) ) !...
           call init_subarray( &
 	    all_borne_echange%R1%RECV(ouest), (/ 0, 0 /) , (/ 1, jmax+1 /) )  !...
!  !     !      Ouest2--- Cyril 2014-019 ---- Attention non mis a jour pour R1 -----------
!             call init_subarray( &
!   	    all_borne_echange%R1%Send(ouest2), (/ 3, 3/) , (/ 1, jmax-1 /) )
!             call init_subarray( &
!   	    all_borne_echange%R1%RECV(ouest2), (/ 1, 1 /) , (/ 1, jmax-1 /) )
! !     !      Sud-Ouest
           call init_subarray( &
 	    all_borne_echange%R1%Send(sudouest), (/ 3, 3, 4, 4/) , (/ 4, 4, 3, 4 /) )  !...
           call init_subarray( &
 	    all_borne_echange%R1%RECV(sudouest), (/ 0, 0, 1, 1  /) , (/ 0, 1, 0, 0  /) )  !...
!     !      Sud
           call init_subarray( &
	    all_borne_echange%R1%Send(sud), (/ 1, imax+1 /), (/  4, 4 /) )
           call init_subarray( &
	    all_borne_echange%R1%Recv(sud), (/ 1, imax+1 /), (/  0, 0 /) )  !...
!     !      Sud-Est
           call init_subarray( &
 	    all_borne_echange%R1%Send(sudest), (/ imax-2, imax-2, imax-1, imax-1  /) , (/ 3, 4, 4, 4  /) ) !...
           call init_subarray( &
 	    all_borne_echange%R1%RECV(sudest), (/ imax+1, imax+1, imax+2, imax+2  /) , (/ 0, 0, 0, 1 /) )  !...
      !      Est
           call init_subarray( &
 	    all_borne_echange%R1%Send(est), (/ imax-2, imax-2 /) , (/ 1, jmax+1 /) ) !...
           call init_subarray( &
 	    all_borne_echange%R1%RECV(est), (/ imax+2, imax+2 /) , (/ 1, jmax+1 /) )  !...
!       !      Est2--- Cyril 2014-02-11 --- Attention non mis a jour pour R1 -----------
!            call init_subarray( &
!  	    all_borne_echange%R1%Send(est2), (/ imax-1, imax-1 /) , (/ 1, jmax-1 /) );;;;;;;
!            call init_subarray( &
!  	    all_borne_echange%R1%RECV(est2), (/ imax+1, imax+1 /) , (/ 1, jmax-1 /) );;;;;;;
     !      Nord-Est
           call init_subarray( &
 	    all_borne_echange%R1%Send(nordest), (/ imax-2, imax-2, imax-1, imax-1  /) , (/ jmax-2, jmax-1, jmax-2, jmax-2  /) ) !...
           call init_subarray( &
 	    all_borne_echange%R1%RECV(nordest), (/ imax+1, imax+1, imax+2, imax+2 /) , (/ jmax+2, jmax+2, jmax+1, jmax+2 /) )  !...
!---------------  FIN TYPE R1------------------------------------------------------------------------

    end subroutine init_borne_echange

   
!---------------------------------------------------------------------
      subroutine init_subarray(ST,arrI,arrJ)
    implicit none
    type(subarray),intent(inout) :: ST
    integer,dimension(:) :: arrI
    integer,dimension(:) :: arrJ
    integer :: nbblk,bcl,ii
    nbblk = size(arrI)/2
    allocate(ST%imin(nbblk))
    allocate(ST%imax(nbblk))
    allocate(ST%jmin(nbblk))
    allocate(ST%jmax(nbblk))
    ST%nb=nbblk
    ii=1
    do bcl=1,nbblk
        ST%imin(bcl)=arrI(ii)
	ST%imax(bcl)=arrI(ii+1)
	ST%jmin(bcl)=arrJ(ii)
	ST%jmax(bcl)=arrJ(ii+1)
	ii=ii+2
    enddo	
    end	subroutine init_subarray

    
!---------------------------------------------------------------------
   subroutine create_type_itps(subarr,lb,ub,type, voisin, typevarmpi, itps)
      implicit none 
      type(subarray),intent(in)  :: subarr
      integer,dimension(:),intent(in) :: lb,ub
      integer,intent(out) :: type
      integer,intent(in) :: typevarmpi, voisin
      integer,intent(in),optional :: itps
      integer,dimension(size(lb)) :: start,nbeltfull,nbeltsub
      integer,dimension(10) :: ltype,blocklens, offsets, oldtype
      integer :: bcl,ierr,ndims,rank
      type=0
      ndims = size(lb)
      start = 0
      nbeltfull = ub(:)-lb(:)+1
      nbeltsub  = ub(:)-lb(:)+1
      if ( present(itps) ) then
	nbeltsub(ndims) = 1	
	start(ndims) = itps-lb(ndims)	
      endif
      do bcl=1,subarr%nb
	nbeltsub(1) = subarr%imax(bcl)-subarr%imin(bcl)+1
 	nbeltsub(2) = subarr%jmax(bcl)-subarr%jmin(bcl)+1
 	start(1) = subarr%imin(bcl)-lb(1)
 	start(2) = subarr%jmin(bcl)-lb(1)
!  	! Gestion des domaine du bord...............................
! 	if ((par%tvoisin(ouest) == mpi_proc_null).and. &
! 	  ((voisin == sud).or.(voisin == nord)) ) then
! 	    nbeltsub(1) = nbeltsub(1)+start(1)
! 	    start(1) = 0
! 	endif
! 	if ((par%tvoisin(sud) == mpi_proc_null).and.&
! 	  ((voisin == ouest).or.(voisin == est)) ) then
! 	    nbeltsub(2) = nbeltsub(2)+start(2)
! 	    start(2) = 0
! 	endif
! 	
! 	if ((par%tvoisin(est) == mpi_proc_null).and. &
! 	  ((voisin == sud).or.(voisin == nord)) ) then
! 	    nbeltsub(1) = (ub(1)-lb(1))-start(1)+1
! 	endif
! 	if ((par%tvoisin(nord) == mpi_proc_null).and.&
! 	  ((voisin == ouest).or.(voisin == est)) ) then
! 	    nbeltsub(2) = (ub(2)-lb(2))-start(2)+1
! 	endif
!! 24-06-2014   Bug sur les tours du domaine 
!! quand on elève un proc null il n'est plus forcement sans coin
 	! Gestion des domaine du bord...............................
	if ((par%tvoisin(ouest) == mpi_proc_null)) then
	   if ((voisin == sud).and.(par%tvoisin(sudouest) == mpi_proc_null )) then
		nbeltsub(1) = nbeltsub(1)+start(1)
 	        start(1) = 0
 	   endif
	   if ((voisin == nord).and.(par%tvoisin(nordouest) == mpi_proc_null )) then
		nbeltsub(1) = nbeltsub(1)+start(1)
 	        start(1) = 0
 	   endif
 	endif
	if ((par%tvoisin(est) == mpi_proc_null)) then
	   if ((voisin == sud).and.(par%tvoisin(sudest) == mpi_proc_null )) nbeltsub(1) = (ub(1)-lb(1))-start(1)+1
	   if ((voisin == nord).and.(par%tvoisin(nordest) == mpi_proc_null )) nbeltsub(1) = (ub(1)-lb(1))-start(1)+1
	endif
	if ((par%tvoisin(nord) == mpi_proc_null)) then
	   if ((voisin == ouest).and.(par%tvoisin(nordouest) == mpi_proc_null )) nbeltsub(2) = (ub(2)-lb(2))-start(2)+1
	   if ((voisin == est).and.(par%tvoisin(nordest) == mpi_proc_null )) nbeltsub(2) = (ub(2)-lb(2))-start(2)+1
	endif
	if ((par%tvoisin(sud) == mpi_proc_null)) then
	   if ((voisin == ouest).and.(par%tvoisin(sudouest) == mpi_proc_null )) then
		nbeltsub(2) = nbeltsub(2)+start(2)
 	        start(2) = 0
 	   endif
	   if ((voisin == est).and.(par%tvoisin(sudest) == mpi_proc_null )) then
		nbeltsub(2) = nbeltsub(2)+start(2)
 	        start(2) = 0
 	   endif
 	 endif
!! 24-06-2014   FIN.......

!	!write(1000+par%rank,*) "nbeltsub=",nbeltsub
!  	! Gestion des domaine du bord FIN...............................
	call MPI_Type_create_subarray( ndims, &  ! number of array dimensions (positive integer
	  nbeltfull,                       &  ! number of elements of type oldtype in each dimension of the full array (array of positive integers)
	  nbeltsub,                        &  ! number of elements of type oldtype in each dimension of the subarray (array of positive integers)
	  start,                           &  ! starting coordinates of the subarray in each dimension (array of nonnegative integers)
	  MPI_ORDER_FORTRAN,               &  ! array storage order flag (state)
	  typevarmpi,                      &  ! array element datatype (handle)
	  ltype(bcl),                      &
	  ierr)
	call mpi_type_commit(ltype(bcl),ierr)
     enddo
     type =  ltype(1)
     !! Regroupement des lignes
     blocklens = 1
     offsets = 0
     oldtype = ltype
     call MPI_Type_struct(subarr%nb,blocklens,offsets,oldtype,type, ierr)
     call mpi_type_commit(type,ierr)
   end subroutine create_type_itps
 !---------------------------------------------------------------------
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                Echanges Dissocies         !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------
! Finalisation des échanges
!---------------------------------------------------------------------
  subroutine loc_wait()
    integer :: bcl,bcl2
    !
!    do bcl =1,parnbreq
!       write(200+par%rank,*) bcl, partabreq(bcl)
!    enddo

     CALL MPI_WAITALL(parnbreq, partabreq(1:parnbreq), partstatus(:,1:parnbreq), IERR)
!    do bcl =1,parnbreq
!       call mpi_wait(partabreq(bcl),partstatus(:,bcl), ierr)
!    enddo
!    call mpi_barrier(par%comm2d,ierr)
!     write(2000+par%rank,*) "MPI_WAITALL Ok, ierr=",ierr
    parnbreq=0
    parnbappel=0
!---------------------------------------------------------------------
!          FIN LOC_WAIT
!---------------------------------------------------------------------
  end subroutine loc_wait

!---------------------------------------------------------------------
! Envoie/Reception entre 2 procs
!---------------------------------------------------------------------
      subroutine echange_voisin4dr8(tab,type,ivoisin)
      implicit none  
      double precision,dimension(:,:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
!    write(200+par%rank,*) "echange_voisin4dr8 type=",type,"ivoisin=",ivoisin
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*1000
       tag2 = tagRecv(ivoisin)+type*1000
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
      end subroutine echange_voisin4dr8
!---------------------------------------------------------------------
!          FIN sendredcv4dr8
!---------------------------------------------------------------------
      subroutine echange_voisin4dr4(tab,type,ivoisin)
      implicit none  
      real,dimension(:,:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
      end subroutine echange_voisin4dr4
!---------------------------------------------------------------------
      subroutine echange_voisin3dr8(tab,type,ivoisin)
      implicit none  
      double precision,dimension(:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin) +type*100
       tag2 = tagRecv(ivoisin) +type*100
!       !write(1000+par%rank,*) "::typeR=",typeR
!       !write(1000+par%rank,*) "::typeS=",typeS
       partag1=tag1
       partag2=tag2
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
      end subroutine echange_voisin3dr8
!---------------------------------------------------------------------
!          FIN echange_voisin3dr8
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine echange_voisin3dr4(tab,type,ivoisin)
      implicit none  
      real,dimension(:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin) +type*100
       tag2 = tagRecv(ivoisin) +type*100
       partag1=tag1
       partag2=tag2
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
      end subroutine echange_voisin3dr4
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine echange_voisin2dr4(tab,type,ivoisin)
      implicit none
      real,dimension(:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), &
            tag1, par%comm2d, partabreq(parnbreq), ierr)
      endif
      end subroutine echange_voisin2dr4
!---------------------------------------------------------------------
!          FIN echange_voisin2dr4
!---------------------------------------------------------------------
      subroutine echange_voisin3di4(tab,type,ivoisin)
      implicit none  
      integer,dimension(:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
     !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
      end subroutine echange_voisin3di4
!---------------------------------------------------------------------
!          FIN echange_voisin3di4
!---------------------------------------------------------------------
      subroutine echange_voisin3di2(tab,type,ivoisin)
      implicit none
      integer(kind=1),dimension(:,:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
     !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), &
            tag1, par%comm2d, partabreq(parnbreq), ierr)
      endif
      end subroutine echange_voisin3di2
!---------------------------------------------------------------------
!          FIN echange_voisin3di2
!---------------------------------------------------------------------
      subroutine echange_voisin2dr8(tab,type,ivoisin)
      implicit none
      double precision,dimension(:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), &
            tag1, par%comm2d, partabreq(parnbreq), ierr)
      endif
      end subroutine echange_voisin2dr8
!---------------------------------------------------------------------
!          FIN echange_voisin2dr8
!---------------------------------------------------------------------
      subroutine echange_voisin2di4(tab,type,ivoisin)
      implicit none
      integer,dimension(:,:) :: tab
      integer,intent(in) :: type, ivoisin
      integer            :: typeS,typeR, tag1, tag2
      !
      if ( par%tvoisin(ivoisin) /= mpi_proc_null) then
       typeS =  l_subtype_echange(type)%subtype_s(ivoisin)
       typeR =  l_subtype_echange(type)%subtype_r(ivoisin)
       tag1 = tagSend(ivoisin)+type*100
       tag2 = tagRecv(ivoisin)+type*100
       parnbreq = parnbreq+1
       call MPI_IRECV( tab,  &
            1, typeR, par%tvoisin(ivoisin), &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS,  par%tvoisin(ivoisin), &
            tag1, par%comm2d, partabreq(parnbreq), ierr)
      endif
      end subroutine echange_voisin2di4
!---------------------------------------------------------------------
!          FIN echange_voisin2di4
!---------------------------------------------------------------------


!---------------------------------------------------------------------
! Envoie/Reception entre 2 procs
!---------------------------------------------------------------------
      subroutine sendredcv4dr8(tab,PosS,PosR,type,voisin,lb,ub,tag1,tag2)
      implicit none  
      integer,dimension(4),intent(in) :: lb,ub
      double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
      integer,dimension(4) :: PosS,PosR
      integer,intent(in) :: type, voisin, tag1, tag2
      !
      if (voisin /= mpi_proc_null) then
       parnbreq = parnbreq+1
       call MPI_IRECV( tab(PosR(1),PosR(2),PosR(3),PosR(4)),  &
            1, type, voisin, &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab(PosS(1),PosS(2),PosS(3),PosS(4)),  &
            1, type, voisin, & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
      endif
!---------------------------------------------------------------------
!          FIN sendredcv4dr8
!---------------------------------------------------------------------
      end subroutine sendredcv4dr8
!---------------------------------------------------------------------
! Envoie/Reception entre 2 procs
!---------------------------------------------------------------------
  subroutine sendredcv3dr8(tab,PosS,PosR,type,voisin,lb,ub,tag1,tag2)
    implicit none  
    integer,dimension(3),intent(in) :: lb,ub
    double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
    integer,dimension(3) :: PosS,PosR
    integer,intent(in) :: type, voisin, tag1, tag2
    !
    if (voisin /= mpi_proc_null) then
       parnbreq = parnbreq+1
       call MPI_IRECV( tab(PosR(1),PosR(2),PosR(3)),  &
            1, type, voisin, &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab(PosS(1),PosS(2),PosS(3)),  &
            1, type, voisin, & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
    endif
!---------------------------------------------------------------------
!          FIN sendredcv3dr8
!---------------------------------------------------------------------
  end subroutine sendredcv3dr8

!---------------------------------------------------------------------
! Envoie/Reception entre 2 procs
!---------------------------------------------------------------------
  subroutine sendredcv4dr8_sub(tab,typeS,typeR,voisin,lb,ub,tag1,tag2)
    implicit none  
    integer,dimension(4),intent(in) :: lb,ub
    double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: tab
    integer,dimension(4) :: PosS,PosR
    integer,intent(in) :: typeS,typeR, voisin, tag1, tag2
    !
    if (voisin /= mpi_proc_null) then
       parnbreq = parnbreq+1 
       call MPI_IRECV( tab,  &
            1, typeR, voisin, &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS, voisin, & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
    endif
  end subroutine sendredcv4dr8_sub
!---------------------------------------------------------------------
!          FIN sendredcv4dr8
!---------------------------------------------------------------------
  subroutine sendredcv3dr8_sub(tab,typeS,typeR,voisin,lb,ub,tag1,tag2)
    implicit none  
    integer,dimension(3),intent(in) :: lb,ub
    double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
    integer,dimension(3) :: PosS,PosR
    integer,intent(in) :: typeS,typeR, voisin, tag1, tag2
    !
    if (voisin /= mpi_proc_null) then 
       parnbreq = parnbreq+1
       call MPI_IRECV( tab,  &
            1, typeR, voisin, &
            tag2, par%comm2d, partabreq(parnbreq), ierr)
       parnbreq = parnbreq+1
       call MPI_ISSEND(tab,  &
            1, typeS, voisin, & 
            tag1, par%comm2d, partabreq(parnbreq), ierr)  
    endif
  end subroutine sendredcv3dr8_sub
!---------------------------------------------------------------------
!          FIN sendredcv3dr8
!---------------------------------------------------------------------
!!$ subroutine sendredcv2dr8_sub_itps(tab,typeS,typeR,voisin,lb,ub,tag1,tag2,itps)
!!$    implicit none  
!!$    integer,dimension(3),intent(in) :: lb,ub
!!$    double precision,dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) :: tab
!!$    integer,dimension(3) :: PosS,PosR
!!$    integer,intent(in) :: typeS,typeR, voisin, tag1, tag2, itps
!!$    !
!!$    if (voisin /= mpi_proc_null) then
!!$       parnbreq = parnbreq+1
!!$       call MPI_IRECV( tab(lb(1),lb(2),itps),  &
!!$            1, typeR, voisin, &
!!$            tag2, par%comm2d, partabreq(parnbreq), ierr)
!!$       parnbreq = parnbreq+1
!!$       call MPI_ISSEND(tab(lb(1),lb(2),itps),  &
!!$            1, typeS, voisin, & 
!!$            tag1, par%comm2d, partabreq(parnbreq), ierr)  
!!$    endif
!!$  end subroutine sendredcv2dr8_sub_itps
! ---------------------------------------------------------------------
!          FIN sendredcv2dr8
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Creation des types MPI d'echange
!---------------------------------------------------------------------
!!$  subroutine create_mpi_type()
!!$    implicit none
!!$    integer :: type_colonne, type_ligne
!!$    integer :: imax,jmax, kmax
!!$    ! 1 : Les echange de plan
!!$      ! type double precison
!!$    imax=par%imax
!!$    jmax=par%jmax
!!$    kmax=par%kmaxp1+1
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$! pour les tableau   (0:imax+1,0:jmax+1,0:kmax+1) 
!!$!+++++ Colonne tab(ii,:,kk,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_CONTIGUOUS(imax-2,mpi_double_precision,type_colonne_1_0,ierr)
!!$!    call MPI_TYPE_CONTIGUOUS(imax+2,mpi_double_precision,type_colonne_1_0,ierr)
!!$    !! Colonne tous les j=1,jmax, tableau avec indice i constant
!!$    call MPI_TYPE_COMMIT(type_colonne_1_0,ierr)
!!$    call MPI_TYPE_CONTIGUOUS(imax,mpi_double_precision,type_colonne_1_1,ierr)
!!$    call MPI_TYPE_COMMIT(type_colonne_1_1,ierr)
!!$    call MPI_TYPE_CONTIGUOUS(imax+2,mpi_double_precision,type_colonne_0_x,ierr)
!!$    call MPI_TYPE_COMMIT(type_colonne_0_x,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ ligne tab(:,jj,kk,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_VECTOR(jmax-2,1,imax+2, mpi_double_precision , type_ligne_1_0, ierr)
!!$!    call MPI_TYPE_VECTOR(jmax+2,1,imax+2, mpi_double_precision , type_ligne_1_0, ierr)
!!$    !! Ligne d'une tableau tous les i=1, imax, avec indice j constant
!!$    call MPI_TYPE_COMMIT(type_ligne_1_0,ierr)
!!$    call MPI_TYPE_VECTOR(jmax,1,imax+2, mpi_double_precision , type_ligne_1_1, ierr)
!!$    call MPI_TYPE_COMMIT(type_ligne_1_1,ierr)
!!$    call MPI_TYPE_VECTOR(jmax+2,1,imax+2, mpi_double_precision , type_ligne_0_y, ierr)
!!$    call MPI_TYPE_COMMIT(type_ligne_0_y,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ Plan (oik) tab(:,jj,::,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_VECTOR(kmax,imax-2,(imax+2)*(jmax+2), mpi_double_precision , type_plan_oik_1_0, ierr)
!!$!    call MPI_TYPE_VECTOR(kmax,imax+2,(imax+2)*(jmax+2), mpi_double_precision , type_plan_oik_1_0, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_1_0,ierr)
!!$    call MPI_TYPE_VECTOR(kmax,imax,(imax+2)*(jmax+2), mpi_double_precision , type_plan_oik_1_1, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_1_1,ierr)
!!$    call MPI_TYPE_VECTOR(kmax,imax+2,(imax+2)*(jmax+2), mpi_double_precision , type_plan_oik_0_x, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_0_x,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ Plan (ojk) tab(ii,:,:,itps) +++++++++++++++++++++++++++
!!$    call mpi_type_hvector(kmax,1,((imax+2)*(jmax+2)*8),type_ligne_1_0,type_plan_ojk_1_0,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_1_0,ierr)
!!$    call mpi_type_hvector(kmax,1,((imax+2)*(jmax+2)*8),type_ligne_1_1,type_plan_ojk_1_1,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_1_1,ierr)
!!$    call mpi_type_hvector(kmax,1,((imax+2)*(jmax+2)*8),type_ligne_0_y,type_plan_ojk_0_y,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_0_y,ierr)
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ ppur les coins
   !!$!+++++ Coin tab(ii,jj,:,itps) +++++++++++++++++++++++++++
!!$    call mpi_type_vector(kmax,1,(imax+2)*(jmax+2),mpi_double_precision,type_coin_1,ierr)
!!$    call MPI_TYPE_COMMIT(type_coin_1,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$
!!$
!!$
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$! pour les tableau   (-1:imax+2,-1:jmax+2,0:kmax+1) 
!!$!+++++ Colonne tab(ii,:,kk,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_CONTIGUOUS(imax-2,mpi_double_precision,type_colonne_2_0,ierr)
!!$    call MPI_TYPE_COMMIT(type_colonne_2_0,ierr)
!!$    call MPI_TYPE_CONTIGUOUS(imax,mpi_double_precision,type_colonne_2_1,ierr)
!!$    call MPI_TYPE_COMMIT(type_colonne_2_1,ierr)
!!$    call MPI_TYPE_CONTIGUOUS(imax+2,mpi_double_precision,type_colonne_2_2,ierr)
!!$    call MPI_TYPE_COMMIT(type_colonne_2_2,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ ligne tab(:,jj,kk,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_VECTOR(jmax-2,1,imax+4, mpi_double_precision , type_ligne_2_0, ierr)
!!$    !! Ligne d'une tableau tous les i=1, imax, avec indice j constant
!!$    call MPI_TYPE_COMMIT(type_ligne_2_0,ierr)
!!$    call MPI_TYPE_VECTOR(jmax,1,imax+4, mpi_double_precision , type_ligne_2_1, ierr)
!!$    call MPI_TYPE_COMMIT(type_ligne_2_1,ierr)
!!$    call MPI_TYPE_VECTOR(jmax+2,1,imax+4, mpi_double_precision , type_ligne_2_2, ierr)
!!$    call MPI_TYPE_COMMIT(type_ligne_2_2,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ Plan (oik) tab(:,jj,::,itps) +++++++++++++++++++++++++++
!!$    call MPI_TYPE_VECTOR(kmax,imax-2,(imax+4)*(jmax+4), mpi_double_precision , type_plan_oik_2_0, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_2_0,ierr)
!!$    call MPI_TYPE_VECTOR(kmax,imax,(imax+4)*(jmax+4), mpi_double_precision , type_plan_oik_2_1, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_2_1,ierr)
!!$    call MPI_TYPE_VECTOR(kmax,imax+2,(imax+4)*(jmax+4), mpi_double_precision , type_plan_oik_2_2, ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_oik_2_2,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ Plan (ojk) tab(ii,:,:,itps) +++++++++++++++++++++++++++
!!$    call mpi_type_hvector(kmax,1,((imax+4)*(jmax+4)*8),type_ligne_2_0,type_plan_ojk_2_0,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_2_0,ierr)
!!$    call mpi_type_hvector(kmax,1,((imax+4)*(jmax+4)*8),type_ligne_2_1,type_plan_ojk_2_1,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_2_1,ierr)
!!$    call mpi_type_hvector(kmax,1,((imax+4)*(jmax+4)*8),type_ligne_2_2,type_plan_ojk_2_2,ierr)
!!$    call MPI_TYPE_COMMIT(type_plan_ojk_2_2,ierr)
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$!+++++ ppur les coins
   !!$!+++++ Coin tab(ii,jj,:,itps) +++++++++++++++++++++++++++
!!$    call mpi_type_vector(kmax,1,(imax+4)*(jmax+4),mpi_double_precision,type_coin_2,ierr)
!!$    call MPI_TYPE_COMMIT(type_coin_2,ierr)
!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!$
!!$!---------------------------------------------------------------------
!!$! Fin cCreation des types MPI d'echange
!!$!---------------------------------------------------------------------
!!$  end subroutine create_mpi_type
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      subroutine echange4dr8_voisin_itps(typeC,tab,lb,ub,itps,voisin)                !04-07-10    
!!$        character(len=2),intent(in) :: typeC
!!$        integer,dimension(4),intent(in) :: lb,ub
!!$        integer,intent(in) :: itps, voisin
!!$        double precision,                                                 &
!!$             dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)),  &
!!$             intent(inout) :: tab
!!$        !
!!$        integer,dimension(2,4) :: pos
!!$        integer :: C1voisin, C2voisin
!!$        integer,dimension(4) :: PosS,PosR
!!$        integer,dimension(4) :: PosC1S,PosC1R
!!$        integer,dimension(4) :: PosC2S,PosC2R
!!$        integer :: type,direction,directionC1,directionC2
!!$        integer :: tag1, tag2, imax, jmax, kmaxp1
!!$        integer :: type_plan_oik, type_plan_ojk, type_coin
!!$        integer :: debi,debj,fini,finj,tagr,nco,bcl
!!$        !
!!$        imax=par%imax
!!$        jmax=par%jmax
!!$        kmaxp1=par%kmaxp1
!!$        debi = 2
!!$        debj = 2
!!$        tagr=0
!!$        ! le bon type
!!$        if ( (lb(1)==-1).and.(ub(1)==imax+2).and. &
!!$             (lb(2)==-1).and.(ub(2)==jmax+2).and. &
!!$             (lb(3)==0).and.(ub(3)==kmaxp1) ) then
!!$           type_plan_oik=type_plan_oik_2_0
!!$           type_plan_ojk=type_plan_ojk_2_0
!!$           type_coin=type_coin_2
!!$        endif
!!$        if ( (lb(1)==0).and.(ub(1)==imax+1).and. &
!!$             (lb(2)==0).and.(ub(2)==jmax+1).and. &
!!$             (lb(3)==0).and.(ub(3)==kmaxp1) ) then
!!$           type_plan_oik=type_plan_oik_1_0
!!$           type_plan_ojk=type_plan_ojk_1_0
!!$           type_coin=type_coin_1
!!$        endif
!!$        pos = get_pos(typeC)
!!$!Le cas des echanges x et y (x0 et y0)
!!$        if ((typeC(1:1) == 'x')) then
!!$           if ((voisin == par%tvoisin(ouest)).or.(voisin == par%tvoisin(est))) then
!!$              type_plan_ojk=type_plan_ojk_0_y
!!$!              type_plan_ojk=type_plan_ojk_1_1
!!$              debi=0
!!$              debj=0
!!$              tagr=1
!!$              !pos = get_pos("z1")
!!$          endif
!!$        endif
!!$        if ((typeC(1:1) == 'y')) then
!!$           if ((voisin == par%tvoisin(nord)).or.(voisin == par%tvoisin(sud))) then
!!$              type_plan_oik=type_plan_oik_0_x
!!$!              type_plan_ojk=type_plan_oik_1_1
!!$              debi=0
!!$              debj=0
!!$              tagr=1
!!$              !pos = get_pos("z1")
!!$           endif
!!$        endif
!!$
!!$        ! les positions des echanges
!!$        if (voisin == par%tvoisin(ouest)) then
!!$           direction = ouest
!!$           PosR = (/ pos(ouestest,4), debj, lb(3), itps /)
!!$           PosS = (/ pos(ouestest,1), debj, lb(3), itps /)
!!$           type=type_plan_ojk
!!$           !COINS
!!$           directionC1=nordouest
!!$           C1voisin = par%tvoisin(nordouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,2), lb(3), itps /) !.........
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,3), lb(3), itps /)
!!$        endif
!!$        if (voisin == par%tvoisin(est)) then
!!$           direction = est
!!$           PosR = (/ pos(ouestest,2), debj, lb(3), itps /)
!!$           PosS = (/ pos(ouestest,3), debj, lb(3), itps /)
!!$           type=type_plan_ojk
!!$           !COINS
!!$           directionC1=sudest
!!$           C1voisin = par%tvoisin(sudest)
!!$           PosC1R = (/ pos(ouestest,2),  pos(nordsud,4), lb(3), itps /)
!!$           PosC1S = (/ pos(ouestest,3),  pos(nordsud,1), lb(3), itps /) !........
!!$        endif
!!$        if (voisin == par%tvoisin(nord)) then
!!$           direction = nord
!!$           PosR = (/ debi, pos(nordsud,2), lb(3), itps /)
!!$           PosS = (/ debi, pos(nordsud,3), lb(3), itps /)
!!$           type=type_plan_oik
!!$           !COINS
!!$           directionC1=nordest
!!$           C1voisin = par%tvoisin(nordest)
!!$           PosC1R = (/  pos(ouestest,2), pos(nordsud,2), lb(3), itps /)
!!$           PosC1S = (/  pos(ouestest,3), pos(nordsud,3), lb(3), itps /)
!!$        endif
!!$        if (voisin == par%tvoisin(sud)) then
!!$           direction = sud
!!$           PosR = (/ debi, pos(nordsud,4), lb(3), itps /)
!!$           PosS = (/ debi, pos(nordsud,1), lb(3), itps /)
!!$           type=type_plan_oik
!!$           !COINS
!!$           directionC1=sudouest
!!$           C1voisin = par%tvoisin(sudouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,4), lb(3), itps /)
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,1), lb(3), itps /)
!!$        endif
!!$        
!!$        tag1 = tagSend(direction)!+tagr
!!$        tag2 = tagRecv(direction)!+tagr
!!$        ! Echange les frontières
!!$        CALL sendredcv4dr8(tab,PosS,PosR,type, voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        ! Echange les coins.....................................................
!!$      end subroutine echange4dr8_voisin_itps
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      subroutine echange3dr8_voisin(typeC,tab,lb,ub,voisin)                !04-07-10    
!!$        character(len=2),intent(in) :: typeC
!!$        integer,dimension(3),intent(in) :: lb,ub
!!$        integer,intent(in) :: voisin
!!$        double precision,                                                 &
!!$             dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),  &
!!$             intent(inout) :: tab
!!$        !
!!$        integer,dimension(2,4) :: pos,pos2
!!$        integer :: C1voisin, C2voisin
!!$        integer,dimension(3) :: PosS,PosR
!!$        integer,dimension(3) :: PosC1S,PosC1R
!!$        integer,dimension(3) :: PosC2S,PosC2R
!!$        integer :: type,direction,directionC1,directionC2
!!$        integer :: tag1, tag2, imax, jmax, kmaxp1
!!$        integer :: type_plan_oik, type_plan_ojk, type_coin
!!$        integer :: debi,debj,fini,finj,tagr,nco,bcl
!!$        !
!!$
!!$        imax=par%imax
!!$        jmax=par%jmax
!!$        kmaxp1=par%kmaxp1
!!$        debi = 2
!!$        debj = 2
!!$        tagr=0
!!$        ! le bon type
!!$        if ( (lb(1)==-1).and.(ub(1)==imax+2).and. &
!!$             (lb(2)==-1).and.(ub(2)==jmax+2).and. &
!!$             (lb(3)==0).and.(ub(3)==kmaxp1) ) then
!!$           type_plan_oik=type_plan_oik_2_0
!!$           type_plan_ojk=type_plan_ojk_2_0
!!$!!Methode 2
!!$           if ((typeC(2:2) == '1')) then
!!$              type_plan_oik=type_plan_oik_2_1
!!$              type_plan_ojk=type_plan_ojk_2_1
!!$              debi=1
!!$              debj=1
!!$           endif
!!$           if ((typeC(2:2) == '2')) then
!!$              type_plan_oik=type_plan_oik_2_2
!!$              type_plan_ojk=type_plan_ojk_2_2
!!$              debi=0
!!$              debj=0
!!$           endif
!!$           type_coin=type_coin_2
!!$        endif
!!$        if ( (lb(1)==0).and.(ub(1)==imax+1).and. &
!!$             (lb(2)==0).and.(ub(2)==jmax+1).and. &
!!$             (lb(3)==0).and.(ub(3)==kmaxp1) ) then
!!$           type_plan_oik=type_plan_oik_1_0
!!$           type_plan_ojk=type_plan_ojk_1_0
!!$!!Methode 2
!!$           if ((typeC(2:2) == '1')) then
!!$              type_plan_oik=type_plan_oik_1_1
!!$              type_plan_ojk=type_plan_ojk_1_1
!!$              debi=1
!!$              debj=1
!!$           endif
!!$           type_coin=type_coin_1
!!$        endif
!!$        pos = get_pos(typeC)
!!$!Le cas des echanges x et y (x0 et y0)
!!$        if ((typeC(1:1) == 'x')) then
!!$           if ((voisin == par%tvoisin(ouest)).or.(voisin == par%tvoisin(est))) then
!!$              type_plan_ojk=type_plan_ojk_0_y
!!$!              type_plan_ojk=type_plan_ojk_1_1
!!$              debi=0
!!$              debj=0
!!$              tagr=1
!!$              !pos = get_pos("z1")
!!$          endif
!!$        endif
!!$        if ((typeC(1:1) == 'y')) then
!!$           if ((voisin == par%tvoisin(nord)).or.(voisin == par%tvoisin(sud))) then
!!$              type_plan_oik=type_plan_oik_0_x
!!$!              type_plan_oik=type_plan_oik_1_1
!!$              debi=0
!!$              debj=0
!!$              tagr=1
!!$              !pos = get_pos("z1")
!!$           endif
!!$        endif
!!$
!!$        ! les positions des echanges
!!$        if (voisin == par%tvoisin(ouest)) then
!!$           direction = ouest
!!$           PosR = (/ pos(ouestest,4), debj, lb(3) /)
!!$           PosS = (/ pos(ouestest,1), debj, lb(3) /)
!!$         type=type_plan_ojk
!!$           !COINS
!!$           directionC1=nordouest
!!$           C1voisin = par%tvoisin(nordouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,2), lb(3) /) !.........
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,3), lb(3) /)
!!$        endif
!!$        if (voisin == par%tvoisin(est)) then
!!$           direction = est
!!$           PosR = (/ pos(ouestest,2), debj, lb(3) /)
!!$           PosS = (/ pos(ouestest,3), debj, lb(3) /)
!!$           type=type_plan_ojk
!!$           !COINS
!!$           directionC1=sudest
!!$           C1voisin = par%tvoisin(sudest)
!!$           PosC1R = (/ pos(ouestest,2),  pos(nordsud,4), lb(3) /)
!!$           PosC1S = (/ pos(ouestest,3),  pos(nordsud,1), lb(3) /) !........
!!$        endif
!!$        if (voisin == par%tvoisin(nord)) then
!!$           direction = nord
!!$           PosR = (/ debi, pos(nordsud,2), lb(3) /)
!!$           PosS = (/ debi, pos(nordsud,3), lb(3) /)
!!$           type=type_plan_oik
!!$           !COINS
!!$           directionC1=nordest
!!$           C1voisin = par%tvoisin(nordest)
!!$           PosC1R = (/  pos(ouestest,2), pos(nordsud,2), lb(3) /)
!!$           PosC1S = (/  pos(ouestest,3), pos(nordsud,3), lb(3) /)
!!$        endif
!!$        if (voisin == par%tvoisin(sud)) then
!!$           direction = sud
!!$           PosR = (/ 2, pos(nordsud,4), lb(3) /)
!!$           PosS = (/ 2, pos(nordsud,1), lb(3) /)
!!$          type=type_plan_oik
!!$           !COINS
!!$           directionC1=sudouest
!!$           C1voisin = par%tvoisin(sudouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,4), lb(3) /)
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,1), lb(3) /)
!!$        endif
!!$        
!!$        tag1 = tagSend(direction)!+tagr
!!$        tag2 = tagRecv(direction)!+tagr
!!$        ! Echange les frontières
!!$        CALL sendredcv3dr8(tab,PosS,PosR,type, voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        ! Echange les coins
!!$        tag1 = tagSend(directionC1)
!!$        tag2 = tagRecv(directionC1)
!!$        CALL sendredcv3dr8(tab,PosC1S,PosC1R,type_coin, C1voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        !if ((typeC(1:1) == 'z')) then
!!$           tag1 = tagSend(directionC1)
!!$           tag2 = tagRecv(directionC1)
!!$           CALL sendredcv3dr8(tab,PosC1S,PosC1R,type_coin, C1voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        !endif
!!$       
!!$      end subroutine echange3dr8_voisin
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      subroutine echange3dr8_voisin_itps(typeC,tab,lb,ub,itps,voisin)                !04-07-10    
!!$        character(len=2),intent(in) :: typeC
!!$        integer,dimension(3),intent(in) :: lb,ub
!!$        integer,intent(in) :: itps, voisin
!!$        double precision,                                                 &
!!$             dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),  &
!!$             intent(inout) :: tab
!!$        !
!!$        integer,dimension(2,4) :: pos
!!$        integer :: C1voisin, C2voisin
!!$        integer,dimension(3) :: PosS,PosR
!!$        integer,dimension(3) :: PosC1S,PosC1R
!!$        integer,dimension(3) :: PosC2S,PosC2R
!!$        integer :: type,direction,directionC1,directionC2
!!$        integer :: tag1, tag2, imax, jmax, kmaxp1
!!$        integer :: type_plan_oik, type_plan_ojk, type_coin
!!$        integer :: debi,debj,fini,finj,tagr,nco,bcl
!!$        !
!!$        imax=par%imax
!!$        jmax=par%jmax
!!$        kmaxp1=par%kmaxp1
!!$        debi = 2
!!$        debj = 2
!!$        ! le bon type
!!$        if ( (lb(1)==-1).and.(ub(1)==imax+2).and. &
!!$             (lb(2)==-1).and.(ub(2)==jmax+2) ) then
!!$           type_plan_oik=type_colonne_2_0
!!$           type_plan_ojk=type_ligne_2_0
!!$!!Methode 2
!!$           if ((typeC(2:2) == '1')) then
!!$              type_plan_oik=type_colonne_2_1
!!$              type_plan_ojk=type_ligne_2_1
!!$              debi=1
!!$              debj=1
!!$           endif
!!$           if ((typeC(2:2) == '2')) then
!!$              type_plan_oik=type_colonne_2_2
!!$              type_plan_ojk=type_ligne_2_2
!!$              debi=0
!!$              debj=0
!!$           endif
!!$           type_coin=mpi_double_precision  
!!$        endif
!!$        if ( (lb(1)==0).and.(ub(1)==imax+1).and. &
!!$             (lb(2)==0).and.(ub(2)==jmax+1) ) then
!!$           type_plan_oik=type_colonne_1_0
!!$           type_plan_ojk=type_ligne_1_0
!!$           type_coin=mpi_double_precision
!!$!!Methode 2
!!$           if ((typeC(2:2) == '1')) then
!!$              type_plan_oik=type_colonne_1_1
!!$              type_plan_ojk=type_ligne_1_1
!!$              debi=1
!!$              debj=1
!!$           endif
!!$        endif
!!$        pos = get_pos(typeC)
!!$!Le cas des echanges x et y (x0 et y0)
!!$        
!!$        if ((typeC(1:1) == 'x')) then        
!!$           type_plan_oik=type_plan_oik_0_x
!!$           type_plan_ojk=type_plan_ojk_0_y
!!$           debi=0
!!$           debj=0
!!$           tagr=1
!!$        endif
!!$
!!$      ! les positions des echanges
!!$        if (voisin == par%tvoisin(ouest)) then
!!$           direction = ouest
!!$           PosR = (/ pos(ouestest,4), debj, itps /)
!!$           PosS = (/ pos(ouestest,1), debj, itps /)
!!$           type=type_plan_ojk
!!$           !COINS
!!$           directionC1=nordouest
!!$           C1voisin = par%tvoisin(nordouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,2), itps /) !.........
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,3), itps /)
!!$        endif
!!$        if (voisin == par%tvoisin(est)) then
!!$           direction = est
!!$           PosR = (/ pos(ouestest,2), debj, itps /)
!!$           PosS = (/ pos(ouestest,3), debj, itps /)
!!$           type=type_plan_ojk
!!$           !COINS
!!$           directionC1=sudest
!!$           C1voisin = par%tvoisin(sudest)
!!$           PosC1R = (/ pos(ouestest,2),  pos(nordsud,4), itps /)
!!$           PosC1S = (/ pos(ouestest,3),  pos(nordsud,1), itps /) !........
!!$        endif
!!$        if (voisin == par%tvoisin(nord)) then
!!$           direction = nord
!!$           PosR = (/ debi, pos(nordsud,2), itps /)
!!$           PosS = (/ debi, pos(nordsud,3), itps /)
!!$           type=type_plan_oik
!!$           !COINS
!!$           directionC1=nordest
!!$           C1voisin = par%tvoisin(nordest)
!!$           PosC1R = (/  pos(ouestest,2), pos(nordsud,2), itps /)
!!$           PosC1S = (/  pos(ouestest,3), pos(nordsud,3), itps /)
!!$        endif
!!$        if (voisin == par%tvoisin(sud)) then
!!$           direction = sud
!!$           PosR = (/ debi, pos(nordsud,4), itps /)
!!$           PosS = (/ debi, pos(nordsud,1), itps /)
!!$           type=type_plan_oik
!!$           !COINS
!!$           directionC1=sudouest
!!$           C1voisin = par%tvoisin(sudouest)
!!$           PosC1R = (/ pos(ouestest,4), pos(nordsud,4), itps /)
!!$           PosC1S = (/ pos(ouestest,1), pos(nordsud,1), itps /)
!!$        endif
!!$        
!!$        tag1 = tagSend(direction)
!!$        tag2 = tagRecv(direction)
!!$        ! Echange les frontières
!!$        write(5000+par%rank,*) typeC(1:1),".",int(voisin,kind=1)," <=> Env :",int(PosS(1),kind=1),int(PosS(2),kind=1), &
!!$             " Rec : ",int(PosR(1),kind=1),int(PosR(2),kind=1)," e:",int(tab(posS(1),posS(2),posS(3)),1), &
!!$             " debi,debj=",int(debi,kind=1),int(debj,kind=1)
!!$        CALL sendredcv3dr8(tab,PosS,PosR,type, voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        ! Echange les coins
!!$        !if ((typeC(1:1) == 'z')) then
!!$           tag1 = tagSend(directionC1)
!!$           tag2 = tagRecv(directionC1)
!!$           CALL sendredcv3dr8(tab,PosC1S,PosC1R,type_coin, C1voisin,lbound(tab),ubound(tab),tag1,tag2)
!!$        !endif
!!$      
!!$      end subroutine echange3dr8_voisin_itps
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine petite_pause(tps)
      implicit none
      integer,intent(in) :: tps
      logical :: fin=.false.
      real :: time0,time1
      call cpu_time(time0)
      do while(.not.fin)
         call cpu_time(time1)
         if ( (time1-time0) > tps) fin=.true.
      end do
      end subroutine petite_pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      end module module_parallele
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else
!!!!!!   VERSION SEQUENTIELLE 

!---------------------------------------------------------------------
!Version NON PARALLELE du module parallele
         module module_parallele
      implicit none
      integer, parameter :: ndims=2
      integer, parameter :: ouest=1,est=2,nord=3,sud=4,haut=5,bas=6
      integer, parameter :: sudouest=7,sudest=8,nordouest=9,nordest=10,ouest2=11,est2=12 
      integer, parameter :: ouestest=1,nordsud=2
      integer, parameter :: nbvoisin=10
!---------------------------------------------------------------------
!     Informations paralleles contenues dans une structure
!--------------------------------------------------------------------
!
      type infopar
         integer                                    ::  rank            !RANG de PROCESSUS
         integer                                    ::  subcycle=1      !Nombre de cycles dans un pas de temps principal
         integer                                    ::  nbdom           !NOMBRE DE DOMAINES
         integer                                    ::  dimax,djmax     !NOMBRE DE DOMAINES
         integer                                    ::  iglb,jglb,imax,jmax       !TAILLE GLOBALE
         integer,dimension(2)                       ::  coords          !COORDONNEES DANS LA GRILLE DE PROCESSUS
         integer,dimension(2)                       ::  timax           !OFFSET DEBUT,FIN
         integer,dimension(2)                       ::  tjmax           !DANS LA NUMEROTATION GLOBALE
         integer,dimension(:,:),allocatable         ::  gtimax          !INDICES DEBUT,FIN
         integer,dimension(:,:),allocatable         ::  gtjmax          !DANS LA NUMEROTATION GLOBALE
         integer,dimension(4)                       ::  tvoisin         !LE NUMERO DES VOISINS DANS L'ORDRE(O,E,S,N)
         integer,dimension(:,:),allocatable         ::  gtvoisin        !LE NUMERO DES VOISINS DANS L'ORDRE(O,E,S,N)

      end type infopar
      integer                         :: ierr
      character(len=6)                :: dom_c            !CHAINE DE CARACTER REPRODUISANT LES COORDONNEES
      type (infopar)                  :: par
      integer, parameter   :: ki4=selected_int_kind(4)
      integer, parameter   :: ki8=selected_int_kind(8)
      integer, parameter   :: kr4=selected_real_kind(6,37)
      integer, parameter   :: kr8=selected_real_kind(15,307)
      integer,parameter    :: mpi_proc_null=-2

      contains
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      subroutine initialise_parallel(dimax,djmax,iglb,jglb,             &
           imax,jmax,kmaxp1,ipbc,jpbc)
!---------------------------------------------------------------------
!      INITIALISATION DES PAREMTRES DE PARALELLISATION :
!         RANK            :-> Numero du processus courant
!         NBDOM           :-> Nombre de Domaines
!         COMM2D          :-> Communicateur intra grille de processus
!         COORDS          :-> Coordonees de process dans la grille
!         TMECO           :-> Offset pour retrouver les indices globaux
!         TNECO           :-> Offset pour retrouver les indices globaux
!         TVOISIN         :-> Numero des processus voisins O,E,N,S
!---------------------------------------------------------------------
       !IMPLICIT NONE
       integer,intent(in)    :: dimax,djmax
       integer,intent(in)    :: iglb,jglb,kmaxp1,imax,jmax
       logical,intent(in)    :: ipbc,jpbc                              !15-10-10
!      integer               :: imax,jmax
!      integer,intent(out)   :: imax,jmax
       integer,dimension(:,:),allocatable         ::  gtimax          !INDICES DEBUT,FIN
       integer,dimension(:,:),allocatable         ::  gtjmax          !DANS LA NUMEROTATION GLOBALE
!
       integer :: nbdom
       integer,dimension(ndims) :: dims
       integer,dimension(2)     :: coords,timax,tjmax
       logical,dimension(ndims) :: periods
       character :: c,c1,c2
       character(len=6)         :: dom_c

! On initialise toutes les varible utiles a la version non parallele.
       par%rank=0
       par%nbdom=1
       par%dimax=1
       par%djmax=1
       par%iglb=iglb
       par%jglb=jglb
       par%coords(1)  = 0
       par%coords(2)  = 0

       par%timax(1) = 0
       par%timax(2) = iglb

       par%tjmax(1) = 0
       par%tjmax(2) = jglb

       par%tvoisin = -2

!      imax=iglb
!      jmax=jglb

      allocate(par%gtimax(0:0,1))
      allocate(par%gtjmax(0:0,1))
      par%gtimax(0,1)=0
      par%gtjmax(0,1)=0

!---------------------------------------------------------------------
!          FIN INIT
!---------------------------------------------------------------------
      end subroutine initialise_parallel
!---------------------------------------------------------------------

      end module module_parallele

#endif
