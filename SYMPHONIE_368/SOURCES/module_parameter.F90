      module module_parameter
      implicit none                                      !
!______________________________________________________________________
! SYMPHONIE ocean model
! release 354 - last update: 12-09-22
!______________________________________________________________________
! Version Date      Description des modifications
! S26     17-09-13  ifb et dim_airsea passent dans module_principal
!         23-07-14  suppression nfichier
!         09-12-14  dim_river lu dans notebook_river
!         03-12-15  par defaut numbyears=200
!         16-05-15  nbincomax=24
!         03-09-18  dbefore=-1,dnow=0,dafter=1
! v253    01-05-19  any2=8
! v269    05-12-19  any2=9
! v287    17-08-20  ip_i4_p,ip_i8_p,ip_r16_p
! v292    16-11-20  les dimensions des drifters sont declarees dans module_principal 
!                   et fixees dans module_drifter
! v296    12-02-21  numbyears=250 nbincomax=30
! v354    12-09-22  numbyears=350
!...............................................................................
!    _________                    .__                  .__             ! (°o°)
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      !
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................

      integer, parameter :: ip_i4_p  = selected_int_kind(9)   ! integer 4
      integer, parameter :: ip_i8_p  = selected_int_kind(15)  ! integer 8
      integer, parameter :: ip_r16_p = selected_real_kind(18) ! real 16


      integer numbyears,dim_river

!________________________________________________________
! The principal parameters concerning the grid size
      integer iglb,jglb,kmax,nbdom_imax,nbdom_jmax,imax,jmax,kmaxp1  &
             ,nbdom
      logical iperiodicboundary,jperiodicboundary
!________________________________________________________


!________________________________________________________
! ---------------- TIME STEPPING       ------------------
      integer before,bef,now,after,aft,before2ext,before2,before3 &
!            ,dbefore,dnow,dafter                                 &
             ,beforebefore
      parameter(before=0,bef=before,now=1,after=2,aft=after       &
!              ,dbefore=-1,dnow=0,dafter=1                        &
               ,before2ext=2,before2=-1,before3=2,beforebefore=2  )

!________________________________________________________
                                                        !
!________________________________________________________
!   DIMENSIONS PARAMETRISATION DE L'ETAT DE REFERENCE
      integer dim_ref
      parameter(dim_ref=6)

!________________________________________________________
! NOMBRE D'ANNEES MAXIMUM DANS UNE SIMULATION
      parameter(numbyears=350) !12-09-22
!________________________________________________________

!_______________________________________________________________________________
! Nombre de variables dans le modele de transport par advection diffusion
! STRADA: DIM_BIO et nombre de points source pour des tests de dispersion
! de polluants: DIM_SOURCE
      integer dim_timebio,bioaver_onoff,bdt1,bdt2
      parameter(                                                        &
        dim_timebio=1                            & ! 0 si schema leapfrog, 1 si schema forward
       ,bioaver_onoff=0                          & ! 1 si calcul de moyenne pour modeles filles
       ,bdt1=2*dim_timebio                       & ! Dimension du temps min
       ,bdt2=2+bioaver_onoff                     & ! Dimension du temps max
       )
!_______________________________________________________________________________

!_______________________________________________________________________________
! Nombre maxi d'inconnues pour resolution d'un (petit) systeme d'equation
! lineaire (surdetermine.F)
      integer nbincomax
      parameter(nbincomax=30)
!_______________________________________________________________________________

!_______________________________________________________________________________
! Nombre Maximum de bouees lagrangiennes
!     integer nbomax,nbobuffermax,nbobuffermax_c
!     parameter(nbomax=14000,nbobuffermax=140,nbobuffermax_c=14)
!_______________________________________________________________________________

!_______________________________________________________________________________
! Nombre maximum de dates dans le fichier notebook_dateoutput
      integer dim_dof
      parameter(dim_dof=10)
!_______________________________________________________________________________

!_______________________________________________________________________________
! Nombre maximum cumulé d'emplacements disponibles dans notebook_graph
! & dimension maximum pour titres de variables visualisées
      integer kmax_grh,dimgrh_titre
      parameter(kmax_grh=95                                             &
               ,dimgrh_titre=100) 

!_______________________________________________________________________________

!_______________________________________________________________________________
! Activer dimension du tableau de moyenne temporelle du courant
      integer mi_onoff,m_mi,n_mi,nr_mi
      parameter(                                                        &
        mi_onoff=0 )
!_______________________________________________________________________________


!_______________________________________________________________________________
! Dimensions des tableaux impliqués dans la filière d'imbrication de modèles
! Symphonie Maman - Symphonie Filles:
      integer nest_max,nest_m,nest_n,nest_r                             &
                      ,nest_dim0,nest_dim1,nest_dim2                    &
                      ,nest_dim3,nest_dim4,nest_dim5                    &
                      ,nest_dim6,nest_dim7
      parameter(nest_max=1                                              &
               ,nest_dim0=0                                             &
               ,nest_dim1=1                                             &
               ,nest_dim2=1                                             &
               ,nest_dim3=1                                             &
               ,nest_dim4=1                                             &
               ,nest_dim5=1                                             &
               ,nest_dim6=1                                             &
               ,nest_dim7=1                                             &
               )
!_______________________________________________________________________________

!_______________________________________________________________________________
! Dimension des variables de modele_wave
      integer                                                           &
        onoff_wave      & ! (1 ou 0) active ou desactive le modele de marée
       ,imax_w          & ! Dimension i des tableaux
       ,jmax_w          & ! Dimension j des tableaux
       ,fmax_w          &
       ,any1                                                            &
       ,any2
      parameter(                                                        &
        onoff_wave=0                                                    &
!      ,imax_w=onoff_wave*imax+(1-onoff_wave)*1                         &
!      ,jmax_w=onoff_wave*jmax+(1-onoff_wave)*1                         &
       ,fmax_w=onoff_wave*30+(1-onoff_wave)*1                           &
       ,any1=0                                                          &
       ,any2=9                                                          &
               )
!_______________________________________________________________________________

!     integer ipt,imt,jpt,jmt,ipu,imu,jpu,jmu,ipv,imv,jpv,jmv,kpt,kmt,kpw,kmw
      integer ipt,imt,jpt,jmt,ipu,imu                ,jpv,jmv,kpt,kmt,kpw,kmw
      parameter( &
      ipt=1  ,   & ! i forward  intrement relative to "t" points
      imt=0  ,   & ! i backward intrement relative to "t" points
      jpt=1  ,   & ! j forward  intrement relative to "t" points
      jmt=0  ,   & ! j backward intrement relative to "t" points
      ipu=0  ,   & ! i forward  intrement relative to "u" points          , : !06-12-09
      imu=1  ,   & ! i backward intrement relative to "u" points
!     jpu=0  ,   & ! j forward  intrement relative to "u" points
!     jmu=1  ,   & ! j backward intrement relative to "u" points
!     ipv=0  ,   & ! i forward  intrement relative to "v" points
!     imv=1  ,   & ! i backward intrement relative to "v" points
      jpv=0  ,   & ! j forward  intrement relative to "v" points
      jmv=1  ,   & ! j backward intrement relative to "v" points
      kpt=1  ,   & ! k upward   intrement relative to "t" points
      kmt=0  ,   & ! k downward intrement relative to "t" points
      kpw=0  ,   & ! k upward   intrement relative to "w" points
      kmw=1      ) ! k downward intrement relative to "w" points

!--------------------------------------------------------------------------------
! Coefficients of the Jackett et al. (2006) equation of state:
      double precision                         &
           c1_jmfwg                           &
          ,c2_jmfwg                           &
          ,c3_jmfwg                           &
          ,c4_jmfwg                           &
          ,c5_jmfwg                           &
          ,c6_jmfwg                           &
          ,c7_jmfwg                           &
          ,c8_jmfwg                           &
          ,c9_jmfwg                           &
          ,c10_jmfwg                          &
          ,c11_jmfwg                          &
          ,c12_jmfwg                          &
          ,c14_jmfwg                          &
          ,c15_jmfwg                          &
          ,c16_jmfwg                          &
          ,c17_jmfwg                          &
          ,c18_jmfwg                          &
          ,c19_jmfwg                          &
          ,c20_jmfwg                          &
          ,c21_jmfwg                          &
          ,c22_jmfwg                          &
          ,c23_jmfwg                          &
          ,c24_jmfwg                          &
          ,c25_jmfwg
      parameter(                              &
           c1_jmfwg=  9.9984085444849347d2    &
          ,c2_jmfwg=  7.3471625860981584d0    &
          ,c3_jmfwg= -5.3211231792841769d-2   &
          ,c4_jmfwg=  3.6492439109814549d-4   &
          ,c5_jmfwg=  2.5880571023991390d0    &
          ,c6_jmfwg= -6.7168282786692355d-3   &
          ,c7_jmfwg=  1.9203202055760151d-3   &
          ,c8_jmfwg=  1.1798263740430364d-2   &
          ,c9_jmfwg=  9.8920219266399117d-8   &
          ,c10_jmfwg= 4.6996642771754730d-6   &
          ,c11_jmfwg=-2.5862187075154352d-8   &
          ,c12_jmfwg=-3.2921414007960662d-12  &
          ,c14_jmfwg= 7.2815210113327091d-3   &
          ,c15_jmfwg=-4.4787265461983921d-5   &
          ,c16_jmfwg= 3.3851002965802430d-7   &
          ,c17_jmfwg= 1.3651202389758572d-10  &
          ,c18_jmfwg= 1.7632126669040377d-3   &
          ,c19_jmfwg=-8.8066583251206474d-6   &
          ,c20_jmfwg=-1.8832689434804897d-10  &
          ,c21_jmfwg= 5.7463776745432097d-6   &
          ,c22_jmfwg= 1.4716275472242334d-9   &
          ,c23_jmfwg= 6.7103246285651894d-6   &
          ,c24_jmfwg=-2.4461698007024582d-17  &
          ,c25_jmfwg=-9.1534417604289062d-18 )
!--------------------------------------------------------------------------------
! Coefficients of the Wright (JAOT 1997) equation of state:
! TABLE 1 column: "Extended formula fit over the reduced range"
      real*4                                &
           a0_w97                           &
          ,a1_w97                           &
          ,a2_w97                           &
          ,b0_w97                           &
          ,b1_w97                           &
          ,b2_w97                           &
          ,b3_w97                           &
          ,b4_w97                           &
          ,b5_w97                           &
          ,c0_w97                           &
          ,c1_w97                           &
          ,c2_w97                           &
          ,c3_w97                           &
          ,c4_w97                           &
          ,c5_w97
      parameter(                            &
           a0_w97= 7.057924e-4              &
          ,a1_w97= 3.480336e-7              &
          ,a2_w97=-1.112733e-7              &
          ,b0_w97= 5.790749e+8              &
          ,b1_w97= 3.516535e+6              &
          ,b2_w97=-4.002714e+4              &
          ,b3_w97= 2.084372e+2              &
          ,b4_w97= 5.944068e+5              &
          ,b5_w97=-9.643486e+3              &
          ,c0_w97= 1.704853e+5              &
          ,c1_w97= 7.904722e+2              &
          ,c2_w97=-7.984422                 &
          ,c3_w97= 5.140652e-2              &
          ,c4_w97=-2.302158e+2              &
          ,c5_w97=-3.079464  )

      end module module_parameter
