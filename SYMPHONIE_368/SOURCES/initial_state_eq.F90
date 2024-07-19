      subroutine initial_state_eq
!______________________________________________________________________
! S model
! release S.26  - last update: 28-11-14
!______________________________________________________________________

      use module_principal
      use module_parallele !#MPI
      implicit none
#ifdef synopsis
       subroutinetitle='initial_state_eq'
       subroutinedescription= &
      'Computes some useful constants relative to the sea-water ' &
      //' density (mean sea water density rho,...)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!_______________________________________________________________________________
! Version Date      Description des modifs
!         02/10/02: creation
!         03/10/02: appel à density pour mise à jour de RHP_Z
!         04/11/03: passage d'un argument dans CALL DENSITY
!         01/01/04: blindage concernant le reset des parametres et de RHP
!         13/11/04: lecture d'un notebook_eqstate
!         28/01/05: ajout blindage: puissance 1.5 appliquée à salinite negative
!                   n'existe pas...
!         03/06/05: l'equation d'etat est linearisée autour de T et S
!                  representatives des zones où la variabilité est forte
!         25/07/05: Ajout d'un cas où les coef restent aux valeurs fixées dans
!                   set_parameters.F
!         19/11/05: lecture des parametres de l'equation d'etat lineaire dans
!                   notebook_eqstate
!         14/02/06: debug du point precedent: cas EQS_STATE3=1 etait mal traite
!         10/05/06: fonctions compatibles avec double precision
!         19/07/06: la densite de reference n'est plus la densite moyenne
!                   mais une densite ponderee par sigma...
!         29/04/09: notebook_eqstate devient NOMFICHIER(17)
!         01-06-09: parallelisation !#MPI
! 2009.3  30-09-09: suppression sighl_z
!         01-10-09: Utilisation des nouveaux facteurs d'echelle verticale
!         05-10-09: ajout d'un "ifdef parallele"
! 2010.8  09-05-10  seul le proc 0 ecrit messages
! 2010.11 19-07-10  prevoir test de conservation de la parallelisation
!         20-07-10  ajout securite en cas de grille 100% masquee
!         23-07-10  routine density renommee equation_of_state
!         25-07-10  En cas d'equation d'etat lineaire ne pas permettre d'introduire
!                   l'effet de compression
! S.26    01-12-13  notebook_eq_state lu dans set_parameters.F90
!         17-05-14  ecriture rho a l'ecran
!         19-06-14  Meme si l'effet de compression est dans l'EOS, le PGF est
!                   calcule de telle sorte que la constante rho doit etre calculee
!                   uniquement sur la base de la densite potentielle a z=0m
!         18-07-14  Conservation mpi: troncature 4eme decimale de la moyenne de rhp
!         28-11-14  remplacer sigma_w par son equivalent
!_______________________________________________________________________________

!cccccIF(EQS_STATE3.GE.2)GOTO 1000                                     !19/11/05
      if(eos_linear.ge.1)goto 1000                                     !14/02/06

!.........................................................................
! Calcul de la densité moyenne pondérée (poid fort aux fortes profondeurs)
      rho=0.
      call equation_of_state('at initial state ',now)   !23-07-10
      sum0=0.
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
      sum5=0.
      sum6=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

!     x0=  mask_c(i,j,k)*dsig_c(i,j,k)*h_z(i,j)
      x0=  mask_t(i,j,k)*dz_t(i,j,k,1)                                  & !01-10-09
        *mask_i_w(i)*mask_j_w(j) !#MPI ne pas sommer 2 fois zone commune 01-06-09
!     x1=x0*sighl_z(i,j,k)
!     x1=x0*0.5*(sigma_w(i,j,k)+sigma_w(i,j,k+1))                       !30-09-09
      x1=x0*(depth_t(i,j,k)+h_w(i,j))/hz_w(i,j,1) !remplacer sigma_w par son equivalent (z+h)/(h+ssh)!28-11-14

      x2=x0*exp( depth_t(i,j,k)/300.)

      sum0=sum0+x1
      sum1=sum1+x2
      sum2=sum2+x1*tem_t(i,j,k,1)
      sum3=sum3+x2*tem_t(i,j,k,1)
      sum4=sum4+x1*sal_t(i,j,k,1)
      sum5=sum5+x2*sal_t(i,j,k,1)
      sum6=sum6+x1*(rhp_t(i,j,k)     & ! potential density
!                  +anyv3d(i,j,k,0)  & ! part of density with compression terms
                   )

      enddo
      enddo
      enddo

!............................................................................
#ifdef parallele
! Special parallelisation: ajouter les sommes de chaque sous-domaine
      call mpi_allreduce(sum0,sum0glb,1,mpi_double_precision,  & !#MPI 01-06-09
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum5,sum5glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum6,sum6glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      sum0=sum0glb ; sum1=sum1glb ; sum2=sum2glb ; sum3=sum3glb ;
      sum4=sum4glb ; sum5=sum5glb ; sum6=sum6glb
#endif
!............................................................................

! T0 S0 RHO sont calculés pour ameliorer le rendement du couplage      !03/06/05
! interne/externe autrement dit la "moyenne" est pondérée par H
! c.a.d. priorité aux zones profondes puisque c'est la que le critere
! CFL est le plus contraignant.
! On remarque egalement la ponderation par SIGHL_Z, cela provient
! de ce que l'integrale verticale de la pression est equivalente
! a l'integrale de G*RHP*H*SIGHL*DSIG. Autrement pour bien tuer
! tout residu barotrope ce n'est pas la densite moyenne simple
! qu'il faut prendre comme reference mais une moyenne ponderee par
! la coordonnee sigma (en gros poid fort aux densites des couches
! de surfaces puisqu'elles jouent plus dans l'integrale du gradient
! barocline qui force le mode externe).
       if(sum0==0) then !20-07-10
         eos_linear=1
         goto 1000
       endif

       t0=sum2/sum0
       s0=sum4/sum0
      rho=sum6/sum0

      rho=nint(rho*1.d4)/1.d4                        !18-07-14
      if(par%rank==0)write(6,'(a,e26.18)')'rho=',rho !18-07-14

! Par contre on ne linearise pas autour des valeurs moyennes T0 et S0
! car elles ne sont pas les plus representatives des zones "dynamiques"
! où nous avons les plus forts courants. On linearire autour de valeurs
! plus proches de la surface, T0SURF & S0SURF. La consequence est que
! notre equation d'etat linearisee s'ecarte de la verite dans l'absolue
! mais les variations de densite en fonction des variations de T et S
! sont neamoins mieux representées autours des valeurs de T et S represenatives
! des zones "dynamiques" et donc les gradients de pressions sont plus
! fideles...
      t0surf=sum3/sum1
      s0surf=sum5/sum1

!..................................................
! Retrouver les coefs d'une relation lineaire en regardant
! les variations de densite données par EOS80 autour de T0
! et S0.
      tem_0d=t0surf
      sal_0d=max(s0surf-0.1,zero)
!     call density(6)
      call equation_of_state('at a single point',now) !23-07-10
      rho_tmp=rho_0d
      sal_0d=s0surf+0.1
!     call density(6)
      call equation_of_state('at a single point',now) !23-07-10
      alp_s=(rho_0d-rho_tmp)/rho/(s0surf+0.1-max(s0surf-0.1,zero))

      tem_0d=t0surf-0.5
      sal_0d=s0surf
!     call density(6)
      call equation_of_state('at a single point',now) !23-07-10
      rho_tmp=rho_0d
      tem_0d=t0surf+0.5
!     call density(6)
      call equation_of_state('at a single point',now) !23-07-10
      alp_t=-(rho_0d-rho_tmp)/1./rho
!..................................................
 1000 continue

! As potential density is computed at the end of each step, the latter
! has to be initialized now in order to compute the PGF at the beginning
! of the first step:
      call equation_of_state('potential density',now)   !23-07-10


!     call equation_of_state('compression terms',now)   !23-07-10
!     do k=1,kmax
!     do j=1,jmax
!     do i=1,imax
!     if(mask_t(i,j,k)==1)write(66,*)depth_t(i,j,k),rhp_t(i,j,k),anyv3d(i,j,k,0)
!     enddo
!     enddo
!     enddo
!     stop 'rigolo'

!..................................................
! Re-Initialiser proprement la densite:
! Utile seulement pour une visu de l'etat initial car
! on repassera dans density dès la 1ere iteration
! avant le calcul du gradient de pression:

! Cas lineaire:
!     if(eos_author.eq.0) then !eqslin eqslin eqslin >
!      call density(0)
!           Cas prise en compte de la pression:
!           if(eos_comprs.eq.1) then !eqspresprespres >
!             call density(4)
!           endif                    !eqspresprespres >
!     endif                    !eqslin eqslin eqslin >
! Cas non lineaire sans pression:
!     if(eos_author.eq.1.and.eos_comprs.eq.0)call density(3)
! Cas non lineaire + pression:
!     if(eos_author.eq.1.and.eos_comprs.eq.1)call density(2)
!..................................................

!.........................................................................

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine initial_state_eq:'
      write(3,*)
      if(eos_author.eq.0)                                               &
      write(3,*)'equation d''etat lineaire'
      if(eos_author.eq.0.and.eos_linear.eq.1)                           &
      write(3,*)'avec coef fixés dans set_parameters.f'
      if(eos_author/=0)                                               &
      write(3,*)'equation d''etat non lineaire No:',eos_author
      if(eos_comprs.eq.0)                                               &
      write(3,*)'sans prise en compte de la pression'
      if(eos_comprs.eq.1)                                               &
      write(3,*)'avec prise en compte de la pression'
      write(3,*)
      write(3,*)'rho_base   ',rho_base
      write(3,*)'rho        ',rho
      write(3,*)'alp_t_base ',alp_t_base
      write(3,*)'alp_t      ',alp_t
      write(3,*)'alp_s_base ',alp_s_base
      write(3,*)'alp_s      ',alp_s
      write(3,*)'t0_base    ',t0_base
      write(3,*)'t0surf    ',t0surf
      write(3,*)'t0         ',t0
      write(3,*)'s0_base    ',s0_base
      write(3,*)'s0surf    ',s0surf
      write(3,*)'s0         ',s0
      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

      end subroutine initial_state_eq

! comparaisons variables Z:
!     IF(NBDOM_MECO.EQ.2) THEN
!     OPEN(UNIT=3,FILE='h'//DOM_C//'.out')
!     I1=0
!     IF(DOM_C.EQ.'10')I1=193
!     DO K=1,NR-1
!     DO I=-1,MECO+2
!     DO J=-1,NECO+2
!     WRITE(3,*)I+I1,J,K
!    &  ,SAL_Z(I,J,K),TEM_Z(I,J,K)
!     ENDDO
!     ENDDO
!     ENDDO
!     CLOSE(3)
!     ELSE
!     OPEN(UNIT=3,FILE='h1.out')
!     DO K=1,NR-1
!     DO I=-1,197
!     DO J=-1,NECO+2
!     WRITE(3,*)I,J,K
!    &  ,SAL_Z(I,J,K),TEM_Z(I,J,K)
!     ENDDO
!     ENDDO
!     ENDDO
!     CLOSE(3)
!     OPEN(UNIT=3,FILE='h2.out')
!     DO K=1,NR-1
!     DO I=192,MECO+2
!     DO J=-1,NECO+2
!     WRITE(3,*)I,J,K
!    &  ,SAL_Z(I,J,K),TEM_Z(I,J,K)
!     ENDDO
!     ENDDO
!     ENDDO
!     CLOSE(3)
!     ENDIF
!     CALL MPI_BARRIER(PAR%COMM2D,K_OUT)      ! synchro processes
!     stop 'Alors mon cher Alfred?'
