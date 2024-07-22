










      module module_subcycle
!______________________________________________________________________
! S model
! release S26.1 - last update: 10-11-13
!______________________________________________________________________
! Version date     Description des modifications
! S.26    10-11-13 correction bug fontions mpi_allreduce
!...............................................................................
      use module_principal ; use module_parallele
      implicit none

contains

!..........................................................................

      subroutine subcycle_exchange_status(iteration_)
      implicit none
      integer iteration_


      subcycle_exchange=8
      mpi_neighbor_list(1:8)=(/ ouest, est, nord, sud, sudouest, sudest, nordouest, nordest /)

      if(par%tvoisin(sud)/=mpi_proc_null) then !----------->
       if(dte_lp<glob_dte_lp(par%tvoisin(sud))) then !******>
        if(mod(iteration_,2)==0) then                 !>>>>>>>>>>
         subcycle_exchange=5
         mpi_neighbor_list(1:5)=(/ ouest, est, nord, nordouest, nordest /)
        endif                                         !>>>>>>>>>>
       endif                                         !******>
      endif                                    !----------->
      end subroutine subcycle_exchange_status

!..........................................................................

      subroutine subcycle_synchro_status(iteration_)
      implicit none
      integer iteration_
! Certaines actions necessitent que tous les rank soient en action
! telles que l'ecriture d'un fichier netcf ou l'arret du code.
! La variable subcycle_synchro=1 dans ce cas.
! par%subcycle indique le plus grand ratio de pas de temps parmi les proc

      if(mod(iteration_,par%subcycle)==0) then

       subcycle_synchro=1

      else

       subcycle_synchro=0

      endif
      end subroutine subcycle_synchro_status

!..........................................................................

      subroutine subcycle_dte(dteglob_)
      implicit none
      double precision dteglob_         & ! minimum barotropic time step for the whole domain
                      ,dte_lp_prv_      & ! valeur de dte_lp avant recalcul
                      ,instability_loc_ &
                      ,instability_glb_

      allocate(glob_dte_lp(0:nbdom-1))     ; glob_dte_lp=0.d00
      allocate(glob_dte_lp_tmp(0:nbdom-1)) ; glob_dte_lp_tmp=0.d00

      dte_lp_prv_=dte_lp

! Dispatch the biggest par%subcycle (number of cycles in one principal time step)
! Gives the subcycle modulo relative to the principal cycle
      call mpi_allreduce(par%subcycle                                 &
                        ,subcycle_modulo                              &
                        ,1                                            & !10-11-13
                        ,mpi_integer,mpi_max,par%comm2d,ierr)
! Compute the barotropic time step of each subdomain assuming that dteglob_ (the minimum
! barotropic time step imposed by the CFL condition) corresponds to the biggest par%subcycle
      dte_lp=dteglob_*subcycle_modulo/par%subcycle

! If instability_>1, the new time step does not respect the CFL criteria
! Decrease the minimum time step using instability_glb_ as a scale factor restores the CFL condition
! for all subdomains
      instability_loc_=max(dte_lp/dte_lp_prv_,1.d0)
      call mpi_allreduce(instability_loc_                             &
                        ,instability_glb_                             &
                        ,1                                            &
                        ,mpi_double_precision,mpi_max,par%comm2d,ierr)  !10-11-13
      dte_lp=dte_lp/instability_glb_

!     write(6,*)'HELLO PAR%RANK DTE_LP APRES',par%rank,dte_lp
!     write(6,*)'HELLO PAR%RANK ',par%rank,par%subcycle,subcycle_modulo
!     write(6,*)'HELLO',par%rank,instability_loc_,instability_glb_
!     if(par%rank<4)write(6,*)par%rank,real(dte_lp_prv_)              &
!                       ,real(dte_lp*instability_glb_)  &
!                       ,real(dte_lp)

!Dispatch dte_lp to glob_dte_lp
      glob_dte_lp_tmp(par%rank)=dte_lp
      call mpi_allreduce(glob_dte_lp_tmp                              &
                        ,glob_dte_lp                                  &
                        ,par%nbdom                                    &
                        ,mpi_double_precision,mpi_sum,par%comm2d,ierr)


      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!     stop 'mimi'

      deallocate(glob_dte_lp_tmp)

      end subroutine subcycle_dte


      end module module_subcycle
