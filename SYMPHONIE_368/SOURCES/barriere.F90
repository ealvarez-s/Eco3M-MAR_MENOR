      subroutine barriere(code1_loc,code2_loc,txt_loc)
!______________________________________________________________________
!
! S model
! release 2010.25  - last update: 26-03-12
! Contact: sirocco@aero.obs-mip.fr
! Laboratoire d'Aerologie, 14 Avenue Edouard Belin, F-31400 Toulouse
! http://
!
!______________________________________________________________________
      use module_principal
      use module_parallele !#mpi
      implicit none
      integer code1_loc,code2_loc,error_loc,k_loc
!     character txt_loc*60
      character(len=*) txt_loc
!______________________________________________________________________
! Version date     Description des modifications
! 2010.23 21-05-11 Mise en service
!         23-05-11 Etiquette de compilation #ifdef parallele
!         27-05-11 Indiquer à l'écran quelle routine a appelé la subroutine
!                  barriere
! 2010.25 29-01-12 k remplce par k_loc
!         26-03-12 modif declaration txt_loc
!______________________________________________________________________
#ifdef synopsis
       subroutinetitle='barriere'
       subroutinedescription= &
      'checks that all ranks are well synchronised'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(nbdom==1)return

#ifdef parallele

      tab1_code(par%rank)=code1_loc
      tab2_code(par%rank)=code2_loc

      call mpi_allreduce(tab1_code,tab_code_glb,nbdom,mpi_integer,  & !#MPI
                         mpi_sum,par%comm2d,ierr)

      error_loc=0
      do k_loc=0,nbdom-1
       if(tab1_code(par%rank)/=tab_code_glb(k_loc)) then
        error_loc=1
        write(*,*)'erreur barriere code1 '     &
         ,par%rank,tab1_code(par%rank)         &
         ,k_loc,tab_code_glb(k_loc)
       endif
      enddo
      if(error_loc==1)then
       write(*,*)'Barriere appelée depuis subroutine:'   &     !27-05-11
                 ,trim(txt_loc)
       stop 'STOP in subroutine barriere.F90'
      endif


      call mpi_allreduce(tab2_code,tab_code_glb,nbdom,mpi_integer,  & !#MPI
                         mpi_sum,par%comm2d,ierr)


      error_loc=0
      do k_loc=0,nbdom-1
       if(tab2_code(par%rank)/=tab_code_glb(k_loc)) then
        error_loc=1
        write(*,*)'erreur barriere code2 '     &
          ,par%rank,tab2_code(par%rank),k_loc  &
          ,tab_code_glb(k_loc)
       endif
      enddo
      if(error_loc==1)then
       write(*,*)'Barriere appelée depuis subroutine:'   &     !27-05-11
                 ,trim(txt_loc)
       stop 'STOP in subroutine barriere.F90'
      endif

#endif

      end subroutine barriere
