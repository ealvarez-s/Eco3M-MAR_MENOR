










      subroutine io_switch(txt_)
!______________________________________________________________________
! S model
! release S26 - last update: 15-07-15
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      implicit none
      character txt_*1,txt_spy_*60
!______________________________________________________________________________
! Version date      Description des modifications
! 2010.10 19-06-10  ecrire ou lire l'interrupteur
! 2010.23 21-05-11  Appel à routine fortran barriere.F90
!         27-05-11  variable texte txt_spy_ permet de savoir qui
!                   appelle la barriere
! 2010.25 09-03-12  Seul proc0 lit interrupteur
! S26     09-03-13  Ne pas passer iteration3d en argument pour subcycling
!         28-06-13  subcycle
!         17-01-15  les arguments de la fonction mpi_bcast sont mpi_integer
!         15-07-15  interrupteur ecrit dans repertoire tmpdirname
!______________________________________________________________________________

      if(subcycle_synchro==0)return                                  !28-06-13
      txt_spy_='io_switch'

!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes
      if(txt_=='r')k2=-1
      if(txt_=='w')k2=-2
      call barriere(60,k2,txt_spy_) !27-05-11

      if(txt_=='r') then !rrrrrrrrrrrrrrrrr>
       if(par%rank==0) then !#mpi>>>>>>>>>>>>>
        open(unit=3,file=trim(tmpdirname)//'interrupteur') !15-07-15
        read(3,'(i1,9x,i1,9x,i1)')kstop,kpvwave,give_chanel9
        close(3)
       endif                !#mpi>>>>>>>>>>>>>
      call mpi_bcast(kstop       ,1,mpi_integer,0,par%comm2d,ierr) !17-01-15
      call mpi_bcast(kpvwave     ,1,mpi_integer,0,par%comm2d,ierr) !17-01-15
      call mpi_bcast(give_chanel9,1,mpi_integer,0,par%comm2d,ierr) !17-01-15
      endif                 !rrrrrrrrrrrrrrrrr>

      if(txt_=='w') then !wwwwwwwwwwwwwwwww>
       if(par%rank==0) then !#mpi>>>>>>>>>>>>>
        open(unit=3,file=trim(tmpdirname)//'interrupteur') !15-07-15
        write(3,'(i1,9x,i1,9x,i1)')kstop,kpvwave,give_chanel9
        close(3)
       endif                !#mpi>>>>>>>>>>>>>
      endif                 !wwwwwwwwwwwwwwwww>

!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes
      if(txt_=='r')k2=-3
      if(txt_=='w')k2=-4
      call barriere(60,k2,txt_spy_) !27-05-11

      end subroutine io_switch
