










      subroutine check_mpi
!______________________________________________________________________
! S model
! release S26 - last update: 27-11-27
! S26     27-11-15   ajout check_mpi_1dv_obc check_mpi_1dv_ts 
!______________________________________________________________________
      implicit none
      end subroutine check_mpi
!........................................................................

      subroutine check_mpi_1dv_internal
      use module_principal
      use module_parallele
      implicit none

      i1=imax/2 ; j1=jmax/2

      do k0=-1,2
      do k=1,kmax

      do j=1,jmax
      do i=1,imax+1
       if(vel_u(i,j,k,k0)/=vel_u(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour vel_u en i j k t=',i,j,k,k0
        write(6,*)vel_u(i,j,k,k0),vel_u(i1,j1,k,k0)
        stop 'check_mpi_1dv_internal'
       endif
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       if(vel_v(i,j,k,k0)/=vel_v(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour vel_v en i j k t=',i,j,k,k0
        write(6,*)vel_v(i,j,k,k0),vel_v(i1,j1,k,k0)
        stop 'check_mpi_1dv_internal'
       endif
      enddo
      enddo

      enddo ! k
      enddo ! k0


      end subroutine check_mpi_1dv_internal

!..........................................................................

      subroutine check_mpi_1dv_external
      use module_principal
      use module_parallele
      implicit none

      i1=imax/2 ; j1=jmax/2

      do k0=-1,2

      do j=1,jmax
      do i=1,imax+1
       if(velbar_u(i,j,k0)/=velbar_u(i1,j1,k0)) then
        write(6,*)'modele 1dv non regulier pour velbar_u en i j t=',i,j,k0
        write(6,*)velbar_u(i,j,k0),velbar_u(i1,j1,k0)
        stop 'check_mpi_1dv_external'
       endif
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       if(velbar_v(i,j,k0)/=velbar_v(i1,j1,k0)) then
        write(6,*)'modele 1dv non regulier pour velbar_v en i j t=',i,j,k0
        write(6,*)velbar_v(i,j,k0),velbar_v(i1,j1,k0)
        stop 'check_mpi_1dv_external'
       endif
      enddo
      enddo

      enddo ! k0

      do k0= 1,2
      do j=1,jmax
      do i=1,imax
       if(ssh_w(i,j,k0)/=ssh_w(i1,j1,k0)) then
        write(6,*)'modele 1dv non regulier pour ssh_w en i j t=',i,j,k0
        write(6,*)ssh_w(i,j,k0),ssh_w(i1,j1,k0)
        stop 'check_mpi_1dv_external'
       endif
      enddo
      enddo
      enddo

      end subroutine check_mpi_1dv_external

!........................................................................

      subroutine check_mpi_1dv_obc !27-11-15
      use module_principal
      use module_parallele
      implicit none

      i1=imax/2 ; j1=jmax/2

      do k0=0,2
      do k=1,kmax

      do j=1,jmax
      do i=1,imax
       if(temobc_t(i,j,k,k0)/=temobc_t(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour temobc en i j k t=',i,j,k,k0
        write(6,*)temobc_t(i,j,k,k0),temobc_t(i1,j1,k,k0)
        stop 'check_mpi_1dv_obc'
       endif
      enddo
      enddo

      do j=1,jmax
      do i=1,imax
       if(salobc_t(i,j,k,k0)/=salobc_t(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour salobc en i j k t=',i,j,k,k0
        write(6,*)salobc_t(i,j,k,k0),salobc_t(i1,j1,k,k0)
        stop 'check_mpi_1dv_obc'
       endif
      enddo
      enddo

      do j=1,jmax
      do i=1,imax+1
       if(velobc_u(i,j,k,k0)/=velobc_u(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour velobc_u en i j k t=',i,j,k,k0
        write(6,*)velobc_u(i,j,k,k0),velobc_u(i1,j1,k,k0)
        stop 'check_mpi_1dv_obc'
       endif
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       if(velobc_v(i,j,k,k0)/=velobc_v(i1,j1,k,k0)) then
        write(6,*)'modele 1dv non regulier pour velobc_v en i j k t=',i,j,k,k0
        write(6,*)velobc_v(i,j,k,k0),velobc_v(i1,j1,k,k0)
        stop 'check_mpi_1dv_obc'
       endif
      enddo
      enddo

      enddo ! k
      enddo ! k0


      i1=imax/2 ; j1=jmax/2

      do k0=0,2

      do j=1,jmax
      do i=1,imax+1
       if(velbarobc_u(i,j,k0)/=velbarobc_u(i1,j1,k0)) then
        write(6,*)'modele 1dv non regulier pour velbarobc_u en i j t=',i,j,k0
        write(6,*)velbarobc_u(i,j,k0),velbarobc_u(i1,j1,k0)
        stop 'check_mpi_1dv_external'
       endif
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
       if(velbarobc_v(i,j,k0)/=velbarobc_v(i1,j1,k0)) then
        write(6,*)'modele 1dv non regulier pour velbarobc_v en i j t=',i,j,k0
        write(6,*)velbarobc_v(i,j,k0),velbarobc_v(i1,j1,k0)
        stop 'check_mpi_1dv_external'
       endif
      enddo
      enddo

      enddo ! k0


      end subroutine check_mpi_1dv_obc

!..........................................................................

      subroutine check_mpi_1dv_ts(txt_)
      use module_principal
      use module_parallele
      implicit none
      character*1 txt_

      i1=imax/2 ; j1=jmax/2

      do k0=-1,2
      do k=1,kmax

      do j=1,jmax
      do i=1,imax

       if(tem_t(i,j,k,k0)/=tem_t(i1,j1,k,k0)) then
        write(6,'(a,a)')'modele 1dv non regulier pour tem_t etape:',txt_
        write(6,*)'i j k t=',i,j,k,k0
        write(6,*)tem_t(i,j,k,k0),tem_t(i1,j1,k,k0)
        stop 'check_mpi_1dv_ts'
       endif

       if(sal_t(i,j,k,k0)/=sal_t(i1,j1,k,k0)) then
        write(6,'(a,a)')'modele 1dv non regulier pour sal_t etape:',txt_
        write(6,*)'i j k t=',i,j,k,k0
        write(6,*)sal_t(i,j,k,k0),sal_t(i1,j1,k,k0)
        stop 'check_mpi_1dv_ts'
       endif

      enddo
      enddo

      enddo ! k
      enddo ! k0


      end subroutine check_mpi_1dv_ts

!..........................................................................

      subroutine check_mpi_1dv_airsea
      use module_principal
      use module_parallele
      implicit none

      i1=imax/2 ; j1=jmax/2
      do j=1,jmax
      do i=1,imax

       if(slhf_w(i,j,1)/=slhf_w(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier pour slhf_w en i j =',i,j
        write(6,*)slhf_w(i,j,1),slhf_w(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif

       if(snsf_w(i,j,1)/=snsf_w(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier pour snsf_w en i j =',i,j
        write(6,*)snsf_w(i,j,1),snsf_w(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif

       if(sshf_w(i,j,1)/=sshf_w(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier pour sshf_w en i j =',i,j
        write(6,*)sshf_w(i,j,1),sshf_w(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif

       if(precipi_w(i,j,1)/=precipi_w(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier pour precipi_w en i j =',i,j
        write(6,*)precipi_w(i,j,1),precipi_w(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif

      enddo
      enddo

      do j=2,jmax-1
      do i=2,imax
       if(wstress_u(i,j,1)/=wstress_u(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier wstress_u en i j =',i,j
        write(6,*)wstress_u(i,j,1),wstress_u(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
       if(wstress_v(i,j,1)/=wstress_v(i1,j1,1)) then
        write(6,*)'modele 1dv non regulier wstress_v en i j =',i,j
        write(6,*)wstress_v(i,j,1),wstress_v(i1,j1,1)
        stop 'check_mpi_1dv_airsea'
       endif
      enddo
      enddo

      end subroutine check_mpi_1dv_airsea

!........................................................................

      subroutine check_mpi_1dv_bio(txt_)
      use module_principal
      use module_parallele
      implicit none
      character*1 txt_

      i1=imax/2 ; j1=jmax/2

      do k0=1,vbmax

      do k=1,kmax
      do j=1,jmax
      do i=1,imax

       if(bio_t(i,j,k,k0)/=bio_t(i1,j1,k,k0)) then
        write(6,'(a,a)')'modele 1dv non regulier pour bio_t etape:',txt_
        write(6,*)'i j k vb=',i,j,k,k0
        write(6,*)bio_t(i,j,k,k0),bio_t(i1,j1,k,k0)
        stop 'check_mpi_1dv_bio'
       endif

      enddo
      enddo

      enddo ! k

      enddo ! k0

      end subroutine check_mpi_1dv_bio
