      subroutine AtmDeposition
      use module_principal
      use module_biology
      use ModuleDeclaration
      use module_parallele
      use module_dust
      use module_depnit 
      use module_depammo 
      implicit none
      double precision :: convP,convN


! Update of the deposition
      if  (idust==1)  call    dust_driver(2) 
      if(idepnit==1)  call  depnit_driver(2)
      if(idepammo==1) call depammo_driver(2)
      

!      if(par%rank==1) then
!       print*,'idepammo=',idepammo
!       stop'dans AtmDeposition'
!      endif


      convN=1./24/3600/14.
      if(ichoicedust.eq.2) then       
! pour regcm
      convP=1./24/3600/31. ! on convertit de mg/m2/d en mmol/m2/s
      elseif(ichoicedust.eq.1.or.ichoicedust.eq.3) then
! pour aladin
       convP=1000000/31. ! on convertit de kg/m2/s en mmol/m2/s
      endif 
!      convN=1./24/3600/14.
!      print*,'conv=',conv
!      stop'verif'

      do i=1,imax
      do j=1,jmax

      if(idepnit==1) then
      AtmDepNit(i,j)=(ddnit_w(i,j,1)+wdwnit_w(i,j,1)+wdrnit_w(i,j,1))*convN
       fluxbio_w(i,j,initrate,2) = AtmDepNit(i,j)
      endif

      if(idepammo==1) then
!     if(par%rank==0) print*,'passe Dep Ammo'
      AtmDepAmmo(i,j)=(ddammo_w(i,j,1)+wdwammo_w(i,j,1)+wdrammo_w(i,j,1))*convN
       fluxbio_w(i,j,iammonium,2) = AtmDepAmmo(i,j)
      endif

      if  (idust==1) then
       if(ichoicedust.le.2) then      
      AtmDepPho(i,j)=(dd_w(i,j,1)+wdw_w(i,j,1)+wdr_w(i,j,1))*psolublefr*dustpfr*convP
       elseif(ichoicedust.eq.3) then
      AtmDepPho(i,j)=(dd_w(i,j,1)+wdr_w(i,j,1))*psolublefr*dustpfr*convP
       endif
       fluxbio_w(i,j,iphosphate,2) = AtmDepPho(i,j)
!      print*,'passe AtmDep',dd_w(i,j,1),wdw_w(i,j,1),wdr_w(i,j,1)
!      print*,'passe AtmDep2',AtmDepPho(i,j),psolublefr,dustpfr,conv 
      endif

      enddo
      enddo

      return
      end subroutine AtmDeposition      
