
! LIGNES 243-291
!***************************************************************************
! CONDITIONS LIMITES: FLEUVES
! DEBUT
      if(ichoix==1) then
!***************************************************************************

!-----------------------------------------------------------------------
! Modifs pour Tung et Nadia (debut)
! HERE ADD A MODIFICATION OF THE CONCENTRATION AT THE SOURCE POINT OF THE RIVER
! time1 = Start time in seconds and time2 = End time in seconds
      time1=0.  
      do vb=1,vbmax      ! Tracer number loop
       time2=31536000.   
       if(elapsedtime_now>time1.and.elapsedtime_now<time2)then
          river_bio(vb,1:nriver,1)=100.
       else
          river_bio(vb,1:nriver,1)=0.
       endif
       time1=time2
       time2=time2+31536000.
      enddo
! Modifs pour Tung et Nadia (fin)
!-----------------------------------------------------------------------

      do kr=1,nriver

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
!     if(river_dom(kr) == par%rank) then !rrrrrrrrrrrrrrrrrr>
      if(rivertrc_inout(kr)==1)        then !rrrrrrrrrrrrrrrrrr>          !01-11-10
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      i=iriver(kr,1)
      j=jriver(kr,1)
      if(riverdir(kr)/=0) then !------>
          k1=kmin_w(i,j) ; k2=kmax
      else                     !------>
               k1=kmax+1 ; k2=kmax+1    !07-11-12
      endif                    !------>

      do vb=1,vbmax

       do 90 k=k1,k2 ! kmin_w(i,j),kmax !07-11-12

      if(flag_nemoffline/=1) then         ! 25-12-14
       bio_t(i,j,k,vb)=river_bio(vb,kr,1)                          !04/04/06
      else
       bio_t(i,j,k,vb)=0.
      endif

   90  continue

       enddo

!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!05-06-09
! Si la riviere est dans le sous-domaine alors appliquer la C.L.
      endif                              !rrrrrrrrrrrrrrrrrr>
!#MPI Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       enddo

!***************************************************************************
! CONDITIONS LIMITES: FLEUVES
! FIN
      return
      endif
!***************************************************************************
