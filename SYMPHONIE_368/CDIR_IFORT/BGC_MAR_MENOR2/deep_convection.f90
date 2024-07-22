










      subroutine deep_convection
!...............................................................................
! SYMPHONIE ocean model
! release 253 - last update: 05-05-19
!...............................................................................
      use module_principal ; use module_parallele
      implicit none
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
!...............................................................................
!...............................................................................
! Version date      Description des modifications
! v253    05-05-19  kmerged_t remplacE par kmerged_t
!...............................................................................
! https://docs.google.com/presentation/d/19f3SSwAQ4flJZG0lJAP-axdwOeSqtoRbucJojyTE_Ac/edit#slide=id.g57193dca64_0_11

      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax+1)==1) then !mskmskmskmsk>

       if(rhp_t(i,j,kmax)>rhp_t(i,j,kmax-1))then !instcond.instcond.instcond>

         rhp_t(i,j,kmax+1)  =rhp_t(i,j,kmax)   ! memoriser rhp_surf
         tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2) ! memoriser Tsurf
         sal_t(i,j,kmax+1,2)=sal_t(i,j,kmax,2) ! memoriser Ssurf

         do k=kmax-1,kmerged_t(i,j)+1,-1

! Deplacer la colonne vers le haut pour combler le vide du niveau de surface:
! T(k+1)=(T(k+1)*(dz(k+1)-dz(kmax))+T(k)*dz(kmax))/dz(k+1)
         x2=dz_t(i,j,kmax,2)/dz_t(i,j,k+1,2)
         x1=1.-x2
         tem_t(i,j,k+1,2)=x1*tem_t(i,j,k+1,2)+x2*tem_t(i,j,k,2)
         sal_t(i,j,k+1,2)=x1*sal_t(i,j,k+1,2)+x2*sal_t(i,j,k,2)
         rhp_t(i,j,k+1  )=x1*rhp_t(i,j,k+1  )+x2*rhp_t(i,j,k  )

          if(rhp_t(i,j,kmax+1)<rhp_t(i,j,k-1))then              ! stablev.stablev.stablev>

! Melanger l'eau de surface avec l'eau du niveau d'equilibre statique:
! T(k)=(T(k)*(dz(k)-dz(kmax))+T(kmax)*dz(kmax))/dz(k)
          x2=dz_t(i,j,kmax,2)/dz_t(i,j,k,2)
          x1=1.-x2
          tem_t(i,j,k,2)=x1*tem_t(i,j,k,2)+x2*tem_t(i,j,kmax+1,2)
          sal_t(i,j,k,2)=x1*sal_t(i,j,k,2)+x2*sal_t(i,j,kmax+1,2)
          rhp_t(i,j,k  )=x1*rhp_t(i,j,k  )+x2*rhp_t(i,j,kmax+1)

! le rearrangement s'arrete là, sortir de la boucle:
           goto 1000


          endif                                         ! stablev.stablev.stablev>


         enddo

! Cas particulier où l'eau de surface est plus dense que celle du fond:
         k=kmerged_t(i,j)

         x2=dz_t(i,j,kmax,2)/dz_t(i,j,k+1,2)
         x1=1.-x2
         tem_t(i,j,k+1,2)=x1*tem_t(i,j,k+1,2)+x2*tem_t(i,j,k,2)
         sal_t(i,j,k+1,2)=x1*sal_t(i,j,k+1,2)+x2*sal_t(i,j,k,2)
         rhp_t(i,j,k+1  )=x1*rhp_t(i,j,k+1  )+x2*rhp_t(i,j,k  )

         x2=dz_t(i,j,kmax,2)/dz_t(i,j,k,2)
         x1=1.-x2
         tem_t(i,j,k,2)=x1*tem_t(i,j,k,2)+x2*tem_t(i,j,kmax+1,2)
         sal_t(i,j,k,2)=x1*sal_t(i,j,k,2)+x2*sal_t(i,j,kmax+1,2)
         rhp_t(i,j,k  )=x1*rhp_t(i,j,k  )+x2*rhp_t(i,j,kmax+1)

 1000  endif                                     !instcondi.instcond.instcond>

      endif                          !mskmskmskmsk>
      enddo
      enddo

      if(iteration3d/=0)return

! A l'initialisation vérifier que la maille de surface est la plus petite de sa colonne:
      flag_stop=0
      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmax+1)==1) then !>>>>>>>>>
       do k=kmerged_t(i,j),kmax-1
        if(dsig_t(i,j,k)<0.999*dsig_t(i,j,kmax)) then !m°v°m>
         flag_stop=1
         write(10+par%rank,*)'L''algo de convection profonde impose que' &
        ,' la couche de surface soit toujours la plus fine, '            &
        ,' or cette condition n''est pas rempli en i j (glob) k='        &
        ,i+par%timax(1),j+par%tjmax(1),k
         write(10+par%rank,*)'dsig(kmax)',dsig_t(i,j,kmax)
         write(10+par%rank,*)'dsig(k   )',dsig_t(i,j,k)
         write(10+par%rank,*)'h_w(i,j)'  ,h_w(i,j)
        endif                                  !m°v°m>
      enddo

      endif                          !>>>>>>>>>
      enddo
      enddo
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'Err subroutine deep_convection. See fortxxx files'

      end subroutine deep_convection
