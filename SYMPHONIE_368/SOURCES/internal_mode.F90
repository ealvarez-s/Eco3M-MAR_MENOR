     subroutine internal_mode
!______________________________________________________________________
! SYMPHONIE ocean model
! release 347 - last update: 01-05-22
!______________________________________________________________________
      use module_principal ; use module_parallele ; use module_q
      use module_s         ; use module_webcanals ; use module_offline
      implicit none

!$ Driving routine for the internal mode velocity subroutines

!.....................................................................
! Version date      Description des modifications
!         02/07/06: l'appel à vishor est deplacé si la diffusion est
!                   appliquee aux vitesses aux temps t-1
!         22/07/06: L'appel a la  routine advection_vel disparait car
!                   l'advection est regroupee dans uv_upd.F
!         21/12/06: conditions aux limites dans les fleuves calculees
!                   dans obc_river. Ceci permet de rappeler la condition
!                   aux limites sur les rivieres quand des operations
!                   faisant intervenir le masque annule l'effet de la
!                   riviere comme dans couple_mode.F
! 2009.3  03-10-09  internal_mode.F90 devient internal_mode.F
!         04-10-09  - uv_upd.F devient momentum_equations.F
!                   - omega_upd.F devient omega.F
!         08-10-09  advection verticale lagrangienne: appel à
!                   cellbox_thickness
!         19-10-09  des arguments passés dans routine cellbox_thickness
!         16-11-09  application de la condition aux limites des fleuves
!                   aux vitesses du schema "forward"
! 2010.7  15-02-10  momentum_equation compte 2 parties: avant et apres
!                   le mode externe
! 2010.12 27-09-10  suppression arg "key" dans appel à routine omega
! 2010.14 16-12-10  routines de turbulence passées de model_3d.F90 dans
!                   internam_mode.F90
! 2010.22 22-04-11  Controle lagrangien de la pression de fond en preparation
! S26.1   04-09-12  Procedure sans time splitting
!         14-05-13  Corriger vel(2) avec dernier instantanee velbar pour continuite
!                   temporelle du terme de la force de stokes dans module_external
!         18-01-14  Possibilite (si activee) d'appeler routine couple_modes_3v qui
!                   produit ub flux advectif base sur 3 vitesses en equilibre
!                   avec la variation de ssh entre t-1 et t+1
!         08-05-14  ajout obc_int_mpi (nouveaux echanges)
!         05-11-14  retour a version anterieure a aout 2014
!         26-02-15  call obc_int_mpi deplace avant call obc_river(2,2,2)
!         15-07-15  modification des arguments passes dans obc_int_mpi
!         01-10-15  ajout calcul cpu_seconds
!         28-11-15  suppression i1d
!         04-01-16  appel A call time_step_dtimax 
!         19-01-16  possibilite zonal restoring pour cas test jet barocline
!         28-01-16  nouveau time stepping LF-FB correction de vel(2) avec 
!                   instantanE
!         03-02-16  Ne plus faire l'echange "u2" sur vel(1) dans couple_mode
!                   impose de faire echange "u12" sur vel(2) dans internal_mode
!         28-03-16  Correction de vel(t+1) avec la valeur instantannEe
!                   de velbar dans la cas du schema FB - implications mpi
!         26-04-16  ajout call turbulence_driver  !26-04-16
!         10-08-16  ajout appels vertmix_coarser2_xxx pour le moment commentEs
!         04-09-16  Pour eviter singularite wetdry le LF ne corrige plus vel(t) 
!                   en meme temps que veldydz, en contre partie de quoi
!                   on corrige vel(t+1) avec la valeur instantanee de velbar
!         29-09-16  correction du courant moyen instantannE de vel(2) en
!                   prenant compte d'un profil dans la couche de fond
!         08-04-17  s-z coordinate
!         18-04-17  Apres modification de obc_river, suppression de call obc_river(2,2,2)
!         17-10-17  couche fusionnee: correction barotrope homogene
!                   reatblie (en attendant mieux) et suppression routine couche ecrasee
!         23-10-17  ajout subroutine internal_mode_nosplitting 
!         11-10-17  modifs sur la chronologie des variables en mode sans time-splitting
!         17-08-18  cas flag_nh3d_uv > 1 pour T,S NH
!         03-09-17  call couple_modes_scalars dans cette routine avant turbulence dont
!                   l'advection utilise la mise A jour (par couple_modes_scalars) des
!                   champs advectant
!         06-09-18  si flag_timesplitting_adve_uv/=0 uppression du "couple_mode et de omega des vitesses"
!                   Pour le cas NH 3D le time-stepping deplace le calcul de tken apres les vitesses
!                   et du coup couple_mode_scalar est calculE juste avant la turbulence
!         09-10-18  if(flag_nh3d==1)call q_wavebreaking3d 
!         12-01-18  if(iteration3d>nh2d_graph_spinup)call graph_out 
!         17-10-18  if(flag_nh3d==1.or.flag_nh3d==2)call q_wavebreaking3d 
!         19-01-19  call wetdry_mask_nosplitting_uv 
! v252    18-04-19  Les canaux envoient u,t,s,ssh A la mer
! v256    06-06-19  if(ioffline==1)call offline_inout(1) est placE dans internal_mode juste aprEs couple_mode 
!                   et omega. Pourquoi? Parce qu'avec le wetdrying no6 les flux horizontaux peuvent Etre ralentis avec 
!                   comme consequence qu'ils ne sont plus egaux aux fluxbar_sum. Il faut donc archiver dans les fichiers
!                   offline les vitesses 3D et omega_surf dans la partie "momentum" pour avoir la coherence entre fluxbar_sum
!                   et les autres champs dynamiques.
! v259    01-10-19  call obc_ssh_driver
!                   call wave_ww3_spec_sete_driver
! v280    02-05-20  call wave_hamilton_ebersole   
! v347    01-05-22  call q_wavebreaking3d est appelE avant obc_int
!...............................................................................
!    _________                    .__                  .__                     !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      itime=0
      if(iteration2d_max_now==1) stop ' Err 74 iteration2d_max_now==1 ' ! Ne pas passer pas internal_mode si iteration2d_max_now=1

! NOTER QUE LE CAS STANDARD EST flag_timesplitting_adve_uv=0
! the external_mode current is substituted to the depth averaged 3D current:
! Note: ne se fait pas si l'advection de u,v se fait avec les champs advectant de l'etape advection des scalaires
! NOTER QUE LE CAS STANDARD EST flag_timesplitting_adve_uv=0
      if(flag_timesplitting_adve_uv==0)call couple_modes !06-09-18

! computes velocities(now) at river mouth (masked) points:
      call obc_river(2,1,1)  ! "leap frog" and "forward" velocities 

! computes the thickness of cell boxes:
      call cellbox_thickness(1,0,0) ; call s_cpu('cellbox_thickness',0)


! computes the omega velocity (for the eulerian part of vertical advection):
! Note: ne se fait pas si l'advection de u,v se fait avec les champs advectant de l'etape advection des scalaires
      if(flag_timesplitting_adve_uv==0)call omega !06-09-18

! time-averaging for post-processing files
      if(ioffline==1)call offline_inout(1) !06-06-19 ! voir notes dans entEte expliquant le deplacement ici

!     if(flag_timesplitting_adve_uv==0) then ! pas obligatoire de le maintenir ici mais comme ca on reste identique A la version precedente
!         call turbulence_driver  ; call s_cpu('turbulence',0) !03-09-18
!     endif

! computes the 3D momentum equations:
      call momentum_equations('after external mode ')              ; call s_cpu('momentum_equations',0)

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
!#endif


! Avec le schema FB on corrige la moyenne de vel(t+1) avec velbar
! instantannee. L'option est differente en Leap-Frog "POM like" 
! oU (voir couple_mode.F90) on corrige vel(t) avec la moyenne temporelle
! du mode externe

      if(flag_timesplitting==1.and.flag_1dv==0)  then !23-10-17
!        call couple_modes_zprofil_instant !29-09-16
         call couple_modes_instant(2,2) !17-10-17
      endif                          

! Zonal restoring (baroclinic jet case)
       if(relaxtype_ts==40)call momentum_zonal_restoring !19-01-16

! Les extremitEs des canaux envoient vel_u(2)
       if(nbcanal>0)call webcanals_gt1_to_gt2_vel

! open boundary conditions:
      call obc_int
! mpi boundary conditions: !08-05-14 !26-02-15
! Echanges au temps t_=2 de type=1  (u1 et v1)
      if(timestep_type==timestep_leapfrog)call obc_int_mpi(2,1)  !28-03-16
      if(timestep_type==timestep_forwbckw)call obc_int_mpi(2,12) ! (t_,obctype_) avec: !03-02-16
      call s_cpu('obc_int_all',0)
                             ! t_=1 ou 2
                             ! obctype=1  pour (u1,v1) 
                             ! obctype=2  pour (u2,v2) 
                             ! obctype=12 pour (u1,v1) & (u2,v2) 



!     if(flag_nh3d/=0) then ! nh pressure>
!       stop 'passe pas par breaking'
!       call q_wavebreaking3d
! Cette ligne est commentee depuis le choix d'ordre chrono (3D uniquement)
! suivant: q(0)  u(0)  q(1)  u(1)  q(2) 
! et deplacee dans internal_mode_nosplitting
!!!!    call q_update_nhp
! details dans https://docs.google.com/document/d/1fDW9hoFog-DtjuAVMhkiqz5FhrIqKdC5JHDJqy-FgFs/edit
!       call s_cpu('nh pressure & wavebreaking',0)
!       if(iteration3d>nh2d_graph_spinup) then !pmx>
!        if(mod(iteration3d,nh2d_graph_period)==0)call graph_out
!       endif                                  !pmx>
!     endif                 ! nh pressure>


! compute fluxes and omega
! Note: la routine suivante contient call omega
      call couple_modes_scalars !06-09-18
! note: update de veldxdz, omega aussitot vel(2) connue pour pouvoir
! continuer avec l'advection de la turbulence. Les termes de
! production par cisaillement de vitesse utiliseront vel(2)

! compute turbulence
!     if(flag_timesplitting_adve_uv==1) then !>>>
       call turbulence_driver  ; call s_cpu('turbulence',0) !06-09-18
!     endif                                  !>>>

! q(0) u(1) q(1)  u(2) ===> q(2)
        if(flag_nh3d/=flag_nh3d_none)call q_update_nhp !07-09-18
! Attention que le moveforward integree dans q_update_nhp fait qu'au sortir
! de cette routine nous avons q(1)=q(2)...

      end subroutine internal_mode

!..................................................................

      subroutine internal_mode_nosplitting !23-10-17
      use module_principal ; use module_parallele ; use module_q ; use module_external_mode
      use module_wave
      implicit none
! Details de l'alternace FB des variables
! https://docs.google.com/document/d/1fDW9hoFog-DtjuAVMhkiqz5FhrIqKdC5JHDJqy-FgFs/edit

      if(nriver/=0) stop 'Err 152 internal_mode_nosplitting nriver/=0'

!     call wave_ww3_spec_sete_driver !01-10-19
!     call wave_hamilton_ebersole    !02-05-20
!     call wave_hamilton_ebersole_jmax
!     call wave_ricore_jeq1      ! cas monochromatique
!     call wave_ricore_jeq1_freq ! cas plusieurs frequences


! omega, veldxdz_v, veldydz_u n'etant pas connus en iteration3d=0 on ne peut pas
! calculer correctement l'advection des moments et de la turbulence, en consequence
! de quoi la premiere iteration passe par dessus le calcul de u,v et tke:
      if(iteration3d==0)goto 157

! Turbulence
!     call turbulence_driver  !03-09-18

! computes wetting/drying wetmask_u & wetmask_v before momentum equations
      call wetdry_mask_nosplitting_uv !19-01-19

! computes the 3D momentum equations:
      call momentum_equations('after external mode ')         

! Wave breaking
      if(flag_nh3d==1.or. &
         flag_nh3d==2)call q_wavebreaking3d ! Avant obc_int svp !01-05-22

! open boundary conditions:
      call obc_int

! mpi boundary conditions: !08-05-14 !26-02-15
! Echanges au temps t_=2 de type=1  (u1 et v1)
      if(timestep_type==timestep_leapfrog)call obc_int_mpi(2,1)  !28-03-16
      if(timestep_type==timestep_forwbckw)call obc_int_mpi(2,12) ! (t_,obctype_) avec: !03-02-16
                             ! t_=1 ou 2
                             ! obctype=1  pour (u1,v1) 
                             ! obctype=2  pour (u2,v2) 
                             ! obctype=12 pour (u1,v1) & (u2,v2)

! Former les flux (moyenne pour les iterations t-1 et t)
  157 do j=1,jmax ; do i=1,imax+1
       fluxbar_u(i,j,1)=0.
      enddo       ; enddo
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

!      if(vel_u(i,j,k,after)/=0)stop 'coco'

         veldydz_u(i,j,k,1)=0.5*(        &
             vel_u(i,j,k,after)          &
             *dz_u(i,j,k,dafter)         & 
            +vel_u(i,j,k,now)            &
             *dz_u(i,j,k,dnow)   )       & 
             *dy_u(i,j)                  &
           *mask_u(i,j,k )

       fluxbar_u(i,j,1)= &
       fluxbar_u(i,j,1)+veldydz_u(i,j,k,1)

!      if(j/=1) then
!       if(veldydz_u(i,j,k,1)/=veldydz_u(i,j-1,k,1)) then
!        write(6,*)'iteration3d',iteration3d
!        write(6,*)'i,j',i,j
!        write(6,*)'pouette',veldydz_u(i,j,k,1),veldydz_u(i,j-1,k,1)
!        write(6,*)'pouette',vel_u(i,j,k,after),vel_u(i,j-1,k,after)
!        write(6,*)'pouette',vel_u(i,j,k,now),vel_u(i,j-1,k,now)
!        stop 'charles'
!       endif
!      endif 

      enddo ; enddo ; enddo
!     if(mod(iteration3d,100)==0)write(6,*)'bidouille pouette'

      do j=1,jmax+1 ; do i=1,imax
       fluxbar_v(i,j,1)=0.
      enddo       ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

         veldxdz_v(i,j,k,1)=0.5*(         &
             vel_v(i,j,k,after)           &
             *dz_v(i,j,k,dafter)          & 
            +vel_v(i,j,k,now)             &
             *dz_v(i,j,k,dnow)   )        & 
             *dx_v(i,j)                   & !16-10-15
           *mask_v(i,j,k )

       fluxbar_v(i,j,1)= &
       fluxbar_v(i,j,1)+veldxdz_v(i,j,k,1)

!      if(j/=1) then
!       if(veldxdz_v(i,j,k,1)/=veldxdz_v(i,j-1,k,1)) then
!        write(6,*)'iteration3d',iteration3d
!        write(6,*)'par%rank,i,j,k',par%rank,i,j,k
!        stop 'erreur vel_v'
!       endif
!      endif 

      enddo ; enddo ; enddo

! Conditions limites pour fluxbar (mpi et physique) 
!     call obc_ext_nosplitting
! note: la routine en suivant a ete essaye mais ne fonctionne pas
!     call obc_ext_nosplitting_clamped ! transport imposE

! Surface Elevation
      do j=1,jmax ; do i=1,imax

! attention modifier ce schema implique de changer aussi l'inversion de
! cette equation dans obc_ext_fb
         ssh_int_w(i,j,2)=                                        &
         max(0.001-h_w(i,j),                                      &
         ssh_int_w(i,j,0)                                         &
              -dti_lp*(                         & !dtdtdt>
                        ( fluxbar_u(i+1,j  ,1)                &
                         -fluxbar_u(i  ,j  ,1)                &
                         +fluxbar_v(i  ,j+1,1)                &
                         -fluxbar_v(i  ,j  ,1))/dxdy_t(i,j)   &

!                          +omega_w(i,j,kmax+1,1)             & !09-04-15
                      ))                          !dtdtdt>


!      if(i/=1) then
!       if(ssh_int_w(i,j,2)/=ssh_int_w(i-1,j,2)) then
!        write(6,*)'iter par%rank,i,j,2ssh',iteration3d,par%rank,i,j,ssh_int_w(i,j,2),ssh_int_w(i-1,j,2)
!        write(6,*)'du',fluxbar_u(i+1,j  ,1)-fluxbar_u(i  ,j  ,1) &
!                      ,fluxbar_u(i  ,j  ,1)-fluxbar_u(i-1,j  ,1)  
!        write(6,*)'dv',fluxbar_v(i  ,j+1,1)-fluxbar_v(i  ,j  ,1) &
!                      ,fluxbar_v(i-1,j+1,1)-fluxbar_v(i-1,j  ,1) 
!        stop 'fifi ssh'
!       endif
!      endif


!      if(iteration3d==2.and.par%rank==7.and.j==175)write(6,*)'i du',fluxbar_u(i+1,j  ,1)-fluxbar_u(i  ,j  ,1)

      enddo ; enddo

      call obc_ssh_flather !01-10-19

      if(isnan(ssh_int_w(imax/2,jmax/2,2))) then
       write(6,*)'iteration3d',iteration3d
       stop ' SSH_INT_W(imax/2,jmax/2,2) IS NAN'
      endif

#ifdef checkmpi
! Attention de bien mettre cet appel avant l'appel a obc_ssh qui echange desormais
! la zone de recouvrement pour les besoins du subcycling
      do j=1,jmax ; do i=1,imax
         ssh_w(i,j,2)=ssh_int_w(i,j,2)
      enddo ; enddo
      call external_mode_checkmpi !13-01-11
#endif

! Conditions limites mpi pour ssh_int_w
      call obc_ssh_int              ! this is not a physical boundary condition
      call obc_ssh_int_mpi('z1',2)

!  Thickness of the whole ocean layer:
      call z_thickness(2,2)

!  computes the wetdry mask (for _t & _w points only) !19-01-19
! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
      call wetdry_mask_airseafluxes

! computes the thickness of cell boxes:
      call cellbox_thickness(1,0,0) ; call s_cpu('cellbox_thickness',0)

! computes the omega velocity (for the eulerian part of vertical advection):
      call omega ; call s_cpu('omega',0)

! Turbulence 
      call turbulence_driver  !03-09-18
! note: l'advection de la turbulence utilise les champs mis A jour de l'equation de continuitE.
! CalculEe apres momentum, le terme de production par cisaillement de vitesse utilise vel(2)

! q(0) u(1) q(1)  u(2) ===> q(2)
        call q_update_nhp
! Attention que le moveforward integree dans q_update_nhp fait qu'au sortir
! de cette routine nous avons q(1)=q(2)...
! details dans https://docs.google.com/document/d/1fDW9hoFog-DtjuAVMhkiqz5FhrIqKdC5JHDJqy-FgFs/edit

!       if(mod(iteration3d,100)==0)write(6,*)'BIDOUILE Q=0'
!       q_t=0.


      if(mod(iteration3d,nh2d_graph_period)==0) then !m°v°m>
          if(iteration3d>nh2d_graph_spinup)call graph_out !12-10-18
      endif                                          !m°v°m>

#ifdef bidon
      if(mod(iteration3d,100)==0) then !>>>>>>>
! Deduire la vitesse de phase de la vitesse de deplacement de la
! valeur ssh=0
      sum0=small1
      sum1=0.
      j=jmax/2
      do i=2,imax
       if(ssh_int_w(i,j,2)*ssh_int_w(i-1,j,2)<0.and. &
          ssh_int_w(i,j,0)*ssh_int_w(i-1,j,0)<0.) then
        x0=ssh_int_w(i,j,0)/(ssh_int_w(i-1,j,0)-ssh_int_w(i,j,0))+real(i) ! indice du "zero" au temps 0
        x1=ssh_int_w(i,j,2)/(ssh_int_w(i-1,j,2)-ssh_int_w(i,j,2))+real(i) ! indice du "zero" au temps 0
        sum0=sum0+1.
        sum1=sum1+(x1-x0)*dxb/dti_lp
       endif
      enddo
      write(texte30,'(a,i0)')'tmp/nhdiag_phasespeed',par%rank
       open(unit=3,file=texte30,position='append')
        write(3,*)iteration3d*dti_fw/3600.,sum1/sum0
       close(3)
      endif                            !>>>>>>>
#endif

      end subroutine internal_mode_nosplitting

!..................................................................

!..................................................................
#ifdef bidon
      subroutine internal_mode_wave2dv_generator
      use module_principal
      use module_parallele
      implicit none
      real omega_,u0_,u1_,wave_vector_,dist_


      omega_=2.*pi/43200. ! periode 43200 secondes (onde semi-diurne)
      u0_=0.001*rampe
!     u1_=u0_*sin(omega_*elapsedtime_now)*rampe
      u1_=u0_*sin(omega_*elapsedtime_now)
      wave_vector_=omega_/c_wave_mode(1,1,1)  & ! (i=1,j=1,mode=1)
                  *sqrt(1.-(coriolis_t(1,1)/omega_)**2)
      i0=1 ; i=i0

      if(par%rank==0) then !000000000000000000>
       do k=1,kmax
       do j=1,jmax
        vel_u(i,j,k,2)=u1_*cos(pi*depth_u(i,j,k)/h_u(i,j))
       enddo
       enddo
       i2=2
       dist_=real(i2+par%timax(1)-i0)*dxb ! distance from the wave generation point
       do k=1,kmax
       do j=1,jmax+1
        vel_v(i2,j,k,2)=                                       &
           u0_*cos(pi*depth_v(i2,j,k)/h_v(i2,j))               &
              *cos(omega_*elapsedtime_now-dist_*wave_vector_)  &
              *coriolis_t(1,1)/omega_
       enddo
       enddo
      endif                !000000000000000000>

!     if(int(elapsedtime_now/86400.)-int(elapsedtime_bef/86400.)==1)then!----->
!     if(mod(iteration3d,100)==0) then          !---------------------->
      if(elapsedtime_aft>=elapsedtime_end) then !---------------------->
       if(par%rank==0) then !!!!!!!!>
         write(6,*)'Ecriture apres (jours):',elapsedtime_now/86400.
         write(6,*)'Temps (jours) mode1 traverse domaine='  &
        ,(iglb-2)*dxb/c_wave_mode(1,1,1)                    &
         *sqrt(1.-(coriolis_t(1,1)/omega_)**2) / 86400.
       endif                !!!!!!!!>
       k=kmax
       j=2
!      open(unit=3,file='profil_u.dat')
       open(unit=3,file=dom_c//'.dat')

         do i=3,imax ! commence 3 pour ne pas repeter mpi
          dist_=real(i+par%timax(1)-i0)*dxb ! distance from the wave generation point
          write(3,'(5(e13.6,1x))')                             &
           dist_                                               &
          ,vel_u(i,j,k,2)                                      &
          ,u0_*cos(pi*depth_u(i,j,k)/h_u(i,j))                 &
              *sin(omega_*elapsedtime_now-dist_*wave_vector_)  &
          ,vel_v(i,j,k,2)                                      &
          ,u0_*cos(pi*depth_u(i,j,k)/h_u(i,j))                 &
              *cos(omega_*elapsedtime_now-dist_*wave_vector_)  &
              *coriolis_t(1,1)/omega_


        enddo
       close(3)
      endif                                     !------------------>

      open(unit=3, &
      file=trim(tmpdirname)//''//dom_c//'_t_.dat',position='append')
      i=imax/2 ; j=jmax/2 ; k=kmax
      dist_=real(i+par%timax(1)-i0)*dxb ! distance from the wave generation point
        write(3,'(3(e13.6,1x))')elapsedtime_now/86400.  &
                               ,real(vel_u(i,j,k,1))   &
           ,u0_*cos(pi*depth_u(i,j,k)/h_u(i,j))                 &
               *sin(omega_*elapsedtime_now-dist_*wave_vector_)
      close(3)

!      sum1=0.
!      j=jmax/2
!      k=kmax
!      do i=3,imax ! commence 3 pour ne pas repeter mpi
!      dist_=real(i+par%timax(1)-i0)*dxb ! distance from the wave generation point
!      sum1=sum1+abs(                                           &
!           vel_u(i,j,k,2)                                      &
!          -u0_*cos(pi*depth_u(i,j,k)/h_u(i,j))                 &
!              *sin(omega_*elapsedtime_now-dist_*wave_vector_)  )
!      enddo
!#ifdef parellele
!       call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,     &
!            mpi_sum,par%comm2d ,ierr)
!       sum1=sum1glb
!#endif


      end subroutine internal_mode_wave2dv_generator
#endif
!..................................................................
