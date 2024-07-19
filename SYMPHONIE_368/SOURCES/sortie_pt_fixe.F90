      subroutine sortie_pt_fixe

! simplification par rapport a version emilie
! ecriture du point du model_1D seulement
! ATTENTION : TEMPS EN JOUR
! LE FORMAT DU TEMPS A CHANGE



      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='sortie_pt_fixe'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      iptbio=max0(min0(iptbio,imax),1) ! NE SURTOUT PAS OTER LE MIN0
      jptbio=max0(min0(jptbio,jmax),1) ! SINON DEPASSEMENT DE MEMOIRE ASSURE!
      write(6,*)iptbio,jptbio,'iptbio jptbio'

      if(flag3dbio.eq.1)then      !°°°°°°°°°°°°°°°°°°°°°°°
       open(unit=34,file=                                               &
       dirgraph(1:lname4)//'bio1d.res'                                  &
        ,position='append')

      do k=kmin_w(iptbio,jptbio),kmax
	
       write(34,100)elapsedtime_now*un/86400.                    &
        , depth_t(iptbio,jptbio,k)                                      &
        ,tem_t(iptbio,jptbio,k,1)                                       &
        ,(bio_t(iptbio,jptbio,k,vb),vb=1,vbmax),                        &
              kh_w(iptbio,jptbio,k)

      enddo
      close(34)
      else                     !°°°°°°°°°°°°°°°°°°°°°°°
       open(unit=34,file=                                               &
       dirgraph(1:lname4)//'bio3d1.res'                                 &
        ,position='append')

      iptbio=min0(52,imax) ! NE SURTOUT PAS OTER LE MIN0
      jptbio=min0(41,jmax) ! SINON DEPASSEMENT DE MEMOIRE ASSURE!

      do k=kmin_w(iptbio,jptbio),kmax

       write(34,100)elapsedtime_now*un/86400.                    &
        , depth_t(iptbio,jptbio,k)                                      &
        ,tem_t(iptbio,jptbio,k,1)                                       &
        ,(bio_t(iptbio,jptbio,k,vb),vb=1,vbmax),                        &
              kh_w(iptbio,jptbio,k)

      enddo
      close(34)
       open(unit=34,file=                                               &
       dirgraph(1:lname4)//'bio3d2.res'                                 &
        ,position='append')

      iptbio=min0(52,imax) ! NE SURTOUT PAS OTER LE MIN0
      jptbio=min0(41,jmax) ! SINON DEPASSEMENT DE MEMOIRE ASSURE!
      do k=kmin_w(iptbio,jptbio),kmax

       write(34,100)elapsedtime_now*un/86400.                    &
        , depth_t(iptbio,jptbio,k)                                      &
        ,tem_t(iptbio,jptbio,k,1)                                       &
        ,(bio_t(iptbio,jptbio,k,vb),vb=1,vbmax),                        &
              kh_w(iptbio,jptbio,k)

      enddo
      close(34)
       open(unit=34,file=                                               &
       dirgraph(1:lname4)//'bio3d3.res'                                 &
        ,position='append')

      iptbio=min0(52,imax) ! NE SURTOUT PAS OTER LE MIN0
      jptbio=min0(41,jmax) ! SINON DEPASSEMENT DE MEMOIRE ASSURE!
      do k=kmin_w(iptbio,jptbio),kmax

       write(34,100)elapsedtime_now*un/86400.                    &
        , depth_t(iptbio,jptbio,k)                                      &
        ,tem_t(iptbio,jptbio,k,1)                                       &
        ,(bio_t(iptbio,jptbio,k,vb),vb=1,vbmax),                        &
              kh_w(iptbio,jptbio,k)

      enddo
      close(34)
       open(unit=34,file=                                               &
       dirgraph(1:lname4)//'bio3d4.res'                                 &
        ,position='append')

      iptbio=min0(52,imax) ! NE SURTOUT PAS OTER LE MIN0
      jptbio=min0(41,jmax) ! SINON DEPASSEMENT DE MEMOIRE ASSURE!

      do k=kmin_w(iptbio,jptbio),kmax

       write(34,100)elapsedtime_now*un/86400.                    &
        , depth_t(iptbio,jptbio,k)                                      &
        ,tem_t(iptbio,jptbio,k,1)                                       &
        ,(bio_t(iptbio,jptbio,k,vb),vb=1,vbmax),                        &
              kh_w(iptbio,jptbio,k)

      enddo
      close(34)

      endif                    !°°°°°°°°°°°°°°°°°°°°°°°

 100  format(f7.2,1x,f11.4,1x,f13.5,1x,30(e11.4))

      end subroutine sortie_pt_fixe
