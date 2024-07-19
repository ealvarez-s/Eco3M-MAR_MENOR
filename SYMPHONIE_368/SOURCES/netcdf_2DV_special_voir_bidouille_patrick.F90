      subroutine netcdf_main(txt2_   )
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26 - last update: 22-07-18
!______________________________________________________________________
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      integer :: nf_type_,flag_=0 , loop_ &
                                  , loopmax_=1
      character*2,intent(in) :: txt2_
#ifdef synopsis
       subroutinetitle='netcdf_main'
       subroutinedescription= &
       'Driver of the subroutines writting the netcdf output files'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
! 2010.11 13-07-10  Mise en service
! 2010.20 17-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 30-04-11  ecriture grille globale
! 2010.23 14-05-11  Ne pas appeler la routine netcdf_obc dans certains cas
! 2010.25 04-03-12  ajout variables dimensions
!         17-03-12  commencer a prevoir d'archiver les lon lat en double precision
!         27-03-12  debug varid(count_netcdfcar)
!         03-04-12  Ajout subroutine netcdf_general_attributes
!                   imax jmax et cie remplac*s par ni nj etc...
!         21-06-12  Uilisation de pnetcdf au lieu de netcdf classique
! S26     04-02-13  test sur proc_null remplac* par test sur i et j sinon
!                   bug sur cas periodique
!         06-02-13  ecrire lon lat des points "f"
!         29-09-13  Pourvoir ecrire avec type short
!         01-11-13  traitement des bords du domaine
!         02-05-14  Possibilite d'ecrire des tableaux en double precision
!         26-06-14  Ne plus passer par le masquage des points OBC
!         05-07-14  Cas des procs inactifs elimines. Les espaces vacants du fichiersont
!                   netcdf sont traites par les procs actifs
!         17-07-14  possibilite de double precision pour variables TYX
!         01-10-14  passage par anyv3dint
!         07-10-14  debug suite au point precedent
!         24-10-14  debug suite au point precedent
!         11-10-14  traitement des bords Ha Hg
!         08-12-14  ajout tableaux "f" 3D
!         12-12-14  traitement des bords depth_u depth_v
!         20-12-14  traitement des bords bio
!         09-01-15  traitement des bords kmin_
!         24-01-15  texte80(11)='mask_obc_z1' pour masquer c.l. z1
!         06-02-15  Ne pas masquer Ha Hg sur les bords
!         21-03-15  alarme depassement de memoire
!         16-05-15  suite point precedent: + de message a l'ecran
!         07-11-15  variable err_h masquee obc
!         19-11-15  variable delta_h_hogcm masque obc
!         23-12-15  champs CLA masquEs
!         13-04-16  Sur remarque de Florent Lyard suppression des bords de grille
!                   masquEs pour longitude latitude des noeuds u,v,f
!         05-09-16  modif pour les cas avec plus de coeurs passifs que de coeurs actifs
!         13-01-17  AJouter les auteurs de la simu dans les attributs des fichiers netcdf
!         24-01-17  Ajouter le nom du schema de turbulence, des formules bulk
!         28-02-17  Ajout d'une URL
!         21-03-17  Message ecran d'aide au debugage
!         05-05-17  Attributs generaux supplementaires
!         22-05-17  Attributs generaux supplementaires
!         29-03-18  seul rank0 lit le fichier .....
!         28-04-18  Ajout d'attributs globaux
!         15-06-18  loopmax_ permet d'ecrire le fichier netcdf en plusieurs
!                   etapes pour que, si besoin, tous les coeurs n'ecrivent pas
!                   en meme temps dans le fichier netcdf
!         22-06-18  ajout des modules utilises par la simulation dans les attributs generaux
!         05-07-18  ajout de la valeur de status dans les messages d'erreur
!         22-07-18  mise A jour lien https dans attributs globaux
!...............................................................................
!  _________                    .__                  .__              !m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____        !
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!......................................................................

      texte80(6)='undefined'

      if(texte80(5)(1:2)=='YX') then
         if(txt2_   =='_t')texte80(6)='latitude_t longitude_t'
         if(txt2_   =='_w')texte80(6)='latitude_t longitude_t'
         if(txt2_   =='_u')texte80(6)='latitude_u longitude_u'
         if(txt2_   =='_v')texte80(6)='latitude_v longitude_v'
         if(txt2_   =='_f')texte80(6)='latitude_f longitude_f'
      endif
      if(texte80(5)(1:3)=='TYX') then
         if(txt2_   =='_t')texte80(6)='time latitude_t longitude_t'
         if(txt2_   =='_w')texte80(6)='time latitude_t longitude_t'
         if(txt2_   =='_u')texte80(6)='time latitude_u longitude_u'
         if(txt2_   =='_v')texte80(6)='time latitude_v longitude_v'
         if(txt2_   =='_f')texte80(6)='time latitude_f longitude_f'
      endif
      if(texte80(5)(1:3)=='ZYX') then
         if(txt2_   =='_t')texte80(6)='depth_t latitude_t longitude_t'
         if(txt2_   =='_w')texte80(6)='depth_w latitude_t longitude_t'
         if(txt2_   =='_u')texte80(6)='depth_u latitude_u longitude_u'
         if(txt2_   =='_v')texte80(6)='depth_v latitude_v longitude_v'
         if(txt2_   =='_f')texte80(6)='depth_f latitude_f longitude_f' !08-12-14
      endif
      if(texte80(5)(1:4)=='TZYX') then
      if(txt2_   =='_t')texte80(6)='time depth_t latitude_t longitude_t'
      if(txt2_   =='_w')texte80(6)='time depth_w latitude_t longitude_t'
      if(txt2_   =='_u')texte80(6)='time depth_u latitude_u longitude_u'
      if(txt2_   =='_v')texte80(6)='time depth_v latitude_v longitude_v'
      endif

      ksecu=0
      if(texte80(7)=='real')   then ; ksecu=1 ; nf_type_=nf_real   ; endif
      if(texte80(7)=='double') then ; ksecu=1 ; nf_type_=nf_double ; endif
      if(texte80(7)=='integer')then ; ksecu=1 ; nf_type_=nf_int    ; endif
      if(texte80(7)=='short')  then ; ksecu=1 ; nf_type_=nf_short  ; endif !27-09-13
      if(ksecu==0) then
       write(6,'(a,a)')'Pb avec variable netcdf: ',trim(texte80(1))
       write(6,'(a,a)')'de type ',trim(texte80(7))
       stop ' stop netcdf_defvar type inconnu'
      endif

      count_netcdfvar=count_netcdfvar+1       ! Compteur de variables
      if(count_netcdfvar>dim_varid) then !>>>>>!19-04-15
       write(6,*)'count_netcdfvar,dim_varid',count_netcdfvar,dim_varid !16-05-15
       Stop 'Err count_netcdfvar>dim_varid (defined in notebook_graph)'
      endif

      if(loop_netcdf==0)call netcdf_defvar(txt2_,nf_type_)

      if(loop_netcdf==1) then !111111111111111>

        do loop_=1,loopmax_  !PPPPPPPPPPPP>

            if(mod(par%rank+loop_,loopmax_)==0) then
             flag_write=1
            else
             flag_write=0
            endif

        call netcdf_write(txt2_,nf_type_   &
                         ,3                & ! imax             & ! bidouille patrick
                         ,jmax             &
                         ,par%coords(1)    & ! par%timax(1)     & ! bidouille patrick
                         ,par%tjmax(1)     &
                         ,0)

! Ensuite bouchage (avec valeur de masquage) des domaines sans proc !05-07-14
! Attention, manifestement avec pnetcdf, il n'est pas possible que
! seuls certains proc entrent dans la routine (car les autres se
! feraient indefiniment attendre par le processus de rdv mpi).
! La strategie est donc que tous les procs ecrivent dans les trous
! meme si cela doit impliquer une redondance d'ecriture...
         if(mpi_hole_plugging/='none') then ! hp hp hp hp >
          i1=-9999
          flag_=0
          open(unit=3,file=mpi_hole_plugging)
           read(3,*,end=42)

   87      read(3,*,end=42)k1 
           read(3,*,end=42)i3,i4 ! coordonnees zone bouchage
           read(3,*,end=42)j3,j4 ! coordonnees zone bouchage
           read(3,'(a)',end=42)texte30

           if(par%rank==k1) then
             flag_=1 ; i1=i3 ; i2=i4 ; j1=j3 ; j2=j4 
           endif

           if(texte30(1:7)=='PERFORM') then !pmx> !05-09-16
             if(flag_==0) then
              i1=i3 ; i2=i4 ; j1=j3 ; j2=j4 
             endif
! Avec cet algo tous les proc non designes pour boucher boucheront quand meme
! le dernier trou... Si i1=-9999 c'est qu'il n'y a aucun proc a boucher
             if(i1/=-9999)call netcdf_write(txt2_,nf_type_,i2-i1,j2-j1,i1,j1,1)
             read(3,'(a)',end=42)texte30
             if(texte30(1:7)/=' ......') then !>>>
               write(6,'(a)')trim(texte30)
               stop 'Erreur sur texte30'
             endif                            !>>>
           endif                            !pmx>

           goto 87

 42       close(3)
         endif                              ! hp hp hp hp >

        enddo         !PPPPPPPPPPPP>
      endif                   !111111111111111>

      texte80(3)='none'
      texte80(4)='none'

      end subroutine netcdf_main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine netcdf_defvar(txt2_,nf_type_)
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      integer nf_type_
      character*2,intent(in) :: txt2_
      character*12 dynamic_range_txt_
      integer(kind=MPI_OFFSET_KIND) :: len_
      real :: scalarval(1)
      double precision :: scalarvalr8(1)
      integer :: scalarvali4(1)
#ifdef synopsis
       subroutinetitle='netcdf_defvar'
       subroutinedescription='Creates the netcdf file header'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      cgridshift(1)=0. ; cgridshift(2)=0. ; cgridshift(3)=0.
      if(txt2_   =='_u')cgridshift(1)= 0.5
      if(txt2_   =='_f')cgridshift(1)= 0.5
      if(txt2_   =='_v')cgridshift(2)= 0.5
      if(txt2_   =='_f')cgridshift(2)= 0.5
      if(txt2_   =='_w')cgridshift(3)=-0.5

      ksecu=0
      if(texte80(1)=='time')ksecu=1
      if(texte80(1)=='cumulativetime')ksecu=1

      if(texte80(6)(1:9)=='undefined') then      !>>>>>>>>
       if(texte80(5)=='X') then
            ksecu=1 ; k0=1
            if(txt2_   =='_t')vardim(1)=i_t_dim
            if(txt2_   =='_u')vardim(1)=i_u_dim
            if(txt2_   =='_v')vardim(1)=i_v_dim
            if(txt2_   =='_w')vardim(1)=i_w_dim
            if(txt2_   =='_f')vardim(1)=i_f_dim
       endif
       if(texte80(5)=='Y') then
            ksecu=1 ; k0=1
            if(txt2_   =='_t')vardim(1)=j_t_dim
            if(txt2_   =='_u')vardim(1)=j_u_dim
            if(txt2_   =='_v')vardim(1)=j_v_dim
            if(txt2_   =='_w')vardim(1)=j_w_dim
            if(txt2_   =='_f')vardim(1)=j_f_dim
       endif
       if(texte80(5)=='Z') then
            ksecu=1 ; k0=1
            if(txt2_   =='_t')vardim(1)=k_t_dim
            if(txt2_   =='_u')vardim(1)=k_u_dim
            if(txt2_   =='_v')vardim(1)=k_v_dim
            if(txt2_   =='_w')vardim(1)=k_w_dim
            if(txt2_   =='_f')vardim(1)=k_f_dim
       endif
      endif                                      !>>>>>>>>

      if(texte80(6)(1:22)=='latitude_t longitude_t') then      !>>>>>>>>
       ksecu=1
       k0=2
       vardim(1)=i_t_dim
       vardim(2)=j_t_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:22)=='latitude_u longitude_u') then      !>>>>>>>>
       ksecu=1
       k0=2
       vardim(1)=i_u_dim
       vardim(2)=j_u_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:22)=='latitude_v longitude_v') then      !>>>>>>>>
       ksecu=1
       k0=2
       vardim(1)=i_v_dim
       vardim(2)=j_v_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:22)=='latitude_f longitude_f') then      !>>>>>>>>
       ksecu=1
       k0=2
       vardim(1)=i_f_dim
       vardim(2)=j_f_dim
      endif                                             !>>>>>>>>

      if(texte80(6)(1:27)=='time latitude_t longitude_t') then !>>>>>>>>
       ksecu=1
       k0=3
       vardim(1)=i_t_dim
       vardim(2)=j_t_dim
       vardim(3)=time_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:27)=='time latitude_u longitude_u') then !>>>>>>>>
       ksecu=1
       k0=3
       vardim(1)=i_u_dim
       vardim(2)=j_u_dim
       vardim(3)=time_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:27)=='time latitude_v longitude_v') then !>>>>>>>>
       ksecu=1
       k0=3
       vardim(1)=i_v_dim
       vardim(2)=j_v_dim
       vardim(3)=time_dim
      endif                                             !>>>>>>>>
      if(texte80(6)(1:27)=='time latitude_f longitude_f') then !>>>>>>>>
       ksecu=1
       k0=3
       vardim(1)=i_f_dim
       vardim(2)=j_f_dim
       vardim(3)=time_dim
      endif                                             !>>>>>>>>

      if(texte80(6)(1:30)=='depth_t latitude_t longitude_t')then !>>>>
       ksecu=1
       k0=3 ; vardim(1)=i_t_dim ; vardim(2)=j_t_dim ; vardim(3)=k_t_dim
      endif
      if(texte80(6)(1:30)=='depth_u latitude_u longitude_u')then !>>>>
       ksecu=1
       k0=3 ; vardim(1)=i_u_dim ; vardim(2)=j_u_dim ; vardim(3)=k_u_dim
      endif
      if(texte80(6)(1:30)=='depth_v latitude_v longitude_v')then !>>>>
       ksecu=1
       k0=3 ; vardim(1)=i_v_dim ; vardim(2)=j_v_dim ; vardim(3)=k_v_dim
      endif
      if(texte80(6)(1:30)=='depth_w latitude_t longitude_t')then !>>>>
       ksecu=1
       k0=3 ; vardim(1)=i_w_dim ; vardim(2)=j_w_dim ; vardim(3)=k_w_dim
      endif
      if(texte80(6)(1:30)=='depth_f latitude_f longitude_f')then !>>>>
       ksecu=1
       k0=3 ; vardim(1)=i_f_dim ; vardim(2)=j_f_dim ; vardim(3)=k_w_dim
      endif


      if(texte80(6)(1:35)=='time depth_t latitude_t longitude_t')then !>>>>
       ksecu=1
       k0=4 ; vardim(1)=i_t_dim ; vardim(2)=j_t_dim                ! 4D ; dim1 ; dim2
       vardim(3)=k_t_dim        ; vardim(4)=time_dim                 ! dim3 ; dim4
      endif
      if(texte80(6)(1:35)=='time depth_u latitude_u longitude_u')then !>>>>
       ksecu=1
       k0=4 ; vardim(1)=i_u_dim ; vardim(2)=j_u_dim                ! 4D ; dim1 ; dim2
       vardim(3)=k_u_dim        ; vardim(4)=time_dim                 ! dim3 ; dim4
      endif
      if(texte80(6)(1:35)=='time depth_v latitude_v longitude_v')then !>>>>
       ksecu=1
       k0=4 ; vardim(1)=i_v_dim ; vardim(2)=j_v_dim                ! 4D ; dim1 ; dim2
       vardim(3)=k_v_dim        ; vardim(4)=time_dim                 ! dim3 ; dim4
      endif
      if(texte80(6)(1:35)=='time depth_w latitude_t longitude_t')then !>>>>
       ksecu=1
       k0=4 ; vardim(1)=i_w_dim ; vardim(2)=j_w_dim                ! 4D ; dim1 ; dim2
       vardim(3)=k_w_dim        ; vardim(4)=time_dim                 ! dim3 ; dim4
      endif

      if(ksecu==0)then
       write(6,*)'Numero de variable ',count_netcdfvar
       write(6,*)trim(texte80(1))
       write(6,*)trim(texte80(6))
       stop 'Echec 1 dans netcdf_defvar'
      endif

      if(texte80(3)(1:4)=='none')texte80(3)=texte80(1)
      if(texte80(4)(1:4)=='none')texte80(4)=texte80(1)



      if(    texte80(1)=='time'                &
         .or.texte80(1)=='cumulativetime') then  !***************>

      if(texte80(7)=='double')                          &
      status=nfmpi_def_var(ncid,texte80(1),nf_double       &
                            ,k0,vardim,varid(count_netcdfvar))
      if(status/=0)call netcdf_error(1,1)

      len_= len_trim(texte80(2))
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'units',len_,texte80(2))
      if(status/=0)call netcdf_error(1,2)

      len_=len_trim(texte80(3))
      if(texte80(3)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'long_name',len_,texte80(3))
      if(status/=0)call netcdf_error(1,3)

      len_=len_trim(texte80(4))
      if(texte80(4)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'standard_name',len_,texte80(4))
      if(status/=0)call netcdf_error(1,4)

      if(texte80(1)=='time') then !tttttttttttttt>

        len_=len_trim(texte80(8))
        status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
              ,'time_origin',len_,texte80(8))
        if(status/=0)call netcdf_error(1,8)

        len_=len_trim(texte80(9))
        status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
              ,'calendar',len_,texte80(9))
        if(status/=0)call netcdf_error(1,9)

      endif                       !tttttttttttttt>

      len_=len_trim(texte80(5))
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'content',len_,texte80(5))
      if(status/=0)call netcdf_error(1,5)

      len_=len_trim(texte80(5))
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'axis',len_,texte80(5))                !05-03-10
      if(status/=0)call netcdf_error(1,5)

      len_=len_trim(texte80(6))
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'associate',len_,texte80(6))
      if(status/=0)call netcdf_error(1,6)


      else                              !***************>

      status=nfmpi_def_var(ncid,texte80(1),nf_type_              &
                            ,k0,vardim,varid(count_netcdfvar))
      if(status/=0)call netcdf_error(1,1)

      len_= len_trim(texte80(2))
      if(texte80(2)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'units',len_,texte80(2))
      if(status/=0)call netcdf_error(1,2)

      len_=len_trim(texte80(3))
      if(texte80(3)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'long_name',len_,texte80(3))
      if(status/=0)call netcdf_error(1,3)

      len_=len_trim(texte80(4))
      if(texte80(4)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)           &
            ,'standard_name',len_,texte80(4))
      if(status/=0)call netcdf_error(1,4)

      len_=len_trim(texte80(5))
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'content',len_,texte80(5))                !05-03-10
      if(status/=0)call netcdf_error(1,5)

      len_=len_trim(texte80(5))
      if(texte80(5)=='X'.or.texte80(5)=='Y'.or.texte80(5)=='Z')   &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'axis',len_,texte80(5))                !05-03-10
      if(status/=0)call netcdf_error(1,5)

      len_=len_trim(texte80(6))
      if(texte80(6)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'associate',len_,texte80(6))
      if(status/=0)call netcdf_error(1,6)

      len_=len_trim(texte80(6))
      if(texte80(6)/='undefined') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'coordinates',len_,texte80(6))
      if(status/=0)call netcdf_error(1,12)

      len_=1
      if(texte80(6)/='undefined') then !------------------->

      scalarval(1)=filval
      scalarvalr8(1)=filval

      if(nf_type_==nf_double)   &
       status=nfmpi_put_att_double(ncid,varid(count_netcdfvar) &
              ,'_FillValue',nf_type_,len_,scalarvalr8)
      if(nf_type_==nf_real)   &
       status=nfmpi_put_att_real(ncid,varid(count_netcdfvar) &
              ,'_FillValue',nf_type_,len_,scalarval)
      scalarvali4=nint(scalarval)
      if(nf_type_==nf_short) &
       status=nfmpi_put_att_int(ncid,varid(count_netcdfvar)  &    !27-09-13
              ,'_FillValue',nf_type_,len_,scalarvali4)

      if(status/=0)call netcdf_error(1,11)


      if(nf_type_==nf_short) then !.......>

       scalarval(1)=var_scalefactor
       status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)  &    !27-09-13
              ,'scale_factor',nf_real,len_,scalarval(1))
       if(status/=0)call netcdf_error(1,17)

       scalarval(1)=var_addoffset
       status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)  &    !27-09-13
              ,'add_offset',nf_real,len_,scalarval(1))
       if(status/=0)call netcdf_error(1,18)

!      scalarval(1)=var_validmin
!      status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)  &    !27-09-13
!             ,'valid_min',nf_real,len_,scalarval(1))
!      if(status/=0)call netcdf_error(1,19)

!      scalarval(1)=var_validmax
!      status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)  &    !27-09-13
!             ,'valid_max',nf_real,len_,scalarval(1))
!      if(status/=0)call netcdf_error(1,20)

      len_=len_trim(texte80(10))
      if(texte80(10)/='none') &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'comment',len_,texte80(10))
      if(status/=0)call netcdf_error(1,21)

      endif                       !.......>

      endif                            !------------------->

#ifdef comodo
      if(    texte80(1)(1:3)=='ni_'             &
         .or.texte80(1)(1:3)=='nj_'             &
         .or.texte80(1)(1:3)=='nk_') then       !NNNNNNNNNNNNNNNN>

      if(cgridshift(1)/=0.) then !ccccccccccccccccc>
      status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)         &
            ,'c_grid_axis_shift',nf_real,1,cgridshift(1))
      if(status/=0)call netcdf_error(1,13)
      endif                      !ccccccccccccccccc>

      if(cgridshift(2)/=0.) then !ccccccccccccccccc>
      status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)         &
            ,'c_grid_axis_shift',nf_real,1,cgridshift(2))
      if(status/=0)call netcdf_error(1,14)
      endif                      !ccccccccccccccccc>

      if(cgridshift(3)/=0.) then !ccccccccccccccccc>
      if(index(texte80(5),"Z")/=0) then !->->->
       status=nfmpi_put_att_real(ncid,varid(count_netcdfvar)       &
             ,'c_grid_axis_shift',nf_real,1,cgridshift(3))
       if(status/=0)call netcdf_error(1,15)
      endif                                 !->->->
      endif                      !ccccccccccccccccc>

      if(texte80(1)(1:4)=='ni_t')write(dynamic_range_txt_,'(i0":"i0)')2,iglb-1
      if(texte80(1)(1:4)=='ni_w')write(dynamic_range_txt_,'(i0":"i0)')2,iglb-1
      if(texte80(1)(1:4)=='ni_u')write(dynamic_range_txt_,'(i0"."i0":"i0"."i0)')2,5,iglb-2,5
      if(texte80(1)(1:4)=='ni_v')write(dynamic_range_txt_,'(i0":"i0)')2,iglb-1

      if(texte80(1)(1:4)=='nj_t')write(dynamic_range_txt_,'(i0":"i0)')2,jglb-1
      if(texte80(1)(1:4)=='nj_w')write(dynamic_range_txt_,'(i0":"i0)')2,jglb-1
      if(texte80(1)(1:4)=='nj_u')write(dynamic_range_txt_,'(i0":"i0)')2,jglb-1
      if(texte80(1)(1:4)=='nj_v')write(dynamic_range_txt_,'(i0"."i0":"i0"."i0)')2,5,jglb-2,5

      if(texte80(1)(1:4)=='nk_t')write(dynamic_range_txt_,'(i0":"i0)')1,kmax
      if(te
      if(texte80(1)(1:4)=='nk_u')write(dynamic_range_txt_,'(i0":"i0)')1,kmax
      if(texte80(1)(1:4)=='nk_v')write(dynamic_range_txt_,'(i0":"i0)')1,kmax

      len_=len_trim(dynamic_range_txt_)
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)         &
            ,'c_grid_dynamic_range',len_,dynamic_range_txt_)

      endif                                     !NNNNNNNNNNNNNNNN>
#endif

      len_=2
      if(index(texte80(5),'Z')/=0) &
      status=nfmpi_put_att_text(ncid,varid(count_netcdfvar)          &
            ,'positive',len_,'up')
      if(status/=0)call netcdf_error(1,16)

      endif                             !***************>


      end subroutine netcdf_defvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine netcdf_write(txt2_,nf_type_,imax_,jmax_,partimax1_,partjmax1_,flag_dom_)
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      integer,intent(in) :: imax_,jmax_,flag_dom_,partimax1_,partjmax1_
      character*2,intent(in) :: txt2_
       integer(kind=MPI_OFFSET_KIND) start(4)
       integer(kind=MPI_OFFSET_KIND) edge(4)
       integer(kind=MPI_OFFSET_KIND) mid(4)
       integer nf_type_,kmax_                           !27-09-13
       double precision :: val(1)
#ifdef synopsis
       subroutinetitle='netcdf_write'
       subroutinedescription= &
       'Writes the variable field in the netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Si seul le rank 0 cree varid(count_netcdfvar) en phase de definition
! du fichier netcdf, lors de l'ecriture les autres proc doivent le definir
! a prtir du nom de la variable et de la fonction nf_inq_varid:
!      status=nf_inq_varid(ncid,texte80(1),varid(count_netcdfvar))      !27-03-12
!      if(status/=0)stop 'error netcdf_write get varid(count_netcdfvar)'
!! Plus necessaire ne pnetcdf, tous les proc definissent les varibles...
      status=0

      start(1)=1 ; start(2)=1 ; start(3)=1 ; start(4)=1

      kmax_   =kmax

      if(txt2_   =='_w')kmax_   =kmax+1

! t or w nodes
      if(txt2_   =='_t'.or.txt2_   =='_w') then !ttttttttttttt>
       istr=1                      ; jstr=1
       iend=imax_-2                ; jend=jmax_-2
       if(partimax1_==0)istr=0               !04-02-13
       if(imax_+partimax1_==iglb)iend=imax_+1
       if(partjmax1_==0)jstr=0
       if(jmax_+partjmax1_==jglb)jend=jmax_+1
       start(1)=1+partimax1_+istr ; start(2)=1+partjmax1_+jstr
      endif                                     !ttttttttttttt>
! u nodes
      if(txt2_   =='_u') then !uuuuuuuuuuuuu>
       istr=1                      ; jstr=1
       iend=imax_-2                ; jend=jmax_-2
       if(partimax1_==0)istr=1
       if(imax_+partimax1_==iglb)iend=imax_+1
       if(partjmax1_==0)jstr=0
       if(jmax_+partjmax1_==jglb)jend=jmax_+1
       start(1)=1+partimax1_      ; start(2)=1+partjmax1_+jstr
      endif                   !uuuuuuuuuuuuu>
! v nodes
      if(txt2_   =='_v') then !vvvvvvvvvvvvv>
       istr=1                      ; jstr=1
       iend=imax_-2                ; jend=jmax_-2
       if(partimax1_==0)istr=0
       if(imax_+partimax1_==iglb)iend=imax_+1
       if(partjmax1_==0)jstr=1
       if(jmax_+partjmax1_==jglb)jend=jmax_+1
       start(1)=1+partimax1_+istr ; start(2)=1+partjmax1_
      endif                   !vvvvvvvvvvvvv>
! f nodes
      if(txt2_   =='_f') then !vvvvvvvvvvvvv>
       istr=1                      ; jstr=1
       iend=imax_-2                ; jend=jmax_-2
       if(partimax1_==0)istr=1
       if(imax_+partimax1_==iglb)iend=imax_+1
       if(partjmax1_==0)jstr=1
       if(jmax_+partjmax1_==jglb)jend=jmax_+1
       start(1)=1+partimax1_      ; start(2)=1+partjmax1_
      endif                   !vvvvvvvvvvvvv>


      edge(1)=iend-istr+1         ; edge(2)=jend-jstr+1

      if(flag_write==0) then !oooo>
!       iend=istr ; jend=jstr ; edge(1)=1 ; edge(2)=1 
        iend=istr ; jend=jstr ; edge(1)=0 ; edge(2)=0 
      endif                  !oooo>


      if(texte80(1)=='time') then !tttttttt>
       x1=elapsedtime_now
       start(1)=1 ; edge(1)=1 ; val(1)=x1 ; status=0
       if(par%rank/=0)edge(1)=0
       status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar),start,edge,val)           !16-11-09
       if(status/=0)stop ' stop erreur netcdf_write_par time'
       return
      endif                       !tttttttt>

      if(texte80(1)=='cumulativetime') then !ccccccc>
       x1=ofl_period_now/3600.
       start(1)=1 ; edge(1)=1 ; val(1)=x1 ; status=0
       if(par%rank/=0)edge(1)=0
       status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar),start,edge,val)           !16-11-09
       if(status/=0)stop ' stop erreur netcdf_write_par time'
       return
      endif                                 !ccccccc>

      if(texte80(5)=='X') then !aaaaaaaaaaaaaaaaaaaaaa>
       dimend=1 ; ksecu=1
       if(texte80(7)=='integer') then !---->
        anyv3dint(istr:iend,1,1)=nint(anyv3d(istr:iend,1,1,1))    !01-10-14
        status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar) &
                    ,start(dimend:dimend),edge(dimend:dimend)     &
                                  ,anyv3dint(istr:iend,1,1))      !01-10-14
!                               ,nint(anyv3d(istr:iend,1,1,1)))
       endif                          !---->

       if(texte80(7)=='real') then   !--->
        anyvar3d(istr:iend,1,1)=real(anyv3d(istr:iend,1,1,1))     !01-10-14
        status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar) &
                     ,start(dimend:dimend),edge(dimend:dimend)     &
                               , anyvar3d(istr:iend,1,1) )        !01-10-14
!                              , real( anyv3d(istr:iend,1,1,1) ))
       endif                         !--->
      endif                    !aaaaaaaaaaaaaaaaaaaaaa>
      if(texte80(5)=='Y') then !aaaaaaaaaaaaaaaaaaaaaa>
       dimend=2 ; ksecu=1
       if(texte80(7)=='integer') then !--->
        anyv3dint(1,jstr:jend,1)=nint(anyv3d(1,jstr:jend,1,2))
        status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar) &
                     ,start(dimend:dimend),edge(dimend:dimend)    &
                                , anyv3dint(1,jstr:jend,1) )        !01-10-14
!                               ,nint(anyv3d(1,jstr:jend,1,2)))
       endif                          !--->
       if(texte80(7)=='real') then !--->
        anyvar3d(1,jstr:jend,1)=real(anyv3d(1,jstr:jend,1,2)) !01-10-14
        status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar) &
                     ,start(dimend:dimend),edge(dimend:dimend)     &
                                ,anyvar3d(1,jstr:jend,1) )    !01-10-14
!                               ,real( anyv3d(1,jstr:jend,1,2) ))
       endif                       !--->
      endif                    !aaaaaaaaaaaaaaaaaaaaaa>
      if(texte80(5)=='Z') then !aaaaaaaaaaaaaaaaaaaaaa>
       kend=kmax_    ; dimend=3 ; edge(3)=kend ; ksecu=1

!      if(flag_write==0) then ; kend=1 ; edge(3)=1 ; endif
       if(flag_write==0) then ; kend=1 ; edge(3)=0 ; endif

       if(texte80(7)=='integer') then !--->
        anyv3dint(1,1,1:kend)=nint(anyv3d(1,1,1:kend,3)) !01-10-14
        status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar) &
                     ,start(dimend:dimend),edge(dimend:dimend)    &
                                 ,anyv3dint(1,1,1:kend)) !01-10-14
!                                ,nint(anyv3d(1,1,1:kend,3)))
       endif                          !--->
       if(texte80(7)=='real') then !--->
        anyvar3d(1,1,1:kend)=real(anyv3d(1,1,1:kend,3)) !01-10-14
        status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar) &
                     ,start(dimend:dimend),edge(dimend:dimend)     &
                                , anyvar3d(1,1,1:kend)) !01-10-14
!                               , real( anyv3d(1,1,1:kend,3) ))
       endif                       !--->
      endif                    !aaaaaaaaaaaaaaaaaaaaaa>

      if(texte80(5)(1:2)=='YX') then !aaaaaaaaaaaaaaaaaaaaaa>
      kend=1                                                        !30-04-11
      call netcdf_obc(txt2_) ! kend doit etre connu pour appeler netcdf_obc
      dimend=2
      ksecu=0

      if(nf_type_==nf_real) then
       if(flag_dom_==1)anyvar2d(istr:iend,jstr:jend)=filval
       ksecu=1
       status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar) &
                               ,start(1:dimend),edge(1:dimend)    &
                               ,anyvar2d(istr:iend,jstr:jend))
      endif

      if(nf_type_==nf_double) then
       if(flag_dom_==1)anyv3d(istr:iend,jstr:jend,1,1)=filval
       ksecu=1
       status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar) &
                               ,start(1:dimend),edge(1:dimend)      &
                               ,anyv3d(istr:iend,jstr:jend,1,1))
      endif

      if(nf_type_==nf_short) then
       if(flag_dom_==1) then !---->
        anyv3dint(istr:iend,jstr:jend,1)=nint(filval) !24-10-14
       else                  !---->                   !24-10-14
        anyv3dint(istr:iend,jstr:jend,1)=nint(anyvar2d(istr:iend,jstr:jend)) !01-10-14
       endif                 !---->
       ksecu=1
       status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar)    &
                               ,start(1:dimend),edge(1:dimend)      &
                               ,anyv3dint(istr:iend,jstr:jend,1)) !01-10-14
!                              ,nint(anyvar2d(istr:iend,jstr:jend)) )
       if(status/=0)stop 'Erreur 615 nfmpi_put_vara_int_all'
      endif

      if(ksecu==0)then
      write(*,'(a,a)')'Pb avec variable netcdf ',trim(texte80(1))
      stop ' ERREUR fonction nfmpi_put_var dans netcdf.F90'
      endif

      if(status/=0)write(*,*)'variable numero ',count_netcdfvar
      if(status/=0)write(*,*)trim(texte80(1))

      endif                           !aaaaaaaaaaaaaaaaaaaaaa>

      if(texte80(5)(1:3)=='ZYX') then !aaaaaaaaaaaaaaaaaaaaaa>
      kend=kmax                                                !30-04-11
      if(txt2_   =='_w')kend=kmax+1                            !30-04-11
      call netcdf_obc(txt2_) ! kend doit etre connu pour appeler netcdf_obc
      dimend=3 ; edge(3)=kend
      edge(3) = kend
      ksecu=0

!      if(flag_write==0) then ; kend=1 ; edge(3)=1 ; endif
       if(flag_write==0) then ; kend=1 ; edge(3)=0 ; endif

      if(nf_type_==nf_real ) then
      if(flag_dom_==1)anyvar3d(istr:iend,jstr:jend,1:kend)=filval
       ksecu=1
       status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar)     &
             ,start(1:dimend),edge(1:dimend)                          &
             ,anyvar3d(istr:iend,jstr:jend,1:kend))

      endif

      if(nf_type_==nf_double) then
      if(flag_dom_==1)anyv3d(istr:iend,jstr:jend,1:kend,1)=filval
       ksecu=1
       status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar)   &
             ,start(1:dimend),edge(1:dimend)                          &
             ,anyv3d(istr:iend,jstr:jend,1:kend,1))
      endif

      if(nf_type_==nf_short) then
       if(flag_dom_==1) then !---->
        anyv3dint(istr:iend,jstr:jend,1:kend)=nint(filval)!24-10-14
       else                  !---->
        anyv3dint(istr:iend,jstr:jend,1:kend)=nint(anyvar3d(istr:iend,jstr:jend,1:kend)) !01-10-14
       endif                 !---->
       ksecu=1
       status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar)     &
             ,start(1:dimend),edge(1:dimend)                         &
             ,anyv3dint(istr:iend,jstr:jend,1:kend))     !01-10-14
!            ,nint(anyvar3d(istr:iend,jstr:jend,1:kend)))
      endif

      if(status/=0)write(*,*)'variable numero ',count_netcdfvar
      if(status/=0)write(*,*)trim(texte80(1))

      endif                           !aaaaaaaaaaaaaaaaaaaaaa>

      if(texte80(5)(1:3)=='TYX') then !aaaaaaaaaaaaaaaaaaaaaa>
      kend=1                                                        !30-04-11
      call netcdf_obc(txt2_) ! kend doit etre connu pour appeler netcdf_obc
      dimend=3 ; edge(3)=1
      edge(3) = 1

      if(nf_type_==nf_double) then !17-07-14
      if(flag_dom_==1)anyv3d(istr:iend,jstr:jend,1,1)=filval
      status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar)   &
                                ,start(1:dimend),edge(1:dimend)      &
                                ,anyv3d(istr:iend,jstr:jend,1,1))
      endif

      if(nf_type_==nf_real) then
      if(flag_dom_==1)anyvar2d(istr:iend,jstr:jend)=filval
      status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar)   &
                              ,start(1:dimend),edge(1:dimend)      &
                              ,anyvar2d(istr:iend,jstr:jend))
      endif

      if(nf_type_==nf_short) then
       if(flag_dom_==1) then !---->
        anyv3dint(istr:iend,jstr:jend,1)=nint(filval) !24-10-14
       else                  !---->
        anyv3dint(istr:iend,jstr:jend,1)=nint(anyvar2d(istr:iend,jstr:jend))!01-10-14
       endif                 !---->
      status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar)    &
                              ,start(1:dimend),edge(1:dimend)      &
                              ,anyv3dint(istr:iend,jstr:jend,1)) !01-10-14
!                             ,nint(anyvar2d(istr:iend,jstr:jend)))
      endif

      if(status/=0)write(*,*)'variable numero ',count_netcdfvar
      if(status/=0)write(*,*)trim(texte80(1))

      endif                           !aaaaaaaaaaaaaaaaaaaaaa>

      if(texte80(5)(1:4)=='TZYX') then !bbbbbbbbbbbbbbbbbbbbb>
      kend=kmax                                                !30-04-11
      if(txt2_   =='_w')kend=kmax+1                            !30-04-11
      call netcdf_obc(txt2_) ! kend doit etre connu pour appeler netcdf_obc
      dimend=4 ; edge(3)=kend ; edge(4)=1

!      if(flag_write==0) then ; kend=1 ; edge(3)=1 ; endif
       if(flag_write==0) then ; kend=1 ; edge(3)=0 ; endif

      status=nfmpi_inq_vartype(ncid,varid(count_netcdfvar),nf_type_)
      if(status/=0)stop 'Erreur nfmpi_inq_vartype_all in netcdf.F90'

      if(nf_type_==nf_double) then !17-07-14
      if(flag_dom_==1)anyv3d(istr:iend,jstr:jend,1:kend,1)=filval
      status=nfmpi_put_vara_double_all(ncid,varid(count_netcdfvar)   &
                                ,start(1:dimend),edge(1:dimend)      &
                                ,anyv3d(istr:iend,jstr:jend,1:kend,1))
      endif

      if(nf_type_==nf_real) then
      if(flag_dom_==1)anyvar3d(istr:iend,jstr:jend,1:kend)=filval
      status=nfmpi_put_vara_real_all(ncid,varid(count_netcdfvar)     &
            ,start(1:dimend),edge(1:dimend)                          &
            ,anyvar3d(istr:iend,jstr:jend,1:kend))
      endif

      if(nf_type_==nf_short) then
      if(flag_dom_==1) then !---->
       anyv3dint(istr:iend,jstr:jend,1:kend)=nint(filval) !24-10-14
      else                  !---->                        !24-10-14
       anyv3dint(istr:iend,jstr:jend,1:kend)=nint(anyvar3d(istr:iend,jstr:jend,1:kend))!07-10-14
      endif                 !---->
      status=nfmpi_put_vara_int_all(ncid,varid(count_netcdfvar)      &
            ,start(1:dimend),edge(1:dimend)                          &
           ,anyv3dint(istr:iend,jstr:jend,1:kend))       !01-10-14
!          ,nint(anyvar3d(istr:iend,jstr:jend,1:kend)) )
      endif

      if(status/=0)write(6,*)'variable numero ',count_netcdfvar
      if(status/=0)write(6,*)trim(texte80(1))


      endif                            !bbbbbbbbbbbbbbbbbbbbb>
!!$

      if(status/=0) then
        write(texte30,'(a,i0)')trim(tmpdirname)//'erreur',par%rank
        open(unit=3,file=texte30)
         write(3,'(a,a)')'Erreur status/=0 pour la variable: '  &
                        ,trim(texte80(1))
         write(3,'(a,i0)')'status=',status !05-07-18
         write(3,'(a,a)')'texte80(5)=',trim(texte80(5))
         write(3,'(a,a)')'texte80(1)=',trim(texte80(1))
         write(3,'(a,i0)')      'var num  = ',count_netcdfvar
         write(3,'(a,i0)')      'flag_dom_= ',flag_dom_
         write(3,'(a,i0,1x,i0)')'istr:iend= ',istr,iend
         write(3,'(a,i0,1x,i0)')'jstr:jend= ',jstr,jend
         write(3,'(a,i0)')      'kend     = ',kend
         write(3,'(a,i0)')      'imax_     = ',imax_
         write(3,'(a,i0)')      'jmax_     = ',jmax_
         write(3,'(a,i0)')      'kmax     = ',kmax
         write(3,'(a,i0)')      'dimend   = ',dimend
         write(3,'(a,4(1x,i0))')'start    = ',start(1:dimend)
         write(3,'(a,4(1x,i0))')'edge     = ',edge(1:dimend)
         write(3,'(a,i0)')'start(1)+edge(1)-1=',start(1)+edge(1)-1
         write(3,'(a,i0)')'iglb+1            =',iglb+1
         write(3,'(a,i0)')'start(2)+edge(2)-1=',start(2)+edge(2)-1
         write(3,'(a,i0)')'jglb+1            =',jglb+1
         if(dimend>2) then
         write(3,'(a,i0)')'start(3)+edge(3)-1=',start(3)+edge(3)-1
         write(3,'(a,i0)')'kend              =',kend
         endif
         if(dimend>3) then
         write(3,'(a,i0)')'start(4)+edge(4)-1=',start(4)+edge(4)-1
         endif
        close(3)
        stop ' Stop in netcdf_write_par. Voir fichiers tmp/erreurXXX'
      endif


      end subroutine netcdf_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine netcdf_dim
      use module_parallele
      use module_principal
      use pnetcdf
      implicit none
!      include 'netcdf.inc'
      integer(kind=MPI_OFFSET_KIND) :: iglbp1,iglbp2
      integer(kind=MPI_OFFSET_KIND) :: jglbp1,jglbp2
      integer(kind=MPI_OFFSET_KIND) :: kmaxp0, kmaxpp1
#ifdef synopsis
       subroutinetitle='netcdf_dim'
       subroutinedescription='Dimensions in the netcdf file header'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! bidouille patrick
!     iglbp1=iglb+1
!     iglbp2=iglb+2
      iglbp1=nbdom_imax+1
      iglbp2=nbdom_imax+2

      jglbp1=jglb+1
      jglbp2=jglb+2
      kmaxp0=kmax
      kmaxpp1=kmax+1

! dimensions grille z pour T et S
!      status=nfmpi _def_dim(ncid,'imax_t',iglbp2,i_t_dim)
      status=0
      status=nfmpi_def_dim(ncid,'ni_t',iglbp2,i_t_dim)
      if(status/=0) then
          write(6,*)'status=',status
          stop ' stop erreur netcdf_dim 1'
      endif

!      status=nfmpi_def_dim(ncid,'jmax_t',jglbp2,j_t_dim)
      status=nfmpi_def_dim(ncid,'nj_t',jglbp2,j_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 2'

!      status=nfmpi_def_dim(ncid,'kmax_t',kmax,k_t_dim)
      status=nfmpi_def_dim(ncid,'nk_t',kmaxp0,k_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 3'

      status=nfmpi_def_dim(ncid,'time',nfmpi_unlimited,time_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 4'

! dimensions grille z pour Kz
!      status=nfmpi_def_dim(ncid,'imax_w',iglbp2,i_w_dim)
      status=nfmpi_def_dim(ncid,'ni_w',iglbp2,i_w_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 5'

!      status=nfmpi_def_dim(ncid,'jmax_w',jglbp2,j_w_dim)
      status=nfmpi_def_dim(ncid,'nj_w',jglbp2,j_w_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 6'

!      status=nfmpi_def_dim(ncid,'kmax_w',kmaxp1,k_w_dim)
      status=nfmpi_def_dim(ncid,'nk_w',kmaxpp1,k_w_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 7'

! dimensions grille x pour u
!      status=nfmpi_def_dim(ncid,'imax_u',iglbp1,i_u_dim)
      status=nfmpi_def_dim(ncid,'ni_u',iglbp1,i_u_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 8'

!      status=nfmpi_def_dim(ncid,'jmax_u',jglbp2,j_u_dim)
      status=nfmpi_def_dim(ncid,'nj_u',jglbp2,j_u_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 9'

!      status=nfmpi_def_dim(ncid,'kmax_u',kmax,k_u_dim)
      status=nfmpi_def_dim(ncid,'nk_u',kmaxp0,k_u_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 10'

! dimensions grille y pour v
!      status=nfmpi_def_dim(ncid,'imax_v',iglbp2,i_v_dim)
      status=nfmpi_def_dim(ncid,'ni_v',iglbp2,i_v_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 11'

!      status=nfmpi_def_dim(ncid,'jmax_v',jglbp1,j_v_dim)
      status=nfmpi_def_dim(ncid,'nj_v',jglbp1,j_v_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 12'

!      status=nfmpi_def_dim(ncid,'kmax_v',kmax,k_v_dim)
      status=nfmpi_def_dim(ncid,'nk_v',kmaxp0,k_v_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 13'

! dimensions grille y pour f
!      status=nfmpi_def_dim(ncid,'imax_f',iglbp1,i_f_dim)
      status=nfmpi_def_dim(ncid,'ni_f',iglbp1,i_f_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 14'

!      status=nfmpi_def_dim(ncid,'jmax_f',jglbp1,j_f_dim)
      status=nfmpi_def_dim(ncid,'nj_f',jglbp1,j_f_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 15'

!      status=nfmpi_def_dim(ncid,'kmax_f',kmax  ,k_f_dim)
      status=nfmpi_def_dim(ncid,'nk_f',kmaxp0  ,k_f_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 16'


      end subroutine netcdf_dim


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine netcdf_error(ki1,ki2)
      use module_principal
      implicit none
      integer ki1,ki2
#ifdef synopsis
       subroutinetitle='netcdf_error'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       write(6,*)'status error= ',status !05-07-18

! diag netcdf errors in defining netcdf variables: begin
      if(ki1.eq.1) then

        write(6,'(4a)')'oups! Failed define netcdf :',trim(texte80(1)) &
                                                      ,' ',texte80(3)  !21-03-17
        if(ki2.eq.1)write(6,'(4a)')'in creating: varid. Suggestion:' & !21-03-17
      ,' the variable name ',trim(texte80(1)),' possibly already exists'
        if(ki2.eq.2)write(6,*)'in creating: units'
        if(ki2.eq.3)write(6,*)'in creating: long_name'
        if(ki2.eq.4)write(6,*)'in creating: standard name'
        if(ki2.eq.5)write(6,*)'in creating: axis'
        if(ki2.eq.6)write(6,*)'in creating: associate'
        if(ki2.eq.8)write(6,*)'in creating: time origin'
        if(ki2.eq.9)write(6,*)'in creating: calendar'
        if(ki2.eq.10)write(6,*)'in creating: missing_value'
        if(ki2.eq.11)write(6,*)'in creating: _FillValue'
        if(ki2.eq.12)write(6,*)'in creating: coordinates '
        if(ki2.eq.13)write(6,*)'in creating: C_grid_Oi_axis_shift'
        if(ki2.eq.14)write(6,*)'in creating: C_grid_Oj_axis_shift'
        if(ki2.eq.15)write(6,*)'in creating: C_grid_Ok_axis_shift'
        if(ki2.eq.16)write(6,*)'in creating: positive'
        if(ki2.eq.17)write(6,*)'in creating: scale_factor'
        if(ki2.eq.18)write(6,*)'in creating: add_offset'
        if(ki2.eq.19)write(6,*)'in creating: valid_min'
        if(ki2.eq.20)write(6,*)'in creating: valid_max'
        if(ki2.eq.21)write(6,*)'in creating: comment'


       stop ' stop in netcdf_error'
      endif
!!$! diag netcdf errors in defining netcdf variables: end



      end subroutine netcdf_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.............................................................

      subroutine netcdf_obc(txt2_   )
      use module_principal
      use module_parallele
      implicit none
      character*2,intent(in) :: txt2_
#ifdef synopsis
       subroutinetitle='netcdf_obc'
       subroutinedescription='Lateral Boundary Values in netcdf files'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Indicate now the variables that will be masked on the i=0 ..... boundaries:

      k0=0  ! Pas de masque par defaut

! Variables masquees (filval) dans les "zones fantomes":
      if(txt2_=='_v')                      k0=1  !masque regle generale u
      if(txt2_=='_u')                      k0=1  !masque regle generale v
      if(texte80(1)=='difv')               k0=1  !masque
      if(texte80(1)=='log10difv')          k0=1  !masque
      if(texte80(1)(1:3)=='tem')           k0=1  !masque
      if(texte80(1)(1:3)=='sal')           k0=1  !masque
!     if(texte80(1)(1:2)=='Ha')            k0=1  !masque !06-02-15
!     if(texte80(1)(1:2)=='Hg')            k0=1  !masque
      if(texte80(1)(1:3)=='bio')           k0=1  !masque
      if(texte80(1)(1:2)=='kh')            k0=1  !masque !16-01-15
      if(texte80(1)(1:5)=='err_h')         k0=1  !masque !07-11-15
      if(texte80(1)(1:13)=='delta_h_hogcm')k0=1  !masque !19-11-15
      if(index(texte80(1),'delta_t')/=0)   k0=1  !masque !22-12-15

      if(texte80(11)(1:11)=='mask_obc_z1')k0=1 
         texte80(11)='none' !masque puis reset!24-01-15

! Variables faisant exception aux regles generales u et v
      if(texte80(1)(1:7)=='depth_u')k0=0       !pas de masque
      if(texte80(1)(1:7)=='depth_v')k0=0       !pas de masque
      if(texte80(1)(1:4)=='kmin')   k0=0       !Pas de masque
      if(texte80(1)(1:10)=='latitude_u')k0=0   !pas de masque !13-04-16
      if(texte80(1)(1:10)=='latitude_v')k0=0   !pas de masque !13-04-16
      if(texte80(1)(1:10)=='latitude_f')k0=0   !pas de masque !13-04-16
      if(texte80(1)(1:11)=='longitude_u')k0=0  !pas de masque !13-04-16
      if(texte80(1)(1:11)=='longitude_v')k0=0  !pas de masque !13-04-16
      if(texte80(1)(1:11)=='longitude_f')k0=0  !pas de masque !13-04-16


      if(k0==0)return ! Si pas de masque on sort

       if(par%timax(1)==0) then !----->
         if(kend==1) then
          do j=jstr,jend
           anyvar2d(0,j)=filval
             anyv3d(0,j,1,1)=filval !07-11-15
          enddo
         else
          do k=1,kend
          do j=jstr,jend
           anyvar3d(0,j,k)=filval
             anyv3d(0,j,k,1)=filval !07-11-15
          enddo
          enddo
         endif
       endif                                      !----->

       if(imax+par%timax(1)==iglb) then !----->
       if(txt2_ /='_u') then !?> !01-11-13
         if(kend==1) then
          do j=jstr,jend
           anyvar2d(imax+1,j)=filval
             anyv3d(imax+1,j,1,1)=filval !07-11-15
          enddo
         else
          do k=1,kend
          do j=jstr,jend
           anyvar3d(imax+1,j,k)=filval
             anyv3d(imax+1,j,k,1)=filval !07-11-15
          enddo
          enddo
         endif
       endif                  !?>
       endif                                      !----->

       if(par%tjmax(1)==0) then   !----->
         if(kend==1) then
          do i=istr,iend
           anyvar2d(i,0)=filval
             anyv3d(i,0,1,1)=filval !07-11-15
          enddo
         else
          do k=1,kend
          do i=istr,iend
           anyvar3d(i,0,k)=filval
             anyv3d(i,0,k,1)=filval !07-11-15
          enddo
          enddo
         endif
       endif                                      !----->

       if(jmax+par%tjmax(1)==jglb) then   !----->
       if(txt2_/='_v') then  !?> !01-11-13
         if(kend==1) then
          do i=istr,iend
           anyvar2d(i,jmax+1)=filval
             anyv3d(i,jmax+1,1,1)=filval !07-11-15
          enddo
         else
          do k=1,kend
          do i=istr,iend
           anyvar3d(i,jmax+1,k)=filval
             anyv3d(i,jmax+1,k,1)=filval !07-11-15
          enddo
          enddo
         endif
       endif                  !?>
       endif                                      !----->


      end subroutine netcdf_obc
!.....................................................................

      subroutine netcdf_general_attributes(ncid_)
      use module_parallele
      use module_principal
      use pnetcdf
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) len_
      double precision :: scalarvalr8(1)
      integer :: scalarvalint(1)
#ifdef synopsis
       subroutinetitle='netcdf_general_attributes'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_nemoffline==0) then !pmxpxm> !26-01-17


      len_=len(trim(model_name))
      status=nfmpi_put_att_text(ncid_,nf_global,'Model',len_,model_name) !02-06-11

      len_=len('https://sites.google.com/view/symphonieoceanmodel/home')
      status=nfmpi_put_att_text(ncid_,nf_global,'Description',len_,&
      'https://sites.google.com/view/symphonieoceanmodel/home') !22-07-18


      len_=len(trim(title_for_netcdf_files_txt))
      status=nfmpi_put_att_text(ncid_,nf_global,'Title',len_,title_for_netcdf_files_txt) !03-04-12

      scalarvalr8(1)=dti_fw ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global,'Internal_time_step',nf_double,len_,scalarvalr8)
      scalarvalr8(1)=dte_fw ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global,'External_time_step',nf_double,len_,scalarvalr8)

! ----- Schema de fermeture turbulente
      if(iturbulence==0)texte90='Gaspar'
      if(iturbulence==1)texte90='K-epsilon'
      if(iturbulence==2)texte90='Constant Kz'
      len_=len(trim(texte90))
      status=nfmpi_put_att_text(ncid_,nf_global                    &
                                     ,'Turbulence_closure_scheme ' & !24-01-17
                                     ,len_,texte90)

       if(iturbulence==1) then !>>>
        len_=len(trim(flag_tke_stab_func))
         status=nfmpi_put_att_text(ncid_,nf_global          &
                                     ,'Stability_functions' & !24-01-17
                                     ,len_,flag_tke_stab_func)
       endif                   !>>>

! ----- Parametrage flux de surface
       if(bulk_scheme==bulk_core) texte90='Core' 
       if(bulk_scheme==bulk_moon) texte90='Moon'
       if(bulk_scheme==bulk_coare)texte90='Coare'
       len_=len(trim(texte90))
       status=nfmpi_put_att_text(ncid_,nf_global      &
                                     ,'Bulk_scheme '  & !24-01-17
                                     ,len_,texte90)

      scalarvalr8(1)=relativewind ; len_=1
      status=nfmpi_put_att_double(ncid_,nf_global,'Relative_Wind_Coef ',nf_double,len_,scalarvalr8)

       if(flag_steric_effect==1) then !>>>
        texte90='used'
        len_=len(trim(texte90))
        status=nfmpi_put_att_text(ncid_,nf_global      &
                                     ,'Steric effect parameterization'  & !28-04-18
                                     ,len_,texte90)
       endif                          !>>>

! ----- Maree:
      if(kmaxtide>0) then !pmx> !28-04-18
       scalarvalint(1)=kmaxtide ; len_=1  !22-05-17
       status=nfmpi_put_att_int(ncid_,nf_global                     &
                               ,'Number_of_tidal_components',nf_int &
                               ,len_,scalarvalint)
      endif               !pmx> !28-04-18

      endif                       !pmxpxm> !26-01-17

      if(flag_nemoffline==1) then !bloom> !26-01-17
       texte90='SIROCCO BLOOM ' ; len_=len(trim(texte90))
       status=nfmpi_put_att_text(ncid_,nf_global         &
                                     ,'NEMO Grid cloning ' & !26-01-17
                                     ,len_,texte90)
      endif                       !bloom> !26-01-17

      if(flag_nemoffline/=0.and.flag_nemoffline/=1) &
      stop 'netcdf flag_nemoffline/=0.and.flag_nemoffline/=1'

      if(flag_merged_levels==1) then !05-05-17
       texte90='s-z'
      else   
       texte90='s'
      endif
      len_=len(trim(texte90))
      status=nfmpi_put_att_text(ncid_,nf_global              &
                                    ,'Vertical_coordinate '  & !05-05-17
                                    ,len_,texte90)

      texte90='For any use of this file in a publication, ' &
            //'please inform the authors of the simulation.'
      len_=len(trim(texte90))
      status=nfmpi_put_att_text(ncid_,nf_global                   &
                                     ,'Terms_of_use' &
                                     ,len_,texte90)

      texte90='anonymous'
      if(par%rank==0) then !pmx> !29-03-18
      open(unit=3,file='authors_of_the_simulation') !13-01-17
         k1=1
   11    read(3,'(a)',end=19)texte90(k1:) ; k1=len(trim(texte90))+2
         goto 11
   19 close(3)
      endif                !pmx>
      call mpi_bcast(texte90,200,mpi_char,0,par%comm2d,ierr) 

      len_=len(trim(texte90))
      status=nfmpi_put_att_text(ncid_,nf_global                   &
                                     ,'Authors_of_the_simulation' &
                                     ,len_,texte90)

! Ecrire les versions de fortran, netcdf, mpi utilisees par S pour faire ce run 
! si le script de lancement contient la commande: module list 2> currently_loaded_modulefile
! Noter que le fichier peut contenir plusieurs lignes ce qui explique le protocole de
! lecture un peu etrange A premiere vue...
      texte90='?'
      if(par%rank==0) then !pmx> !22-06-18
       open(unit=3,file='currently_loaded_modulefiles')
         k1=1
  1341   read(3,'(a)',end=1338)texte90(k1:) ; k1=len(trim(texte90))+2
         goto 1341
  1338 close(3)
      endif                !pmx> !22-06-18
      call mpi_bcast(texte90,200,mpi_char,0,par%comm2d,ierr) 
      len_=len(trim(texte90))
      status=nfmpi_put_att_text(ncid_,nf_global &
                               ,'currently_loaded_modulefiles',len_,texte90)

      end subroutine netcdf_general_attributes
