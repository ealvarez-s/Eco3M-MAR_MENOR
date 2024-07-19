       subroutine graph_out_sedim
#ifdef key_sedim
       use module_principal
       use module_parallele
       use module_parameter_sedim
       use comsubstance
       use sedimento, only : dzs, cv_sed, ksma
       use module_biology, only : fluxbio_w,vbmax_eco3ms
       use pnetcdf
       implicit none

!------------------------------------------------------------------!
! déclaration local
!------------------------------------------------------------------!
      integer :: iv
      integer(kind=MPI_OFFSET_KIND) :: iglbp0
      integer(kind=MPI_OFFSET_KIND) :: jglbp0
      integer(kind=MPI_OFFSET_KIND) :: kmaxp0
      integer(kind=MPI_OFFSET_KIND) :: len_
      integer(kind=MPI_OFFSET_KIND),dimension(4) :: start,edge
      integer,dimension(1:nv_tot) :: varsedid, varsedwid, varflxs2wid  &
                                   , varflxw2sid
      integer :: vardzsid
!     real(kind=rsh),dimension(imax,jmax,ksdmin:ksdmax,1) :: anyvarsed
!     real(kind=rsh),dimension(imax,jmax,1)               :: anyvarflxsed
!     real(kind=rsh),dimension(imax,jmax,kmax,1)          :: anyvarsedwater
      real(kind=rsh),dimension(:,:,:,:),allocatable :: anyvarsed
      real(kind=rsh),dimension(:,:,:)  ,allocatable :: anyvarflxsed
      real(kind=rsh),dimension(:,:,:,:),allocatable :: anyvarsedwater
!------------------------------------------------------------------!
! partie excutable 
!------------------------------------------------------------------!

      if(.not.allocated(anyvarsed))     allocate(anyvarsed(imax,jmax,ksdmin:ksdmax,1));anyvarsed=0.
      if(.not.allocated(anyvarflxsed))  allocate(anyvarflxsed(imax,jmax,1))           ;anyvarflxsed=0.
      if(.not.allocated(anyvarsedwater))allocate(anyvarsedwater(imax,jmax,kmax,1))    ;anyvarsedwater=0.

! définition du nom du fichier
      call elapsedtime2date(elapsedtime_now,i5,i6,i7,i3,i2,i1)     !01-04-13
      write(texte80(1),                            &
      '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')  &
      i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

!.....écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04


! Nom des fichiers de variables
      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //'.SEDIM-S.LA.tlse.nc'

!.....puis afficher la date à l'écran:                                     !27/01/03
      if(par%rank==0) then !00000>
       write(6,*)'--------------------------------'
       write(6,*)'SEDIMENTO netcdf file date='
       write(6,'(a33)')texte80(1)(1:33)
       write(6,'(a,a)')'File name: ',trim(texte250)
      endif                !00000>

! creation du fichier de sortie
      status=nfmpi_create(par%comm2d,texte250,nf_clobber,MPI_INFO_NULL,ncid)


      loop_netcdf=0

      iglbp0=iglb
      jglbp0=jglb
      kmaxp0=ksdmax

      istr=1
      jstr=1

      iend=imax
      jend=jmax

      start(1)=1+par%timax(1) ; start(2)=1+par%tjmax(1)
      start(3)=1
      start(4)=1 

      edge(1)=iend
      edge(2)=jend
      edge(3)=ksdmax
      edge(4)=1

      status=0
      status=nfmpi_def_dim(ncid,'ni_t',iglbp0,i_t_dim)
      if(status/=0) then
          write(6,*)'status=',status
          stop ' stop erreur netcdf_dim 1 sedim'
      endif

!      status=nfmpi_def_dim(ncid,'jmax_t',jglbp2,j_t_dim)
      status=nfmpi_def_dim(ncid,'nj_t',jglbp0,j_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 2 sedim'

!      status=nfmpi_def_dim(ncid,'kmax_t',kmax,k_t_dim)
      status=nfmpi_def_dim(ncid,'nk_t',kmaxp0,k_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 3 sedim'

      status=nfmpi_def_dim(ncid,'time',nfmpi_unlimited,time_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 4 sedim'


!Definition des axes
      k0=4 ; vardim(1)=i_t_dim ; vardim(2)=j_t_dim                ! 3D ; dim1 ; dim2
      vardim(3)=k_t_dim ; vardim(4)=time_dim

      do iv=vbmax_eco3ms+1,nv_tot 

        ! concentration dans le sediment
        write(texte80(1),'(a,i0)')'sedim',iv
        texte80(2)='none'         ! variable ; units
        write(texte80(3),'(a,i0)')'sediment_variable_number_',iv
        texte80(4)=texte80(3)
        texte80(7)='real'
        status=nfmpi_def_var(ncid,texte80(1),nf_real,k0,         &
                            (/i_t_dim,j_t_dim,k_t_dim,time_dim/) &
                           ,varsedid(iv))
        if (status/=0)stop 'erreur init var sedim'
        len_=len_trim(texte80(2))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'unit',len_,texte80(2))
        if (status/=0)stop 'erreur init var sedim unit'
        len_=len_trim(texte80(3))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'long_name',len_,texte80(3))
        if (status/=0)stop 'erreur init var sedim long_name'
        len_=len_trim(texte80(4))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'standard_name',len_,texte80(4))
        if (status/=0)stop 'erreur init var sedim standrad_name'

        ! flux du sediment vers la colonne d'eau
        write(texte80(1),'(a,i0)')'flx_s2w_',iv
        texte80(2)='none'         ! variable ; units
        write(texte80(3),'(a,i0)')'flx_s2w_variable_number_',iv
        texte80(4)=texte80(3)
        texte80(7)='real'
        status=nfmpi_def_var(ncid,texte80(1),nf_real,k0-1,         &
                            (/i_t_dim,j_t_dim,time_dim/) &
                           ,varflxs2wid(iv))
        if (status/=0)stop 'erreur init var flxs2w'
        len_=len_trim(texte80(2))
        status=nfmpi_put_att_text(ncid,varflxs2wid(iv),'unit',len_,texte80(2))
        if (status/=0)stop 'erreur init var flxs2w unit'
        len_=len_trim(texte80(3))
        status=nfmpi_put_att_text(ncid,varflxs2wid(iv),'long_name',len_,texte80(3))
        if (status/=0)stop 'erreur init var flxs2w long_name'
        len_=len_trim(texte80(4))
        status=nfmpi_put_att_text(ncid,varflxs2wid(iv),'standard_name',len_,texte80(4))
        if (status/=0)stop 'erreur init var flxs2w standrad_name'

        ! flux de la colonne d'eau vers le sediment
        write(texte80(1),'(a,i0)')'flx_w2s_',iv
        texte80(2)='none'         ! variable ; units
        write(texte80(3),'(a,i0)')'flx_w2s_variable_number_',iv
        texte80(4)=texte80(3)
        texte80(7)='real'
        status=nfmpi_def_var(ncid,texte80(1),nf_real,k0-1,         &
                            (/i_t_dim,j_t_dim,time_dim/) &
                           ,varflxw2sid(iv))
        if (status/=0)stop 'erreur init var flxw2s'
        
        len_=len_trim(texte80(2))
        status=nfmpi_put_att_text(ncid,varflxw2sid(iv),'unit',len_,texte80(2))
        if (status/=0)stop 'erreur init var flxw2s unit'
        len_=len_trim(texte80(3))
        status=nfmpi_put_att_text(ncid,varflxw2sid(iv),'long_name',len_,texte80(3))
        if (status/=0)stop 'erreur init var flxw2s long_name'
        len_=len_trim(texte80(4))
        status=nfmpi_put_att_text(ncid,varflxw2sid(iv),'standard_name',len_,texte80(4))
        if (status/=0)stop 'erreur init var flxw2s standrad_name'

      enddo
      

      write(texte80(1),'(a,i0)')'dzs'
      texte80(2)='m'         ! variable ; units
      write(texte80(3),'(a,i0)')'epaisseur_de_couche'
      texte80(4)=texte80(3)
      texte80(5)='TZYX' ; texte80(7)='real'
      status=nfmpi_def_var(ncid,texte80(1),nf_real,k0,         &
                          (/i_t_dim,j_t_dim,k_t_dim,time_dim/) &
                         ,vardzsid)
      if (status/=0)stop 'erreur init var dzs'
        len_=len_trim(texte80(2))
      status=nfmpi_put_att_text(ncid,vardzsid,'unit',len_,texte80(2))
      if (status/=0)stop 'erreur init var dzs unit'
        len_=len_trim(texte80(3))
      status=nfmpi_put_att_text(ncid,vardzsid,'long_name',len_,texte80(3))
      if (status/=0)stop 'erreur init var dzs long_name'
        len_=len_trim(texte80(4))
      status=nfmpi_put_att_text(ncid,vardzsid,'standard_name',len_,texte80(4))
      if (status/=0)stop 'erreur init var dzs standrad_name'

! fin des déclarations
      status=nfmpi_enddef(ncid)
      if(status/=0)stop 'erreur fin definition netcdf sedimento'

! phase d'ecriture
      filval=-9999.
      filval=0.
      do iv=vbmax_eco3ms+1,nv_tot

        ! concentration dans le sédiment
        anyvarsed(:,:,:,1)=filval 
        do j=1,jmax ; do i=1,imax ; do k1=1,ksma(i,j)
         anyvarsed(i,j,k1,1)=cv_sed(iv,k1,i,j)*mask_t(i,j,kmax) &
                          +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo ; enddo

        status=nfmpi_put_vara_real_all(ncid,varsedid(iv),start(1:4),edge(1:4),anyvarsed(1:imax,1:jmax,1:ksdmax,1))
        if(status/=0) stop 'erreur ecriture var sedim'

        anyvarsed(:,:,:,1)=filval 
        do j=1,jmax ; do i=1,imax ; do k1=1,ksma(i,j)
        anyvarsed(i,j,k1,1)=dzs(k1,i,j)*mask_t(i,j,kmax)          & !20-12-14
                        +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo ; enddo

        status=nfmpi_put_vara_real_all(ncid,vardzsid,start(1:4),edge(1:4),anyvarsed(1:imax,1:jmax,ksdmin:ksdmax,1))
        if(status/=0) stop 'erreur ecriture var dzs'

        ! flux du sediment vers la colonne d'eau
        anyvarflxsed(:,:,1)=filval 
        do j=1,jmax ; do i=1,imax 
         anyvarflxsed(i,j,1)=fluxbio_w(i,j,iv,1)*mask_t(i,j,kmax) &
                          +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo 

        status=nfmpi_put_vara_real_all(ncid,varflxs2wid(iv),(/start(1),start(2),start(4)/),(/edge(1),edge(2),edge(4)/),anyvarflxsed(1:imax,1:jmax,1))
        if(status/=0) stop 'erreur ecriture var fluxbio_w'

        ! flux du sediment vers la colonne d'eau
        anyvarflxsed(:,:,1)=filval 
        do j=1,jmax ; do i=1,imax 
         anyvarflxsed(i,j,1)=flx_w2s(iv,i,j)*mask_t(i,j,kmax) &
                          +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo 

        status=nfmpi_put_vara_real_all(ncid,varflxw2sid(iv),(/start(1),start(2),start(4)/),(/edge(1),edge(2),edge(4)/),anyvarflxsed(1:imax,1:jmax,1))
        if(status/=0) stop 'erreur ecriture var flx_w2s'

      enddo


! fermeture du fichier
      status=nfmpi_close(ncid)
      if(status/=0)stop 'erreur fermeture fichier sedimento water'

! fichier des variables dans la colonne d'eau
! Nom des fichiers de variables
      texte250=dirgraph(1:lname4)//texte90(1:8)              &
                            //'_'//texte90(10:15)            &
                                 //'.WATER-S.LA.tlse.nc'

! creation du fichier de sortie
      status=nfmpi_create(par%comm2d,texte250,nf_clobber,MPI_INFO_NULL,ncid)


      loop_netcdf=0

      iglbp0=iglb
      jglbp0=jglb
      kmaxp0=kmax

      istr=1
      jstr=1

      iend=imax
      jend=jmax

      start(1)=1+par%timax(1)
      start(2)=1+par%tjmax(1)
      start(3)=1
      start(4)=1

      edge(1)=iend
      edge(2)=jend
      edge(3)=kmax
      edge(4)=1

      status=0
      status=nfmpi_def_dim(ncid,'ni_t',iglbp0,i_t_dim)
      if(status/=0) then
          write(6,*)'status=',status
          stop ' stop erreur netcdf_dim 1 water'
      endif

!      status=nfmpi_def_dim(ncid,'jmax_t',jglbp2,j_t_dim)
      status=nfmpi_def_dim(ncid,'nj_t',jglbp0,j_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 2 water'

!      status=nfmpi_def_dim(ncid,'kmax_t',kmax,k_t_dim)
      status=nfmpi_def_dim(ncid,'nk_t',kmaxp0,k_t_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 3 water'

      status=nfmpi_def_dim(ncid,'time',nfmpi_unlimited,time_dim)
      if(status/=0)stop ' stop erreur netcdf_dim 4 water'

!Definition des axes
      k0=4 ; vardim(1)=i_t_dim ; vardim(2)=j_t_dim                ! 3D ; dim1 ; dim2
      vardim(3)=k_t_dim ; vardim(4)=time_dim

      do iv=vbmax_eco3ms+1,nv_tot 

        write(texte80(1),'(a,i0)')'sedim',iv
        texte80(2)='none'         ! variable ; units
        write(texte80(3),'(a,i0)')'sediment_variable_number_',iv
        texte80(4)=texte80(3)
        texte80(7)='real'
        status=nfmpi_def_var(ncid,texte80(1),nf_real,k0,         &
                            (/i_t_dim,j_t_dim,k_t_dim,time_dim/) &
                           ,varsedid(iv))
        if (status/=0)stop 'erreur init var sedim water'
        len_=len_trim(texte80(2))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'unit',len_,texte80(2))
        if (status/=0)stop 'erreur init var sedim unit water'
        len_=len_trim(texte80(3))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'long_name',len_,texte80(3))
        if (status/=0)stop 'erreur init var sedim long_name water'
        len_=len_trim(texte80(4))
        status=nfmpi_put_att_text(ncid,varsedid(iv),'standard_name',len_,texte80(4))
        if (status/=0)stop 'erreur init var sedim standrad_name water'
        
        write(texte80(1),'(a,i0)')'sedim_w',iv
        write(texte80(3),'(a,i0)')'sediment_variable_number_',iv
        status=nfmpi_def_var(ncid,texte80(1),nf_real,k0,         &
                            (/i_t_dim,j_t_dim,k_t_dim,time_dim/) &
                           ,varsedwid(iv))
        if (status/=0)stop 'erreur init var sedim water'
        len_=len_trim(texte80(2))
        status=nfmpi_put_att_text(ncid,varsedwid(iv),'unit',len_,texte80(2))
        if (status/=0)stop 'erreur init var sedim unit water'
        len_=len_trim(texte80(3))
        status=nfmpi_put_att_text(ncid,varsedwid(iv),'long_name',len_,texte80(3))
        if (status/=0)stop 'erreur init var sedim long_name water'
        len_=len_trim(texte80(4))
        status=nfmpi_put_att_text(ncid,varsedwid(iv),'standard_name',len_,texte80(4))
        if (status/=0)stop 'erreur init var sedim standrad_name water'
        

        
      enddo
      
! fin des déclarations
      status=nfmpi_enddef(ncid)
      if(status/=0)stop 'erreur fin definition netcdf sedimento water'

! phase d'ecriture
      filval=-9999.
      filval=0.
      do iv=vbmax_eco3ms+1,nv_tot

        anyvarsed(:,:,:,1)=filval 
        do j=1,jmax ; do i=1,imax ; do k1=1,kmax
         anyvarsedwater(i,j,k1,1)=bio_t(i,j,k1,iv)*mask_t(i,j,kmax)   & 
                          +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo ; enddo

        status=nfmpi_put_vara_real_all(ncid,varsedid(iv),start(1:4),edge(1:4),anyvarsedwater(1:imax,1:jmax,1:kmax,1))
        if(status/=0) stop 'erreur ecriture var sedim water'

        anyvarsed(:,:,:,1)=filval 
        do j=1,jmax ; do i=1,imax ; do k1=1,kmax
         anyvarsedwater(i,j,k1,1)=ws3(k1,iv,i,j)*mask_t(i,j,kmax)   &
                          +(1-mask_t(i,j,kmax))*filval
        enddo ; enddo ; enddo

        status=nfmpi_put_vara_real_all(ncid,varsedwid(iv),start(1:4),edge(1:4),anyvarsedwater(1:imax,1:jmax,1:kmax,1))
        if(status/=0) stop 'erreur ecriture var sedim water'

      enddo

! fermeture du fichier
      status=nfmpi_close(ncid)
      if(status/=0)stop 'erreur fermeture fichier sedimento water'


#endif
       end subroutine graph_out_sedim
