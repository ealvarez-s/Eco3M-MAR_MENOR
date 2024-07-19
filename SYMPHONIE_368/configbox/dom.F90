      program dom
      implicit none
      double precision sum1,sum2,jmax_r8,imax_r8,dxb,dyb,phi0,longi0 &
                      ,latit0,rayonterre,northpole_lon               &
                      ,northpole_lat,southpole_lon,southpole_lat     &
                      ,factor
      integer :: nbdom_imax=0    &
                ,nbdom_jmax=0
      integer i,j,rank,k,i1,i2,j1,j2,i0,j0,kmax,grid_i0,grid_j0
      logical iperiodicboundary,jperiodicboundary
      integer,dimension(:),allocatable :: power2
      integer,dimension(:,:),allocatable :: &
       rank_ij,imax_ij,jmax_ij,subcycle_number
      integer ::           &
                 iglb=0    &
                ,jglb=0    &
                ,iperiodic=0 & ! Par defaut 0, 1 si periodique en i
                ,jperiodic=0 & ! Par defaut 0, 1 si periodique en j
                ,debug=0     &
                ,delta_i,delta_j,sum_i,sum_j
      character directory*100,notebook_name*100,initialgridfile_txt*100 &
               ,vert_axis_conv_direc*2,vert_axis_conv_start*1          &
               ,vert_axis_conv_end*1,hori_axis_conv_start*1            &
               ,hori_axis_conv_end*1                                   &
               ,lonlatfile*100
      logical fplan_grid

      namelist/notebook_grid/iglb,jglb,kmax,nbdom_imax,nbdom_jmax   &
      ,iperiodicboundary,jperiodicboundary,dxb,dyb,phi0,longi0       &
      ,latit0,grid_i0,grid_j0,rayonterre,northpole_lon               &
      ,northpole_lat,initialgridfile_txt,vert_axis_conv_direc        &
      ,vert_axis_conv_start,vert_axis_conv_end,hori_axis_conv_start  &
      ,southpole_lon,southpole_lat,lonlatfile,hori_axis_conv_end     &
      ,fplan_grid

      open(unit=66,file='description_domaine.txt')

      open(unit=3,file='notebook_list')
       read(3,'(40x,a)')directory
       read(3,*)
       read(3,*)
       read(3,'(40x,a)')notebook_name
      close(3)
      notebook_name=trim(directory)//'/'//trim(notebook_name)
      write(6,'(a)')trim(directory)
      write(6,'(a)')trim(notebook_name)

! Lecture du notebook_grid
      open(100,file=trim(notebook_name))
      read(100,nml=notebook_grid)
      close(100)

      if(iperiodicboundary)iperiodic=1
      if(jperiodicboundary)jperiodic=1

      allocate(rank_ij        (-1:nbdom_imax  ,-1:nbdom_jmax ))
      allocate(imax_ij        ( 0:nbdom_imax-1,0:nbdom_jmax-1))
      allocate(jmax_ij        ( 0:nbdom_imax-1,0:nbdom_jmax-1))
      allocate(subcycle_number( 0:nbdom_imax-1,0:nbdom_jmax-1)) 
      allocate(power2         ( 0:nbdom_jmax-1 ))
      rank_ij=-2

! RESET
! Par défaut pas de subcycling soit subcycle_number=1
      do j=0,nbdom_jmax-1
      do i=0,nbdom_imax-1
            subcycle_number(i,j)=1
      enddo
      enddo
      
      do j=0,nbdom_jmax-1

! Definir le nombre de sous-cycle à l'interieur du plus grand pas de temps
! en fonction de l'indice j de la numerotation des sous-domaines.
! Attention à respecter les valeurs possibles: 1 2 4 8 16 32 etc.....

#ifdef bidon
      if(j==0)subcycle_number(0:nbdom_imax-1,j)=1
      if(j==1)subcycle_number(0:nbdom_imax-1,j)=2
      if(j==2)subcycle_number(0:nbdom_imax-1,j)=2
      if(j==3)subcycle_number(0:nbdom_imax-1,j)=2
      if(j==4)subcycle_number(0:nbdom_imax-1,j)=4
      if(j==5)subcycle_number(0:nbdom_imax-1,j)=4
      if(j==6)subcycle_number(0:nbdom_imax-1,j)=8
      if(j==7)subcycle_number(0:nbdom_imax-1,j)=8
#endif

       debug=0
       do k=0,10
        if(subcycle_number(0,j)==2**k) then !---->
          debug=1
          power2(j)=k
        endif                               !---->
       enddo
       if(debug==0) then
        write(6,*)subcycle_number(0:nbdom_imax-1,j), &
        ' est une valeur impossible pour subcycle_number'
        stop' STOP ERREUR'
       endif


      enddo

      factor=1.
      sum1=0.
      do j=0,nbdom_jmax-1
!     sum1=sum1+1./subcycle_number(0,j)
      sum1=sum1+(factor**power2(j))/subcycle_number(0,j)
      enddo

!     imax_r8=(iglb-2)/(nbdom_imax) + 2
      imax_r8=real(iglb-2)/real(nbdom_imax) + 2.
      jmax_r8=(jglb+2*(nbdom_jmax-1))/sum1

      write(6,*)'jmax_r8=',jmax_r8
      write(6,*)'imax_r8=',imax_r8

      rank=0
      do i=0,nbdom_imax-1
      do j=0,nbdom_jmax-1

!      subcycle_number(i,j)=2**power(j)

       imax_ij(i,j)=int(imax_r8)
       jmax_ij(i,j)=int(jmax_r8)

!      jmax_ij(i,j)=nint(jmax_r8/(2**power(j)))
!!     jmax_ij(i,j)=nint(jmax_r8/subcycle_number(i,j))
!      jmax_ij(i,j)=nint( 2.   &
!              +(jglb-2)/sum1*(factor**power2(j))/subcycle_number(0,j) )

!      if(j==0)write(6,*)'i,j,jmax_ij(i,j)',i,j,jmax_ij(i,j)
!      if(j==0)write(6,*)'i,j,imax_ij(i,j)',i,j,imax_ij(i,j)

       rank_ij(i,j)=rank
       rank=rank+1

!      write(6,*)i,j,rank_ij(i,j)

      enddo
      enddo
! A ce stade tous les procs ont les memes dimensions mais
! a cause de l'arrondi à l'entier inferieur la somme sur i des imax_ij(i,j)-2*(nbdom_imax-1) n'est pas iglb
! Le deficit: iglb - [ nbdom_imax*(imax_ij(0,0)-2)+2 ] est reparti, "un a un" sur les premiers domaines:
      delta_i = iglb - ( nbdom_imax*(imax_ij(0,0)-2)+2 )
      delta_j = jglb - ( nbdom_jmax*(jmax_ij(0,0)-2)+2 )

      write(6,*)'delta_i=',delta_i
      write(6,*)'delta_j=',delta_j

      do i=0,delta_i-1
      do j=0,nbdom_jmax-1
       imax_ij(i,j)=imax_ij(i,j)+1
      enddo
      enddo

      do i=0,nbdom_imax-1
      do j=0,delta_j-1
       jmax_ij(i,j)=jmax_ij(i,j)+1
      enddo
      enddo
      
! Verification
      do j=0,nbdom_jmax-1
       sum_i=2
       do i=0,nbdom_imax-1
        sum_i=sum_i+imax_ij(i,j)-2
       enddo
       if(sum_i/=iglb)stop'erreur sum_i'
      enddo
      do i=0,nbdom_imax-1
       sum_j=2
       do j=0,nbdom_jmax-1
        sum_j=sum_j+jmax_ij(i,j)-2
       enddo
       if(sum_j/=jglb)stop'erreur sum_j'
      enddo

      write(6,*)'sum_i iglb = ',sum_i,iglb
      write(6,*)'sum_j jglb = ',sum_j,jglb

! conditions cycliques direction Oi
      if(iperiodic==1) then !>>>>>>>>
       do j=0,nbdom_jmax-1
! OBC: nbdom_imax
        rank_ij(nbdom_imax,j-1)=rank_ij(0,j-1)
        rank_ij(nbdom_imax,j  )=rank_ij(0,j  )
        rank_ij(nbdom_imax,j+1)=rank_ij(0,j+1)
! OBC: -1
        rank_ij(-1,j-1)=rank_ij(nbdom_imax-1,j-1)
        rank_ij(-1,j  )=rank_ij(nbdom_imax-1,j  )
        rank_ij(-1,j+1)=rank_ij(nbdom_imax-1,j+1)
       enddo
      endif                 !>>>>>>>>
      if(jperiodic==1) then !>>>>>>>>
       do i=0,nbdom_imax-1
! OBC: nbdom_jmax
        rank_ij(i-1,nbdom_jmax)=rank_ij(i-1,0)
        rank_ij(i  ,nbdom_jmax)=rank_ij(i  ,0)
        rank_ij(i+1,nbdom_jmax)=rank_ij(i+1,0)
! OBC: -1
        rank_ij(i-1,-1)=rank_ij(i-1,nbdom_jmax-1)
        rank_ij(i  ,-1)=rank_ij(i  ,nbdom_jmax-1)
        rank_ij(i+1,-1)=rank_ij(i+1,nbdom_jmax-1)
       enddo
      endif                 !>>>>>>>>

      write(6,*)'voisin West',rank_ij(-1,0),rank_ij(nbdom_imax-1,0  )


      write(66,'(3i6,a60)')nbdom_imax,nbdom_jmax  &
                          ,nbdom_imax*nbdom_jmax, &
       ' ! Number of sub-domains in each direction & nbdom'
      write(66,*)iglb,jglb,' ! iglb jglb'
      do i=0,nbdom_imax-1
      do j=0,nbdom_jmax-1

       i2=imax_ij(0,j)
       do i1=1,i
       i2=i2+imax_ij(i1,j)-2
       enddo

       j2=jmax_ij(i,0)
       do j1=1,j
       j2=j2+jmax_ij(i,j1)-2
       enddo
       if(j==nbdom_jmax-1)then
         j0=j2
         j2=jglb
         jmax_ij(i,j)=jmax_ij(i,j)-j0+j2
       endif
       if(i==nbdom_imax-1)then
         i0=i2
         i2=iglb
         imax_ij(i,j)=imax_ij(i,j)-i0+i2
       endif

!     write(6,*)'iend-istart=',imax_ij(i,j)
!     write(6,*)'jend-jstart=',jmax_ij(i,j)

      write(66,*)'!---------------------------'
      write(66,*)rank_ij(i,j),'             ! sub-domain order number'
      write(66,*)i,j,' ! sub-domain (i,j) indexes in the mpi space'
      write(66,*)subcycle_number(i,j) &
      ,'             ! number of cycles in one principal time step'
!     write(66,*)i2-imax_ij(i,j),i2,' ! i start i end'
!     write(66,*)j2-jmax_ij(i,j),j2,' ! j start j end'

      write(66,'(i5,1x,i5,1x,i5,a)')     &
       i2-imax_ij(i,j),i2,imax_ij(i,j)   &
      ,'         ! i start i end imax'

      write(66,'(i5,1x,i5,1x,i5,a)')     &
       j2-jmax_ij(i,j),j2,jmax_ij(i,j)   &
      ,'         ! j start j end jmax'

      write(66,'(8(i4,1x),a40)')        &
                rank_ij(i-1,j  )    &
               ,rank_ij(i+1,j  )    &
               ,rank_ij(i  ,j+1)    &
               ,rank_ij(i  ,j-1)    &
               ,rank_ij(i-1,j-1)    &
               ,rank_ij(i+1,j-1)    &
               ,rank_ij(i-1,j+1)    &
               ,rank_ij(i+1,j+1)    &
               ,' ! Neighbors: w e n s ws es wn en'
      i1=i2-2
      j1=j2-2
      enddo
      enddo


      end
!-------------------------------------
!--  Description
!-------------------------------------
! Nombre de domaines
! iglb, jglb
! Numero de domaine  :: Boucle sur le nombre de domaine
! Coordonee du domaine suivant i, suivant j
! imin, imax
! jmin, jmax
! voisin i inferieur
! voisin i superieur
! voisin j superieur 
! voisin j inferieur
! voisin (j inferieur & i inferieur) 
! voisin (j inferieur & i superieur)  
! voisin (j superieur & i inferieur) 
! voisin (j supeireur & i superieur)
!! Dit autrement :  !Ouest, Est, Nord, Sud,sudouest,sudest,nordouest,nordest
!!!! Boucle sur le nombre de domaine continue



