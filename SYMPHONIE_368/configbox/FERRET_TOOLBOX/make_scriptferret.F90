      program make_scriptferret
      implicit none
      character*200 varfile,gridfile,txtline,txtline2
      integer :: varmax=10            &
                ,axisnum,normalindex  &
                ,alongindex1          &
                ,alongindex2          &
                ,loop=0               &
                ,loopmax              &
                ,rangeoption          &
                ,axisunit             &
                ,kount=0              &
                ,plannum              &
                ,klevel
      integer, dimension(:) , allocatable :: varnum, cgridloc
      character*30, dimension(:) , allocatable :: varname
      character*60 command,txtkount
      double precision rangevalue(3)

      allocate(varnum(varmax))   ; varnum=0
      allocate(varname(varmax))  ; varname='none'
      allocate(cgridloc(varmax)) ; cgridloc=1 ! _t=1, _u=2, _v=3
!     allocate(axisunit(varmax)) ; axisunit=1 !  indexes, latitude, longitude


      write(6,*)'Enter grid file: (0=default name: grid.nc)'
      read(5,'(a)')gridfile
      if(trim(gridfile)=='0')gridfile='grid.nc'
      write(6,'(a)')trim(gridfile)

      write(6,*)'Enter variable file: (0=all 2....nc files)'
      read(5,'(a)')varfile

   
      if(trim(varfile)=='0')then
       command='ls 2*nc > file_list'
      else
       command='ls '//trim(varfile)//' > file_list'
      endif
      call system(command)

!     goto 1000

   23 loop=loop+1
      write(6,*)'2DV script processor'
      write(6,*)'Enter variable number:'
      write(6,*)'temperature.. 1'      
      write(6,*)'temref....... 2'      
      write(6,*)'temref-temobc 5'      
      write(6,*)'tem-temobc   16'      
      write(6,*)'salinity..... 3'      
      write(6,*)'salref....... 4'      
      write(6,*)'salref-salobc 6'      
      write(6,*)'sal-salobc   17'      
      write(6,*)'presgrad_u....7'      
      write(6,*)'presgrad_v....8'      
      write(6,*)'u         ....9'      
      write(6,*)'v        ....10'      
      write(6,*)'velobc_u.....11'      
      write(6,*)'velobc_v.....12'      
      write(6,*)'temobc_t.....13'      
      write(6,*)'salobc_t.....14'      
      write(6,*)'KH      .....18'      
      write(6,*)'ssh_w   .....15'      
      write(6,*)'continue..... 0'      
      read(5,*)varnum(loop)

      if(varnum(loop)==1) then
        varname(loop)='tem' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==2) then
        varname(loop)='temref' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==3) then
        varname(loop)='sal' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==4) then
        varname(loop)='salref' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==5)then
        varname(loop)='temref_temobc' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==6) then
        varname(loop)='salref_salobc' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==7) then
         varname(loop)='presgrad_u' ; cgridloc(loop)=2
      endif
      if(varnum(loop)==8) then
         varname(loop)='presgrad_v' ; cgridloc(loop)=3
      endif
      if(varnum(loop)==9) then
         varname(loop)='u' ; cgridloc(loop)=2
      endif
      if(varnum(loop)==10) then
         varname(loop)='v' ; cgridloc(loop)=3
      endif
      if(varnum(loop)==11) then
         varname(loop)='velobc_u' ; cgridloc(loop)=2
      endif
      if(varnum(loop)==12) then
         varname(loop)='velobc_v' ; cgridloc(loop)=3
      endif
      if(varnum(loop)==13) then
         varname(loop)='temobc_t' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==14) then
         varname(loop)='salobc_t' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==15) then
         varname(loop)='ssh_w' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==16)then
        varname(loop)='tem_temobc' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==17)then
        varname(loop)='sal_salobc' ; cgridloc(loop)=1
      endif
      if(varnum(loop)==18)then
        varname(loop)='kh' ; cgridloc(loop)=0
      endif

      if(varnum(loop)/=0)goto 23
      loopmax=loop-1

      write(6,*)'Define plan'
      write(6,*)'2D horizontal...1'
      write(6,*)'2D vertical.....2'
      read(5,*)plannum

      if(plannum==2) then !>>>>>>
       write(6,*)'Define horizontal axis'
       write(6,*)'Oi axis .......1'
       write(6,*)'Oj axis .......2'
       read(5,*)axisnum
       write(6,*)'Define horizontal axis unity'
       write(6,*)'grid points.1'
       write(6,*)'latitude....2'
       write(6,*)'longitude...3'
       read(5,*)axisunit
       write(6,*)'Enter index on normal axis:'
       read(5,*)normalindex
       write(6,*)'Enter min & max index on alongaxis (0,0 otherwise)'
       read(5,*)alongindex1,alongindex2
      else                !>>>>>>
       axisnum=0
       axisunit=0
       normalindex=0
       alongindex1=0 ; alongindex2=0
       write(6,*)'Enter vertical level (or 0)'
       read(5,*)klevel
      endif               !>>>>>>

      write(6,*)'range option? yes=1 no=0'
      read(5,*)rangeoption
      if(rangeoption/=1)goto 1000
      write(6,*)'Enter min max step:'
      read(5,*)rangevalue(1:3)

 1000 continue

      open(unit=10,file='job.jnl')
      open(unit=11,file='file_list')

!     txtline='use "../OFFLINE/grid.nc"'
      txtline='use "'//trim(gridfile)//'"'
      write(10,'(a)')trim(txtline)
  
  201 read(11,*,end=120)varfile
      kount=kount+1

      txtline='use "'//trim(varfile)//'"'
      write(10,'(a)')trim(txtline)
      
      txtline='show d' 
     !write(10,'(a)')trim(txtline)

      do loop=1,loopmax !---------->

       if(axisnum==0) then !mmmm>
         if(klevel/=0)write(txtline,'(a,i0)')'shade/k=',klevel
         if(klevel==0)write(txtline,'(a   )')'shade '
       endif               !mmmm>

       if(axisnum==1) then
         write(txtline,'(a,i0)')'shade/j=',normalindex
         if(alongindex1/=0)write(txtline,'(a,a,i0,a,i0)') &
         trim(txtline),'/i=',alongindex1,':',alongindex2
       endif
       if(axisnum==2) then
        write(txtline,'(a,i0)')'shade/i=',normalindex
         if(alongindex1/=0)write(txtline,'(a,a,i0,a,i0)') &
         trim(txtline),'/j=',alongindex1,':',alongindex2
       endif

       if(rangeoption==1) then !>>>
       write(txtline2,'(3(e13.6,1x))')real(rangevalue)
       txtline2=adjustl(txtline2)

       write(txtline,'(a,a,a,a)')trim(txtline),'/levels=(' &
                          ,trim(txtline2),')'
       endif                   !>>>
!/levels=(-0.5 0.5 0.01)

       write(txtline,'(a,a,a,a)')trim(txtline),' ' &
            ,trim(varname(loop)),'[d=2]'

       if(axisnum/=0) then !0000000000>

       if(axisnum==1) then !1111111111>
        if(axisunit==1) then
         if(cgridloc(loop)==0)write(txtline2,'(a)')', i_index_t[d=1]'
         if(cgridloc(loop)==1)write(txtline2,'(a)')', i_index_t[d=1]'
         if(cgridloc(loop)==2)write(txtline2,'(a)')', i_index_u[d=1]'
         if(cgridloc(loop)==3)write(txtline2,'(a)')', i_index_v[d=1]'
        endif
       endif               !1111111111>

       if(axisnum==2) then !2222222222>
        if(axisunit==1) then
         if(cgridloc(loop)==0)write(txtline2,'(a)')', j_index_t[d=1]'
         if(cgridloc(loop)==1)write(txtline2,'(a)')', j_index_t[d=1]'
         if(cgridloc(loop)==2)write(txtline2,'(a)')', j_index_u[d=1]'
         if(cgridloc(loop)==3)write(txtline2,'(a)')', j_index_v[d=1]'
        endif
       endif               !2222222222>

        if(axisunit==2) then
         if(cgridloc(loop)==0)write(txtline2,'(a)')', latitude_t[d=1]'
         if(cgridloc(loop)==1)write(txtline2,'(a)')', latitude_t[d=1]'
         if(cgridloc(loop)==2)write(txtline2,'(a)')', latitude_u[d=1]'
         if(cgridloc(loop)==3)write(txtline2,'(a)')', latitude_v[d=1]'
        endif

        if(axisunit==3) then
         if(cgridloc(loop)==0)write(txtline2,'(a)')', longitude_t[d=1]'
         if(cgridloc(loop)==1)write(txtline2,'(a)')', longitude_t[d=1]'
         if(cgridloc(loop)==2)write(txtline2,'(a)')', longitude_u[d=1]'
         if(cgridloc(loop)==3)write(txtline2,'(a)')', longitude_v[d=1]'
        endif

       txtline=trim(txtline)//trim(txtline2)

       if(cgridloc(loop)==0)write(txtline,'(a,a)')trim(txtline) &
                                         ,', depth_w[d=1]'
       if(cgridloc(loop)==1)write(txtline,'(a,a)')trim(txtline) &
                                         ,', depth_t[d=1]'
       if(cgridloc(loop)==2)write(txtline,'(a,a)')trim(txtline) &
                                         ,', depth_u[d=1]'
       if(cgridloc(loop)==3)write(txtline,'(a,a)')trim(txtline) &
                                         ,', depth_v[d=1]'

       else                !0000000000>

         if(cgridloc(loop)==0)write(txtline2,'(a)') &
         ', longitude_t[d=1], latitude_t[d=1]'
         if(cgridloc(loop)==1)write(txtline2,'(a)') &
         ', longitude_t[d=1], latitude_t[d=1]'
         if(cgridloc(loop)==2)write(txtline2,'(a)') &
         ', longitude_u[d=1], latitude_u[d=1]'
         if(cgridloc(loop)==3)write(txtline2,'(a)') &
         ', longitude_v[d=1], latitude_v[d=1]'

          txtline=trim(txtline)//trim(txtline2)
       endif               !0000000000>

       write(10,'(a)')trim(txtline)

       write(txtkount,'(i0)')kount
       txtline='frame/file="'//trim(varname(loop)) &
                             //trim(txtkount)//'.gif"'

       write(10,'(a)')trim(txtline)

       write(10,'(a)')'cancel d 2'

      enddo             !---------->

      goto 201 

  120 close(10)


      end
