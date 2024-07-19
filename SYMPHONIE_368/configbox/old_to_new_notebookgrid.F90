      program old_to_new_notebook_grid
      implicit none
      double precision dxb,dyb,phi0,longi0,latit0,rayonterre,pole_lon,pole_lat
      integer i0,j0,typegrid,i1d,iglb,jglb,kmax
      character*1 texte

      open(unit=10,file='notebook_grid')
      read(10,*)dxb
      read(10,*)dyb
      read(10,*)phi0
      read(10,*)longi0
      read(10,*)latit0
      read(10,*)
      read(10,*)i0
      read(10,*)j0
      read(10,*)rayonterre
      read(10,*)typegrid
      read(10,*)i1d
      read(10,*)iglb,jglb,kmax
      read(10,*)pole_lon
      read(10,*)pole_lat
      close(10)


      open(unit=10,file='notebook_grid_s26')
      write(10,'(a)')'&notebook_grid'
      write(10,*)'iglb=',iglb
      write(10,*)'jglb=',jglb
      write(10,*)'kmax=',kmax
      write(10,*)'dxb=',dxb
      write(10,*)'dyb=',dyb
      write(10,*)'phi0=',phi0
      write(10,*)'longi0=',longi0
      write(10,*)'latit0=',latit0
      write(10,*)'grid_i0=',i0
      write(10,*)'grid_j0=',j0
      write(10,*)'rayonterre=',rayonterre
      write(10,*)'northpole_lon=',pole_lon
      write(10,*)'northpole_lat=',pole_lat
      write(10,'(a)')'southpole_lon=9999.'
      write(10,'(a)')'southpole_lat=9999.'
      write(10,'(4a)')'lonlatfile=',achar(39),'nofile',achar(39)
      write(10,'(a)')'nbdom_imax=1'
      write(10,'(a)')'nbdom_jmax=1'
      write(10,'(a)')'iperiodicboundary=.false.'
      write(10,'(a)')'jperiodicboundary=.false.'
      write(10,'(4a)')'initialgridfile_txt=',achar(39),'none',achar(39)
      write(10,'(4a)')'vert_axis_conv_direc=',achar(39),'dw',achar(39)
      write(10,'(4a)')'vert_axis_conv_start=',achar(39),'w',achar(39)
      write(10,'(4a)')'vert_axis_conv_end=',achar(39),'t',achar(39)
      write(10,'(4a)')'hori_axis_conv_start=',achar(39),'t',achar(39)

      write(10,'(a)')'/'
      close(10)


      end
