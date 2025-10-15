#!/usr/bin/env Rscript

#####################################################################
## Sustituir variables Eco3M por variables MedBFM4 en el año promedio
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

year <- as.character(args[1])
index <- as.character(args[2])

if (require(RNetCDF)==FALSE) install.packages("RNetCDF")
if (require(stringr)==FALSE) install.packages("stringr")
if (require(akima)==FALSE) install.packages("akima")
if (require(chron)==FALSE) install.packages("chron")
library(RNetCDF)
library(stringr)
library(akima)
library(chron)

## Read the Eco3M grid
o_dir<-"/home/alvarez/EVA/Eco3M-MAR_MENOR/GLOBMED2/GRAPHIQUES/forcing_marmenor/decoup_marmenor/"
    filename <- paste(o_dir,"grid_xtrW.nc", sep="")
    filenc <- open.nc(filename)
    filerc <- read.nc(filenc)    
    longitude_Eco3M<-filerc$longitude_t
    latitude_Eco3M<-filerc$latitude_t
    depth_Eco3M<-filerc$depth_t
    mask_Eco3M<-filerc$mask_t
    
    
## Variables and units in the MedBFM4 reanalysis    
o_dir<-"/home/alvarez/EVA/Eco3M-MAR_MENOR/MedBFM4/"
    variable<-c("o2", "dissic","talk", "no3","nh4","po4", "chl","phyc")
    datasets<-c("bio",   "car", "car", "nut","nut","nut", "pft", "pft")
    VAR_Eco3M<-c("oxygen", "dic","talk", "nitrate","ammonium","phosphate", "synechl","synec")
    unidades<-c(1,1000,1000,1,1,1,1,1)
    
    # Read one dataset to get dates, latitude, lontitude and depth
    filename <- paste(o_dir,"analysis_forecast/",datasets[1],"_",year,".nc", sep="")
    filenc <- open.nc(filename)
    filerc <- read.nc(filenc)
    tim <- filerc$time
    #### Origen del CMS request -> 01-01 de year  
    fechas2<-chron(c(0:(length(tim)-1)),origin=c(month=1, day=1, year=as.numeric(year)), out.format=c("y-m-d","h:m:s"))
    fechas3<-as.Date(fechas2)
    longitude_MedBFM4 <- filerc$longitude
    latitude_MedBFM4 <- filerc$latitude
    depth_MedBFM4 <- filerc$depth


## Files Eco3M to modify # copied from the mean year
o_dir<-paste0(getwd(),"/",year)
filenames <-list.files(path=paste(o_dir, sep=""),pattern="ECO3M-S.LA.tlse_xtrW.nc",recursive = T,full.names = T)

#     for (q in c(1:length(filenames))){
    for (q in c(as.numeric(index):length(filenames))){
      fechaCHAR<-substring(filenames[q],nchar(filenames[q])-38,nchar(filenames[q])-31)
      YYYY<-substring(fechaCHAR,1,4)
      MM<-substring(fechaCHAR,5,6)
      DD<-substring(fechaCHAR,7,8)
      f<-as.Date(paste(YYYY,MM,DD, sep="-"))
      # que fecha CMS corresponde al archivo Eco3M
      cual<-which(fechas3==f) # la primera (2014-01-01 deberia ser el 366)
      #dim(filerc$o2[,,,cual])
      print(f) # Eco3M
      print(fechas3[cual]) # MedBFM
      
      # read daily file Eco3M
      filename <- filenames[q]
      filenc <- open.nc(filename)
      fileEco3M <- read.nc(filenc)
      #names(fileEco3M)
      
      # vectorizar el grid en el que esta MedBFM4
      longitud<-rep(longitude_MedBFM4,times=length(latitude_MedBFM4))
      latitud<-rep(latitude_MedBFM4,each=length(longitude_MedBFM4))
      profundidad<-rep(depth_MedBFM4,each=length(longitud))
      longitud<-rep(longitud,times=length(depth_MedBFM4))
      latitud<-rep(latitud,times=length(depth_MedBFM4))
      #length(longitud)
      #length(latitud)
      #length(profundidad)
      
      
      for (v in c(1:length(variable))){
        VAR<-variable[v]
        o_dir<-"/home/alvarez/EVA/Eco3M-MAR_MENOR/MedBFM4/"
        # open dataset that contains the variable
        filename <- paste(o_dir,"analysis_forecast/",datasets[v],"_",year,".nc", sep="")
        filenc <- open.nc(filename)
        filerc <- read.nc(filenc)
        # read variable
        cu<-which(names(filerc)==VAR)
        lista<-c(filerc[[cu]][,,,cual])
        valor<-filerc[[cu]][,,,cual]
        #dim(valor)
        #range(valor, na.rm=T)
        #length(lista)
        
        # vectorizar el grid en el que lo quiero (Eco3M)
        largo<-c(longitude_Eco3M)
        ancho<-c(latitude_Eco3M)
        fondo<-(-(c(depth_Eco3M)))
        #length(fondo)
        #range(fondo)
        #depth_Eco3M[1,1,] # los indeces estan al reves
        largo<-rep(largo,times=dim(depth_Eco3M)[3])
        ancho<-rep(ancho,times=dim(depth_Eco3M)[3])     
        #length(largo)
        #length(ancho)
        #dim(fileEco3M$oxygen)
        #range(fileEco3M$oxygen)
        cu<-which(names(fileEco3M)==VAR_Eco3M[v])
        file_Eco3M<-fileEco3M[[cu]]
        valor_Eco3M<-c(fileEco3M[[cu]])
        valor_Eco3M[valor_Eco3M==0]<-NA
     
        
     ## interpolar valores MedBFM4 (valor) al grid Eco3M
        new_depth<-array(data=NA,dim=c(dim(valor)[c(1,2)],dim(depth_Eco3M)[c(3)]))
        #dim(new_depth)
        
        for (i in c(1:dim(valor)[1])){
          for (j in c(1:dim(valor)[2])){
            resta<-longitude_Eco3M-longitude_MedBFM4[i]
            indexes1<-which(resta==min(resta,na.rm=T), arr.ind=T) 
            longitude_Eco3M[indexes1]
            resta<-latitude_Eco3M-latitude_MedBFM4[j]
            indexes2<-which(resta==min(resta,na.rm=T), arr.ind=T) 
            latitude_Eco3M[indexes2]
            if (sum(!is.na(valor[i,j,]))>2){
              column<-approx(x=-depth_MedBFM4, y=valor[i,j,], xout=depth_Eco3M[indexes1[1],indexes2[2],], method = "linear", rule = 2, f = 0, ties = mean)
              new_depth[i,j,]<-column$y} else {new_depth[i,j,]<-rep(NA,length=dim(depth_Eco3M)[c(3)])}
          }} # end loop interpolate only depth
        
        
        # nuevo grid con coordenadas
        new_coordinates<-array(data=NA,dim=dim(depth_Eco3M))
        #dim(new_coordinates)
        longitud<-rep(longitude_MedBFM4,times=length(latitude_MedBFM4))
        latitud<-rep(latitude_MedBFM4,each=length(longitude_MedBFM4))
        largo<-c(longitude_Eco3M)
        ancho<-c(latitude_Eco3M)   
        
        for (h in c(1:dim(depth_Eco3M)[3])){
          valorcito<-c(new_depth[,,h])
          longitud2<-longitud[!is.na(valorcito)]
          latitud2<-latitud[!is.na(valorcito)]
          valorcito<-valorcito[!is.na(valorcito)]
          g<-interp(x=longitud2,y=latitud2,z=valorcito,xo=largo,yo=ancho)
          #resta1<-g$x-longitude_Eco3M
          #resta2<-g$y-latitude_Eco3M
          #indexes<-which(resta1==min(resta1,na.rm=T) & resta2==min(resta2,na.rm=T))
          #new_coordinates[,,h]<-matrix(g$z[indexes],nrow=dim(depth_Eco3M)[1],ncol=dim(depth_Eco3M)[2])
          for (i in c(1:dim(depth_Eco3M)[1])){
            for (j in c(1:dim(depth_Eco3M)[2])){
              index1<-match(longitude_Eco3M[i,j],largo)
              index2<-match(latitude_Eco3M[i,j],ancho)
              new_coordinates[i,j,h]<-g$z[index1,index2]
            }}
        } # end loop interpolate coordinates
        
        #dim(mask_Eco3M)
        new_coordinates[mask_Eco3M==0]<-NA
        #dim(new_coordinates)
        new_coordinates<-new_coordinates*unidades[v]
        assign(x=VAR,value=new_coordinates)
        print(VAR)
      } # end loop v variable
      
      # divide state variables as initialization
      #chl
      diachl<-chl*0.5
      nanochl<-chl*0.35
      synechl<-chl*0.15
      #phyc          
      diac<-diachl*6
      nanoc<-nanochl*6
      synec<-synechl*6
      # n,p,si
      synen<-synec*((0.05+0.2)/2)
      synep<-synec*((0.004+0.019)/2)
      nanon<-nanoc*((0.05+0.2)/2)
      nanop<-nanoc*((0.002+0.019)/2)
      dian<-diac*((0.05+0.2)/2)
      diap<-diac*((0.002+0.019)/2)
      diasi<-diac*((0.05+0.19)/2)   
      
      # modify daily .nc file Eco3M
      filename <- filenames[q]
      filenc <- open.nc(filename,write=T)
      #print.nc(filenc)
      var.put.nc(ncfile=filenc, variable="diachl",  data=diachl, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="nanochl", data=nanochl, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="synechl", data=synechl, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="diac",  data=diac, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="nanoc", data=nanoc, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="synec", data=synec, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="dian",  data=dian, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="nanon", data=nanon, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="synen", data=synen, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="diap",  data=diap, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="nanop", data=nanop, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="synep", data=synep, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="diasi", data=diasi, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="oxygen", data=o2, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="dic", data=dissic, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="talk", data=talk, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="nitrate", data=no3, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="ammonium", data=nh4, start=c(1,1,1,1), count=c(26,28,60,1))
      var.put.nc(ncfile=filenc, variable="phosphate", data=po4, start=c(1,1,1,1), count=c(26,28,60,1))
      close.nc(filenc)
      print(filename)
    } # end loop q fechas