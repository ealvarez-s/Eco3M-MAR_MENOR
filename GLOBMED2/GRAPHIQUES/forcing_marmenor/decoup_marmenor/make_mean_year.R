#!/usr/bin/env Rscript
##########################################################
## Crear año promedio GLOBMED2 con 2014, 2015, 2017 y 2018
##########################################################
args <- commandArgs(trailingOnly = TRUE)

year <- as.character(args[1])
sdir <-as.character(args[2])

if (require(RNetCDF)==FALSE) install.packages("RNetCDF")
if (require(stringr)==FALSE) install.packages("stringr")
library(RNetCDF)
library(stringr)

rdir <-getwd()

# daily files Eco3M 2014
o_dir<-paste0(rdir,"/GLOBMED2_2014")
listica<-list.files(path=paste(o_dir, sep=""),recursive = T, full.names = T)

# colocar ano medio en nuevo_ano
nuevo_ano=as.numeric(year)
o_dir<-paste0(sdir,"/",nuevo_ano,"")
listica2<-list.files(path=paste(o_dir, sep=""),recursive = T, full.names = T)

## Para cada dia, encontrar el dia en todos los años
for (i in 1:length(listica)){
  fecha<-substring(listica[i],nchar(listica[i])-34,nchar(listica[i])-31)
  # find the rest
  filenames <-list.files(path=paste(rdir, sep=""), recursive = T, pattern=paste0(fecha,"_"), full.names = T)
  filenc <- open.nc(filenames[1])
  #print.nc(filenc)  
  filerc <- read.nc(filenc)

        ## Para cada variable
        for (k in c(1:34,36:38)){
        nombre<-names(filerc)[k]
        dimensiones<- dim(filerc[[k]])
        ARRAY<-array(NA,dim=c(dimensiones,length(filenames)))
        
        ## Hacer la media
        for (j in c(1:length(filenames))){
        filenc <- open.nc(filenames[j])
        #print.nc(filenc)  
        filerc <- read.nc(filenc)
        dim(filerc[[k]])
        ARRAY[,,,j]<-filerc[[k]]}
        datos<-apply(ARRAY, MARGIN=c(1,2,3), FUN=mean, na.rm=T)
        
        # Rellenar la variable (usar como template el primer año[2014])
        filenc <- open.nc(listica2[i],write=T)
        #print.nc(filenc)  
        var.put.nc(ncfile=filenc, variable=nombre,  data=datos, start=c(1,1,1,1), count=c(dimensiones,1))
        close.nc(filenc)
        } # end loop k variable
        
  # Cambiar el nombre
  new_name<-str_replace(listica2[i],"/2014",paste0("/",nuevo_ano))
  file.rename(listica2[i],new_name)
  filenc <- open.nc(new_name,write=T)
  
  # Cambiar el tiempo (origin="2013-08-15 00:00:00")
  diferencia<-nuevo_ano-2014
  var.put.nc(ncfile=filenc, variable="time", data=(filerc$time+(diferencia*31622400)))  
  close.nc(filenc)

} # end loop i, cada dia