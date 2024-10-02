#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

rdir <-as.character(args[1])
var <- as.character(args[2])
ano <- as.character(as.double(args[3]))
sdir <-as.character(args[4])

if (require(RNetCDF)==FALSE) install.packages("RNetCDF")
#if (require(abind)==FALSE) install.packages("abind")
library(RNetCDF)
#library(abind)

#rdir <-"C:/Datos/Res_C32_Belich/runs_3D/RUN_reference"
#var <- "zoonanoc"
#ano<-"2016" 
#sdir <-"C:/Datos/Res_C32_Belich/runs_3D/RUN_reference"

  #archivo<-list.files(rdir,pattern="ECO3M-S.LA.tlse.nc")
  archivo<-intersect(list.files(rdir,pattern="ECO3M-S.LA.tlse.nc"),list.files(rdir,pattern=ano))
  eliminar<-list.files(rdir,pattern=".bkp")
  archivo<-setdiff(archivo, eliminar)
  fecha<-substring(archivo,1,8)
  hora<-substring(archivo,10,15)
  nombre<-unique(substring(archivo,1,15))
  tiempo<-rep(NA,length(archivo))
  arrayEd<-array(data=NA, dim=c(432,532,20,length(archivo)),dimnames=NULL)
  
  for (k in 1:length(archivo)){
      archivoEd <-archivo[k]
      filename <- paste(rdir,archivoEd,sep="/")
      filenc <- open.nc(filename)
      #print.nc(filenc)
      filerc <- read.nc(filenc)
      tiempo[k]<-filerc$time
      Ed<-filerc[[which(names(filerc)==var)]]
      arrayEd[,,,k]<-Ed
      print(archivoEd)
  } # end loop k
  save(tiempo,arrayEd,file=paste(sdir,"/",var,"_",ano,"_3D.RData",sep=""))
   