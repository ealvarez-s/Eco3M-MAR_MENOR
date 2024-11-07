#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

rdir <-as.character(args[1])
var <- as.character(args[2])
sdir <-as.character(args[3])

#rdir <-"/home/alvarez/EVA/Eco3M-MAR_MENOR/GLOBMED2/GRAPHIQUES/forcing_marmenor/xtr"
#var <- "silice"
#sdir <-"/home/alvarez/EVA/Eco3M-MAR_MENOR/GLOBMED2/GRAPHIQUES/forcing_marmenor"

if (require(RNetCDF)==FALSE) install.packages("RNetCDF")
#if (require(abind)==FALSE) install.packages("abind")
library(RNetCDF)
#library(abind)

  archivo<-list.files(rdir,pattern="ECO3M-S.LA.tlse_xtrW.nc")
  eliminar<-list.files(rdir,pattern=".gz")
  archivo<-setdiff(archivo, eliminar)
  fecha<-substring(archivo,1,8)
  hora<-substring(archivo,10,15)
  nombre<-unique(substring(archivo,1,15))
  tiempo<-rep(NA,length(archivo))
  arrayEd<-array(data=NA, dim=c(26,28,60,length(archivo)),dimnames=NULL)

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
  save(tiempo,arrayEd,file=paste(sdir,"/",var,"_BC_3D.RData",sep=""))