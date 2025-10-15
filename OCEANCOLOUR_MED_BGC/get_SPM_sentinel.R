#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

site <- as.character(args[1])
name <- as.character(args[2])
lat <- as.character(args[3])
lon <- as.character(args[4])

site <-paste0(site,"/sentinel_hr/")
#name <- "B"
#lati <- "37.7"
#long <- "-0.8"

#if (require(RNetCDF)==FALSE) install.packages("RNetCDF")
#if (require(abind)==FALSE) install.packages("abind")
#if (require(pracma)==FALSE) install.packages("pracma")
library(RNetCDF)
#library(abind)
library(pracma)
library(chron)
library(lubridate)

## GET DATA
###########
    archive<-paste(site,"cmems_obs_oc_med_bgc_tur-spm-chl_nrt_l3-hr-mosaic_P1D-m.nc", sep="")
    filenc <- open.nc(archive)
    tmp <- read.nc(filenc)
    names(tmp)
    
    lon <- tmp$longitude
    lat <- tmp$latitude
    tim <- tmp$time
    
    chl <- tmp$CHL
    spm <- tmp$SPM
    tur <- tmp$TUR
    dimnames(chl)=list(x=lon, y=lat, t=tim)
    dimnames(spm)=list(x=lon, y=lat, t=tim)
    dimnames(tur)=list(x=lon, y=lat, t=tim)

## EXTRACT VALUES AT LOCATIONS (mean & sd)
##############################
    lim_norte<-as.numeric(lati)+0.01
    lim_sur  <-as.numeric(lati)-0.01
    lim_oeste<-as.numeric(long)-0.01
    lim_este <-as.numeric(long)+0.01
    
    # CHL
    datos_chl2<-rep(NA,length=length(tim))
    datos_chl3<-rep(NA,length=length(tim))
    for (k in c(1:length(tim))){
      station <- chl[lon>lim_oeste & lon<lim_este,  lat>lim_sur & lat<lim_norte, k]
      # dim(station)
      if (sum(!is.na(station))>2) {datos_chl2[k]<-mean(station, na.rm=T)
                                   datos_chl3[k]<-sd(station, na.rm=T)}
    }
    fechas2<-chron(c(0:(length(tim)-1)),origin=c(month=1, day=1, year=2020), out.format=c("y-m-d","h:m:s"))
    julianos<- julian(x=month(fechas2), d=day(fechas2), y=year(fechas2), origin.=c(month = 1, day = 1, year = 1998))        
    fechas2<-fechas2[!is.na(datos_chl2)]
    julianos<-julianos[!is.na(datos_chl2)]
    datos_chl3<-datos_chl3[!is.na(datos_chl2)]
    datos_chl2<-datos_chl2[!is.na(datos_chl2)]
    fechas2<-as.Date(fechas2)

    # SPM
    datos_spm2<-rep(NA,length=length(tim))
    datos_spm3<-rep(NA,length=length(tim))
    for (k in c(1:length(tim))){
      station <- spm[lon>lim_oeste & lon<lim_este,  lat>lim_sur & lat<lim_norte, k]
      # dim(station)
      if (sum(!is.na(station))>2) {datos_spm2[k]<-mean(station, na.rm=T)
                                   datos_spm3[k]<-sd(station, na.rm=T)}
    }

    # TUR
    datos_tur2<-rep(NA,length=length(tim))
    datos_tur3<-rep(NA,length=length(tim))
    for (k in c(1:length(tim))){
      station <- tur[lon>lim_oeste & lon<lim_este,  lat>lim_sur & lat<lim_norte, k]
      # dim(station)
      if (sum(!is.na(station))>2) {datos_tur2[k]<-mean(station, na.rm=T)
                                   datos_tur3[k]<-sd(station, na.rm=T)}
    }


## write .obs file
##################
      FECH<-paste(paste(year(fechas2),sprintf("%02d",month(fechas2)),sprintf("%02d",day(fechas2)), sep="-")," 00:00:00", sep="")      
      VALO<-datos_chl2
      SVAL<-datos_spm2
      TVAL<-datos_tur2
      FECH<-FECH[!is.na(VALO)]
      STDV<-STDV[!is.na(VALO)]
      VALO<-VALO[!is.na(VALO)]
      VALO[VALO<0]<-0
      
      sink(paste(site,name,"_Chl_SPM_TUR_mean.obs",sep=""))
      datoss  <- cbind(paste("",format(FECH, width = 10, justify = "left"), sep=""),
                       paste("",format(VALO, width = 10, justify = "right"), sep=""),
		       paste("",format(SVAL, width = 10, justify = "right"), sep=""),
		       paste("",format(TVAL, width = 10, justify = "right"), sep="")) 
      write.table(datoss, file=paste(site,name,"_Chl_SPM_TUR_mean.obs",sep=""),sep = "   ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      sink()