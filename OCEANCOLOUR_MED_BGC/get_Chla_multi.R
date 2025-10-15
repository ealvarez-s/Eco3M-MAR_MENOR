#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

site <- as.character(args[1])
name <- as.character(args[2])
lat <- as.character(args[3])
lon <- as.character(args[4])

site <-paste0(site,"/multi/")
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

## GET RRS DATA
###############
    archive<-paste(site,"cmems_obs-oc_med_bgc-reflectance_my_l3-multi-1km_P1D.nc", sep="")
    filenc <- open.nc(archive)
    tmp <- read.nc(filenc)
    names(tmp)
    lon <- tmp$longitude
    lat <- tmp$latitude
    tim <- tmp$time
    rrs555 <- tmp$RRS555
    rrs670 <- tmp$RRS670

## COMPUTE CHL-a
################
RG=log10(rrs670/rrs555)
Log10Chla=0.353*RG^3+2.132*RG^2+3.905*RG+2.110   #fluorescence
#Log10Chla=0.965*RG^3+4.456*RG^2+6.200*RG+2.619  #spectrophot.
rrs=10^Log10Chla
dimnames(rrs)=list(x=lon, y=lat, t=tim)

## EXTRACT VALUES AT LOCATIONS (mean & sd)
##############################
    lim_norte<-as.numeric(lati)+0.01
    lim_sur  <-as.numeric(lati)-0.01
    lim_oeste<-as.numeric(long)-0.01
    lim_este <-as.numeric(long)+0.01

    datos2<-rep(NA,length=length(tim))
    datos3<-rep(NA,length=length(tim))
    for (k in c(1:length(tim))){
      station <- rrs[lon>lim_oeste & lon<lim_este,  lat>lim_sur & lat<lim_norte, k]
      # dim(station)
      if (sum(!is.na(station))>2) {datos2[k]<-mean(station, na.rm=T)
                                   datos3[k]<-sd(station, na.rm=T)}
    }
    fechas2<-chron(c(0:(length(tim)-1)),origin=c(month=1, day=1, year=1998), out.format=c("y-m-d","h:m:s"))
    julianos<- julian(x=month(fechas2), d=day(fechas2), y=year(fechas2), origin.=c(month = 1, day = 1, year = 1998))        
    fechas2<-fechas2[!is.na(datos2)]
    julianos<-julianos[!is.na(datos2)]
    datos3<-datos3[!is.na(datos2)]
    datos2<-datos2[!is.na(datos2)]
    fechas2<-as.Date(fechas2)

## SMOOTH TIMESERIES
####################
    ## Running averages
    ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
    datos2_rmean <- ma(datos2,n=30) #1month
    ## Savitzky-Golay smoothing
    datos2_savgol <-savgol(datos2, fl=91)

## write .obs file
##################
      FECH<-paste(paste(year(fechas2),sprintf("%02d",month(fechas2)),sprintf("%02d",day(fechas2)), sep="-")," 00:00:00", sep="")      
      VALO<-datos2
      STDV<-datos3
      FECH<-FECH[!is.na(VALO)]
      STDV<-STDV[!is.na(VALO)]
      VALO<-VALO[!is.na(VALO)]
      #data.frame(FECH,VALO)
      VALO[VALO<0]<-0
      
      sink(paste(site,name,"_Chla_mean_sd.obs",sep=""))
      datoss  <- cbind(paste("",format(FECH, width = 10, justify = "left"), sep=""),
                       paste("",format(VALO, width = 10, justify = "right"), sep=""),
		       paste("",format(STDV, width = 10, justify = "right"), sep="")) 
      write.table(datoss, file=paste(site,name,"_Chla_mean_sd.obs",sep=""),sep = "   ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      sink()