if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
library(RNetCDF)


## Get atmospheric data from ERA5
#################################
#Hourly data on single levels (total column): tclw, tcwv, tco3
#https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
# extract directly from the web, through the Download tab
# OR
# do an API request (python code) with the to_download_ERA5 env activated

rdir<-paste0(OS,"Datos/Res_C37_Bottom_RRS/setup_MM/reanalysis-era5-single-levels/")
sdir<-paste0(OS,"Datos/Res_C37_Bottom_RRS/setup_MM/")

archivos<-list.files(rdir, pattern=".nc",recursive = TRUE)

## save era5_singleLevels.dat file
sink(paste(sdir,"era5_singleLevels.dat",sep=""))

for (k in c(1:length(archivos))){
  archivo<-archivos[k]
  filename <- paste(rdir,archivo, sep="")
    filenc <- open.nc(filename)
    #print.nc(filenc)
    filerc <- read.nc(filenc)
    #names(filerc)

    latitude<-filerc$latitude
    longitude<-filerc$longitude
    fechas2<-as.POSIXlt(filerc$valid_time, origin="1969-12-31 23:00:01")
    FECH<-substring(as.character(fechas2),1,19)
    tclw<-apply(filerc$tclw, MARGIN=3,FUN=mean, na.rm=T)
    tcwv<-apply(filerc$tcwv, MARGIN=3,FUN=mean, na.rm=T)
    tco3<-apply(filerc$tco3, MARGIN=3,FUN=mean, na.rm=T)

      datoss <- cbind(FECH,tclw,tco3,tcwv)
      #head(datoss)
      # era5_singleLevels.dat
      write.table(datoss, file=paste(sdir,"era5_singleLevels.dat",sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
}
sink()
