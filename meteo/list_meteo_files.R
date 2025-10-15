#!/usr/bin/env Rscript

rdir <-getwd()
archivo<-list.files(rdir,pattern="previ-Med",full.names = T, recursive = T)
eliminar<-list.files(rdir,pattern=".zip",full.names = T, recursive = T)
archivos<-setdiff(archivo, eliminar)

# SAVE liste_ecmwf_dp2m
sink(paste(rdir,"/liste_ecmwf_dp2m",sep=""))
# no HEADER
write.table(archivos, file=paste(rdir,"/liste_ecmwf_dp2m",sep=""),
            sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
sink()