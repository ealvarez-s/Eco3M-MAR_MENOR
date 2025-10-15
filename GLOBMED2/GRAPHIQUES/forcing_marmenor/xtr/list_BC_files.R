#!/usr/bin/env Rscript

rdir <-getwd()
archivo<-list.files(rdir,pattern="ECO3M-S.LA.tlse_xtrW.nc",full.names = T,recursive = TRUE)
eliminar<-list.files(rdir,pattern=".zip",full.names = T, recursive = TRUE)
archivos<-setdiff(archivo, eliminar)

# SAVE biovarfilelist 
sink(paste(rdir,"/biovarfilelist",sep=""))
# no HEADER
write.table(archivos, file=paste(rdir,"/biovarfilelist",sep=""),
            sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
sink()