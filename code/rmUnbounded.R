
files.names = c("el1sparse.csv", "res_loqo_el1sparse.csv", "el1.csv", "res_loqo_el1.csv",
 "l1ls_res.csv", "mirrorprox_res.csv", "fhtp_res.csv", "fhtp_res_agnostic.csv")



#show the number of unbounded and remove these columns
for(i in 1:length(file.names)){
  dat = read.csv(file.names[i], header=FALSE)

  idx = which(is.na(dat[,4]))
  if(length(idx)>0){
    write(paste0(file.names[i]," contains ",length(idx)," instances of unbounded problems out of ",
     nrow(dat)," instances."), stdout())
    dat = as.matrix(dat[-idx,])
    dat = matrix(as.numeric(dat), nrow=nrow(dat), ncol=ncol(dat))
    write.table(dat,file=file.names[i],quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE)
  }
}


