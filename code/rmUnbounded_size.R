file.names = c("el1sparse_size.csv", "res_loqo_el1sparse_size.csv", "el1_size.csv",
 "res_loqo_el1_size.csv", "l1ls_res_size.csv", "mirrorprox_res_size.csv", "fhtp_res_size.csv")


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

