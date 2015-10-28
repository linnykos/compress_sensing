library(xtable)

alldat = list(0)

files.names = c("el1sparse.csv", "res_loqo_el1sparse.csv", "el1.csv", "res_loqo_el1.csv",
 "l1ls_res.csv", "mirrorprox_res.csv", "fhtp_res.csv", "fhtp_res_agnostic.csv")

for(i in 1:8){
  alldat[[i]] = read.csv(file.names[i], header=FALSE)
}

for(i in 1:8){
  val = sort(unique(alldat[[i]][,1]),decreasing=FALSE)
  tmp = matrix(NA,ncol=ncol(alldat[[i]]),nrow=length(val))
  for(j in 1:length(val)){
    idx = which(alldat[[i]][,1] == val[j])
    tmp[j,] = apply(alldat[[i]][idx,],2,median)
  }
  alldat[[i]] = tmp
}

#construct xtables
levels = c(2,20,50,70,100,150)

for(i in c(1,3,5)){
  k = 0
  tmp = matrix(NA,ncol=3,nrow=8)
  for(j in 1:8){
    tmp[j,] = alldat[[j]][which(alldat[[j]][,1]==levels[i+k]),2:4]
  }
  
  colnames(tmp) = c("Time (Seconds)","Relative Error", "Uniform Error")
  rownames(tmp) = c("Simplex KCS","IPM KCS","Simplex","IPM","l1ls", "Mirror Prox", "FHTP (Oracle)", "FHTP (Agnostic)")
  
  k = 1
  tmp2 = matrix(NA,ncol=3,nrow=8)
  for(j in 1:8){
    tmp2[j,] = alldat[[j]][which(alldat[[j]][,1]==levels[i+k]),2:4]
  }
  
  colnames(tmp2) = c("Time (Seconds)","Relative Error", "Uniform Error")
  rownames(tmp2) = c("Simplex KCS","IPM KCS","Simplex","IPM","l1ls", "Mirror Prox", "FHTP (Oracle)", "FHTP (Agnostic)")
  
  tmp = cbind(tmp,tmp2)
  print(xtable(tmp,floating=FALSE,digits=5))
}
