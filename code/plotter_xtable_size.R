library(xtable)

alldat = list(0)

alldat[[1]] = read.csv("el1sparse_size.csv",header=FALSE) 
alldat[[2]] = read.csv("res_loqo_el1sparse_size.csv",header=FALSE)

alldat[[3]] = read.csv("el1_size.csv",header=FALSE)
alldat[[4]] = read.csv("res_loqo_el1_size.csv",header=FALSE)

alldat[[5]] = read.csv("l1ls_res_size.csv",header=FALSE)
alldat[[6]] = read.csv("mirrorprox_res_size.csv",header=FALSE)
alldat[[7]] = read.csv("fhtp_res_size.csv",header=FALSE)

for(i in 1:7){
  val = sort(unique(alldat[[i]][,1]),decreasing=FALSE)
  tmp = matrix(NA,ncol=ncol(alldat[[i]]),nrow=length(val))
  for(j in 1:length(val)){
    idx = which(alldat[[i]][,1] == val[j])
    tmp[j,] = apply(alldat[[i]][idx,],2,mean)
  }
  alldat[[i]] = tmp
}

#construct xtables
levels = c(102,202,302,402)

for(i in c(1,3)){
  k = 0
  tmp = matrix(NA,ncol=3,nrow=7)
  for(j in 1:7){
    idx = which(alldat[[j]][,1]==levels[i+k])
	if(length(idx)>0) tmp[j,] = alldat[[j]][idx,2:4] 
  }
  
  colnames(tmp) = c("Time (Seconds)","Relative Error", "Uniform Error")
  rownames(tmp) = c("Simplex KCS","IPM KCS","Simplex","IPM","l1ls", "Mirror Prox", "FHTP (Oracle)")
  
  k = 1
  tmp2 = matrix(NA,ncol=3,nrow=7)
  for(j in 1:7){
    idx = which(alldat[[j]][,1]==levels[i+k])
	if(length(idx)>0) tmp2[j,] = alldat[[j]][idx,2:4] 
  }
  
  colnames(tmp2) = c("Time (Seconds)","Relative Error", "Uniform Error")
  rownames(tmp2) = c("Simplex KCS","IPM KCS","Simplex","IPM","l1ls", "Mirror Prox", "FHTP (Oracle)")
  
  tmp = cbind(tmp,tmp2)
  print(xtable(tmp,floating=FALSE,digits=5))
}
